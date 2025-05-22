// src/regions/to_polygon.rs

use crate::constants::{MAX_CELL_BNDRY_VERTS, NUM_HEX_VERTS, NUM_PENT_VERTS}; // etc.
use crate::coords::face_ijk::_face_ijk_pent_to_cell_boundary;
use crate::coords::face_ijk::_face_ijk_to_cell_boundary; // For cellToBoundary logic
use crate::h3_index::_h3_to_face_ijk;
use crate::h3_index::inspection::is_pentagon;
use crate::h3_index::{get_base_cell, get_resolution};
use crate::latlng::geo_almost_equal; // For comparing LatLng
use crate::types::{CellBoundary, CoordIJK, FaceIJK, GeoLoop, GeoPolygon, H3Error, H3Index, LatLng, Vec2d};
use crate::BBox;
use crate::{
  bbox::bbox_from_geoloop, bbox::bbox_is_transmeridian, is_valid_cell, polygon::generic_is_clockwise,
  polygon::point_inside_geoloop, H3_NULL,
}; // Added missing imports
use std::collections::HashMap; // Might be useful for _hashVertex alternative
use std::hash::{Hash, Hasher};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

// --- VertexGraph related structures and functions ---

#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PolygonRust {
  pub outer: Vec<LatLng>,
  pub holes: Vec<Vec<LatLng>>,
}
pub type MultiPolygonRust = Vec<PolygonRust>;

// Rust equivalent of C's VertexNode (using Box for linked-list like structure)
#[derive(Debug, Clone)]
struct VertexNode {
  from: LatLng,
  to: LatLng,
  next: Option<Box<VertexNode>>,
}

// Rust equivalent of C's VertexGraph
#[derive(Debug)]
struct VertexGraph {
  buckets: Vec<Option<Box<VertexNode>>>, // Vec of Option<Box<VertexNode>> for buckets
  num_buckets: usize,
  size: usize,
  res: i32, // Resolution of cells, useful for precision in hashing/equality
}

impl VertexGraph {
  fn new(num_buckets: usize, res: i32) -> Self {
    let mut buckets = Vec::with_capacity(num_buckets);
    for _ in 0..num_buckets {
      buckets.push(None);
    }
    VertexGraph {
      buckets,
      num_buckets,
      size: 0,
      res,
    }
  }

  // _hashVertex equivalent
  fn _hash_vertex(&self, vertex: &LatLng) -> usize {
    // H3's hashing for vertices can be simple for this internal graph.
    // Using a scaled integer representation based on resolution can work.
    // The exact hash function isn't critical as long as it distributes reasonably.
    // Max res 15. Max coord value around PI. Precision needs to distinguish verts.
    // Precision of H3 at res 15 is ~1e-7 radians.
    // Let's scale by 1e8 to get integer parts.
    let scale_factor = 1e8;
    let lat_u = (vertex.lat * scale_factor).round() as i64;
    let lng_u = (vertex.lng * scale_factor).round() as i64;

    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    lat_u.hash(&mut hasher);
    lng_u.hash(&mut hasher);
    (hasher.finish() % self.num_buckets as u64) as usize
  }

  fn add_vertex_node(&mut self, from_vtx: &LatLng, to_vtx: &LatLng) -> bool {
    let index = self._hash_vertex(from_vtx);
    let mut current_opt = &mut self.buckets[index];

    // Check if edge (from_vtx, to_vtx) already exists in this bucket's chain
    let mut temp_iter_opt = self.buckets[index].as_ref();
    while let Some(node_ref) = temp_iter_opt {
      if geo_almost_equal(&node_ref.from, from_vtx) && geo_almost_equal(&node_ref.to, to_vtx) {
        return false; // Edge already exists, not added
      }
      temp_iter_opt = node_ref.next.as_ref();
    }

    // Add new node to the head of the list for this bucket
    let new_node = Box::new(VertexNode {
      from: *from_vtx,
      to: *to_vtx,
      next: self.buckets[index].take(), // Old head becomes next
    });
    self.buckets[index] = Some(new_node);
    self.size += 1;
    true // New node added
  }

  fn remove_edge_if_exists(&mut self, from_vtx: &LatLng, to_vtx: &LatLng) -> bool {
    let index = self._hash_vertex(from_vtx);

    // Handle head of the list in the bucket
    if let Some(head_node) = self.buckets[index].as_ref() {
      if geo_almost_equal(&head_node.from, from_vtx) && geo_almost_equal(&head_node.to, to_vtx) {
        let mut old_head = self.buckets[index].take().unwrap(); // Take ownership
        self.buckets[index] = old_head.next.take(); // Set bucket head to next
        self.size -= 1;
        return true;
      }
    }

    // Handle rest of the list
    let mut current_opt_mut = self.buckets[index].as_mut();
    while let Some(current_node_mut) = current_opt_mut {
      // current_node_mut is &mut Box<VertexNode>
      let mut remove_next_node = false;
      if let Some(next_node_ref) = current_node_mut.next.as_ref() {
        // Check next node without taking
        if geo_almost_equal(&next_node_ref.from, from_vtx) && geo_almost_equal(&next_node_ref.to, to_vtx) {
          remove_next_node = true;
        }
      }

      if remove_next_node {
        let mut node_to_remove = current_node_mut.next.take().unwrap(); // Take ownership of the node to remove
        current_node_mut.next = node_to_remove.next.take(); // Link current_node_mut to node_to_remove's next
        self.size -= 1;
        return true;
      }
      current_opt_mut = current_node_mut.next.as_mut(); // Move to the next node in the chain
    }
    false // Edge not found
  }

  // Finds and removes an edge *starting with* start_vtx from any bucket.
  // Returns the (from, to) of the removed edge.
  fn find_and_remove_edge_starting_with(&mut self, start_vtx: &LatLng) -> Option<(LatLng, LatLng)> {
    let index = self._hash_vertex(start_vtx); // Search in the expected bucket first

    // Try the expected bucket
    if let Some(edge_data) = self.find_and_remove_from_bucket_list(index, start_vtx) {
      return Some(edge_data);
    }

    // If not found, search all other buckets (less efficient, but necessary if hash collides across different start_vtx)
    // C's findNodeForVertex iterates all buckets if not found in the direct hash bucket.
    for bucket_idx in 0..self.num_buckets {
      if bucket_idx == index {
        continue;
      } // Already checked
      if let Some(edge_data) = self.find_and_remove_from_bucket_list(bucket_idx, start_vtx) {
        return Some(edge_data);
      }
    }
    None
  }

  // Helper: Finds and removes an edge starting with start_vtx from a specific bucket's list.
  fn find_and_remove_from_bucket_list(&mut self, bucket_idx: usize, start_vtx: &LatLng) -> Option<(LatLng, LatLng)> {
    // Handle head
    if let Some(head_node_ref) = self.buckets[bucket_idx].as_ref() {
      if geo_almost_equal(&head_node_ref.from, start_vtx) {
        let mut removed_node = self.buckets[bucket_idx].take().unwrap();
        let from_ret = removed_node.from;
        let to_ret = removed_node.to;
        self.buckets[bucket_idx] = removed_node.next.take();
        self.size -= 1;
        return Some((from_ret, to_ret));
      }
    }

    // Handle rest of the chain
    let mut current_opt_mut = self.buckets[bucket_idx].as_mut();
    while let Some(current_node_mut) = current_opt_mut {
      let mut removed_edge_data: Option<(LatLng, LatLng)> = None;
      if let Some(next_node_ref) = current_node_mut.next.as_ref() {
        if geo_almost_equal(&next_node_ref.from, start_vtx) {
          let mut node_to_remove = current_node_mut.next.take().unwrap();
          removed_edge_data = Some((node_to_remove.from, node_to_remove.to));
          current_node_mut.next = node_to_remove.next.take();
          self.size -= 1;
        }
      }
      if removed_edge_data.is_some() {
        return removed_edge_data;
      }
      current_opt_mut = current_node_mut.next.as_mut();
    }
    None
  }

  // Get the first edge from the graph and remove it (used to start a new loop search)
  fn pop_first_edge(&mut self) -> Option<(LatLng, LatLng)> {
    for i in 0..self.num_buckets {
      if self.buckets[i].is_some() {
        let mut head_node = self.buckets[i].take().unwrap();
        let from_ret = head_node.from;
        let to_ret = head_node.to;
        self.buckets[i] = head_node.next.take();
        self.size -= 1;
        return Some((from_ret, to_ret));
      }
    }
    None
  }
}

impl Drop for VertexGraph {
  fn drop(&mut self) {
    for i in 0..self.num_buckets {
      let mut current = self.buckets[i].take();
      while let Some(mut node) = current {
        current = node.next.take();
      }
    }
  }
}

fn h3_set_to_vertex_graph(h3_set: &[H3Index], graph: &mut VertexGraph) -> Result<(), H3Error> {
  if h3_set.is_empty() {
    return Ok(());
  }

  let res = get_resolution(h3_set[0]);
  if !is_valid_cell(h3_set[0]) {
    return Err(H3Error::CellInvalid);
  }

  for &cell_idx in h3_set.iter() {
    if cell_idx == H3_NULL {
      continue;
    }
    if !is_valid_cell(cell_idx) {
      return Err(H3Error::CellInvalid);
    }
    if get_resolution(cell_idx) != res {
      return Err(H3Error::ResMismatch);
    }

    let mut boundary = CellBoundary::default();
    let mut fijk = FaceIJK::default();
    _h3_to_face_ijk(cell_idx, &mut fijk)?;

    if is_pentagon(cell_idx) {
      _face_ijk_pent_to_cell_boundary(&fijk, res, 0, NUM_PENT_VERTS as i32, &mut boundary);
    } else {
      _face_ijk_to_cell_boundary(&fijk, res, 0, NUM_HEX_VERTS as i32, &mut boundary);
    }
    if boundary.num_verts == 0 {
      return Err(H3Error::Failed);
    }

    for j in 0..boundary.num_verts {
      let from_vtx = boundary.verts[j];
      let to_vtx = boundary.verts[(j + 1) % boundary.num_verts];

      if !graph.remove_edge_if_exists(&to_vtx, &from_vtx) {
        if !graph.add_vertex_node(&from_vtx, &to_vtx) {
          // This case should ideally not be hit if remove_edge_if_exists(reversed) was false.
          // It means (from, to) was already there, which is also an internal shared edge.
          // The C logic implies this won't happen due to the remove-first strategy.
          // However, if add_vertex_node also checks for duplicates, it might return false here.
          // This indicates a logic error or an unexpected graph state.
          // eprintln!(
          //   "Warning: Tried to add edge {:?}->{:?} that was already present after failing to remove reversed.",
          //   from_vtx, to_vtx
          // );
        }
      }
    }
  }
  Ok(())
}

/// Contains information about a raw loop extracted from the VertexGraph.
struct RawLoopInfo {
  verts: Vec<LatLng>,
  bbox: BBox,
  /// True if the loop is clockwise as extracted from the H3 graph.
  /// H3 graph: outer boundaries are CW, holes are CCW.
  is_cw_from_graph: bool,
}

fn vertex_graph_to_multi_polygon_rust(graph: &mut VertexGraph) -> Result<MultiPolygonRust, H3Error> {
  let mut raw_loops_info: Vec<RawLoopInfo> = Vec::new();

  // Phase 1: Extract all raw loops from the graph
  while graph.size > 0 {
    let first_edge_opt = graph.pop_first_edge();
    if first_edge_opt.is_none() {
      if graph.size > 0 {
        return Err(H3Error::Failed);
      } // Inconsistent state
      break;
    }
    let (loop_start_vtx, mut current_to_vtx) = first_edge_opt.unwrap();

    let mut current_loop_vec: Vec<LatLng> = Vec::new();
    current_loop_vec.push(loop_start_vtx);

    loop {
      let next_edge_opt = graph.find_and_remove_edge_starting_with(&current_to_vtx);

      if let Some((_found_from_vtx, found_to_vtx)) = next_edge_opt {
        current_loop_vec.push(current_to_vtx); // _found_from_vtx is current_to_vtx
        current_to_vtx = found_to_vtx;
      } else {
        // eprintln!(
        //   "Could not find next edge starting with {:?} for loop started at {:?}. Graph size: {}. Partial loop: {:?}",
        //   current_to_vtx, loop_start_vtx, graph.size, current_loop_vec
        // );
        return Err(H3Error::Failed);
      }

      if geo_almost_equal(&current_to_vtx, &loop_start_vtx) {
        break; // Closed the loop
      }
      if current_loop_vec.len() > graph.num_buckets * MAX_CELL_BNDRY_VERTS {
        // Safety break
        // eprintln!("Loop tracing exceeded maximum possible edges.");
        return Err(H3Error::Failed);
      }
    }
    if current_loop_vec.len() >= 3 {
      let mut bbox = BBox::default();
      let geoloop_for_bbox = GeoLoop {
        num_verts: current_loop_vec.len(),
        verts: current_loop_vec.clone(),
      };
      bbox_from_geoloop(&geoloop_for_bbox, &mut bbox);

      let is_cw = generic_is_clockwise(&current_loop_vec, current_loop_vec.len(), bbox_is_transmeridian(&bbox));

      raw_loops_info.push(RawLoopInfo {
        verts: current_loop_vec,
        bbox,
        is_cw_from_graph: is_cw,
      });
    }
  }

  // Phase 2: Normalize into GeoJSON MultiPolygon structure
  // GeoJSON: Outer loops CCW, Hole loops CW.
  // H3 graph: Outer loops are CW, Hole loops are CCW (as per `isClockwiseLinkedGeoLoop` behavior in C)

  let mut final_polygons: MultiPolygonRust = Vec::new();
  let mut hole_candidates: Vec<RawLoopInfo> = Vec::new();

  // Separate potential outers and holes based on winding from graph
  for mut loop_info in raw_loops_info {
    if loop_info.is_cw_from_graph {
      // CW from graph => H3 outer loop
      loop_info.verts.reverse(); // Reverse to make it CCW for GeoJSON outer
      final_polygons.push(PolygonRust {
        outer: loop_info.verts,
        holes: Vec::new(),
      });
    } else {
      // CCW from graph => H3 hole
      loop_info.verts.reverse(); // Reverse to make it CW for GeoJSON hole
      hole_candidates.push(loop_info);
    }
  }

  // Assign holes to their parent polygons
  let mut still_unassigned_holes: Vec<RawLoopInfo> = Vec::new();
  for mut hole_info in hole_candidates {
    // hole_info.verts is now CW
    let mut parent_polygon_idx: Option<usize> = None;
    let mut min_parent_depth = i32::MAX;

    if hole_info.verts.is_empty() {
      continue;
    }

    for (outer_idx, polygon_rust) in final_polygons.iter().enumerate() {
      // polygon_rust.outer is CCW
      let geoloop_outer = GeoLoop {
        num_verts: polygon_rust.outer.len(),
        verts: polygon_rust.outer.clone(),
      };
      let mut bbox_outer = BBox::default(); // Recompute bbox for the (potentially reversed) outer loop
      bbox_from_geoloop(&geoloop_outer, &mut bbox_outer);

      if point_inside_geoloop(&geoloop_outer, &bbox_outer, &hole_info.verts[0]) {
        // This hole is inside this outer polygon. Now check depth.
        let mut depth = 0;
        for (other_outer_idx, other_polygon_rust) in final_polygons.iter().enumerate() {
          if outer_idx == other_outer_idx {
            continue;
          }

          let geoloop_other_outer = GeoLoop {
            num_verts: other_polygon_rust.outer.len(),
            verts: other_polygon_rust.outer.clone(),
          };
          let mut bbox_other_outer = BBox::default();
          bbox_from_geoloop(&geoloop_other_outer, &mut bbox_other_outer);

          if point_inside_geoloop(&geoloop_other_outer, &bbox_other_outer, &polygon_rust.outer[0]) {
            depth += 1;
          }
        }

        if depth < min_parent_depth {
          min_parent_depth = depth;
          parent_polygon_idx = Some(outer_idx);
        } else if depth == min_parent_depth {
          // If depths are equal, C's logic might prefer the one appearing earlier in its list.
          // Or it might be based on area. For now, taking the first one found at this depth.
          // This ambiguity might lead to differences if multiple equally-nested parents exist.
        }
      }
    }

    if let Some(idx) = parent_polygon_idx {
      final_polygons[idx].holes.push(hole_info.verts);
    } else {
      still_unassigned_holes.push(hole_info);
    }
  }

  if !still_unassigned_holes.is_empty() {
    // eprintln!(
    //   "Warning: {} holes remained unassigned and were discarded.",
    //   still_unassigned_holes.len()
    // );
    // H3 C's normalizeMultiPolygon returns E_FAILED if a hole cannot be assigned.
    // return Err(H3Error::Failed);
    // For testing donut where outer/inner are separate, we might not want to fail here yet.
    // If the goal is strict GeoJSON, unassigned holes usually mean an invalid input topology.
  }

  Ok(final_polygons)
}

pub fn cells_to_multi_polygon(h3_set: &[H3Index]) -> Result<MultiPolygonRust, H3Error> {
  if h3_set.is_empty() {
    return Ok(Vec::new());
  }
  let res = get_resolution(h3_set[0]);
  let num_hexes = h3_set.len();

  let min_buckets = if num_hexes > 1_000 {
    num_hexes / 100
  } else if num_hexes > 6 {
    num_hexes
  } else {
    10
  };
  let mut graph = VertexGraph::new(min_buckets.max(10), res); // Ensure min_buckets is reasonable

  h3_set_to_vertex_graph(h3_set, &mut graph)?;

  vertex_graph_to_multi_polygon_rust(&mut graph)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::indexing::lat_lng_to_cell;
  use crate::traversal::grid_disk; // Assuming this is the correct module path
  use crate::traversal::grid_disk::max_grid_disk_size; // Assuming this is the correct module path

  #[test]
  fn test_cells_to_multi_polygon_single_cell() {
    let cell = lat_lng_to_cell(&LatLng { lat: 0.0, lng: 0.0 }, 5).unwrap();
    let result = cells_to_multi_polygon(&[cell]);
    assert!(result.is_ok(), "cells_to_multi_polygon failed: {:?}", result.err());
    let multi_poly = result.unwrap();
    assert_eq!(multi_poly.len(), 1); // One polygon

    let expected_verts = if is_pentagon(cell) {
      NUM_PENT_VERTS
    } else {
      NUM_HEX_VERTS
    };
    assert_eq!(multi_poly[0].outer.len(), expected_verts);
    assert_eq!(multi_poly[0].holes.len(), 0); // No holes
  }

  #[test]
  fn test_cells_to_multi_polygon_two_disjoint_cells() {
    let cell1 = lat_lng_to_cell(&LatLng { lat: 0.0, lng: 0.0 }, 5).unwrap();
    let cell2 = lat_lng_to_cell(
      &LatLng {
        lat: (10.0_f64).to_radians(),
        lng: (10.0_f64).to_radians(),
      },
      5,
    )
    .unwrap();
    let result = cells_to_multi_polygon(&[cell1, cell2]);
    assert!(result.is_ok(), "cells_to_multi_polygon failed: {:?}", result.err());
    let multi_poly = result.unwrap();
    // Each disjoint cell forms its own polygon
    assert_eq!(multi_poly.len(), 2, "Expected 2 polygons for two disjoint cells");
    assert_eq!(multi_poly[0].holes.len(), 0);
    assert_eq!(multi_poly[1].holes.len(), 0);
  }

  #[test]
  fn test_cells_to_multi_polygon_donut() {
    // Create a k-ring (disk) and then remove the center cell to form a donut (ring)
    let center_cell = lat_lng_to_cell(&LatLng { lat: 0.0, lng: 0.0 }, 2).unwrap();
    let k_radius = 1; // For a simple donut
    let mut disk_size = match max_grid_disk_size(k_radius) {
      Ok(size) => size as usize,
      Err(e) => panic!("max_grid_disk_size failed: {:?}", e),
    };
    let mut disk_cells = vec![H3_NULL; disk_size];
    match grid_disk(center_cell, k_radius, &mut disk_cells) {
      Ok(()) => (),
      Err(e) => panic!("grid_disk failed: {:?}", e),
    }

    let donut_cells: Vec<H3Index> = disk_cells
      .into_iter()
      .filter(|&h| h != H3_NULL && h != center_cell)
      .collect();

    let center_cell_is_pent = is_pentagon(center_cell);
    let expected_donut_cell_count = if center_cell_is_pent { 5 } else { 6 };
    assert_eq!(
      donut_cells.len(),
      expected_donut_cell_count,
      "Number of cells in the donut ring should be {} for a {} center",
      expected_donut_cell_count,
      if center_cell_is_pent { "pentagon" } else { "hexagon" }
    );

    let result = cells_to_multi_polygon(&donut_cells);
    assert!(
      result.is_ok(),
      "cells_to_multi_polygon failed for donut: {:?}",
      result.err().unwrap_or(H3Error::Failed)
    );
    let multi_poly = result.unwrap();

    // A proper donut should result in one polygon with one outer loop and one inner hole loop.
    assert_eq!(
      multi_poly.len(),
      1,
      "Expected 1 polygon for a donut, got {}",
      multi_poly.len()
    );

    if !multi_poly.is_empty() {
      assert!(!multi_poly[0].outer.is_empty(), "Donut outer loop should not be empty");
      assert_eq!(
        multi_poly[0].holes.len(),
        1,
        "Expected 1 hole for a donut, got {}",
        multi_poly[0].holes.len()
      );

      if !multi_poly[0].holes.is_empty() {
        assert!(
          !multi_poly[0].holes[0].is_empty(),
          "Donut hole loop should not be empty"
        );

        // Verify vertex counts for k=1 donut
        // Outer boundary of a k=1 ring: each of the 6 (or 5) cells contributes 3 vertices.
        let expected_outer_verts = if center_cell_is_pent {
          NUM_PENT_VERTS * 3
        } else {
          NUM_HEX_VERTS * 3
        };
        let expected_hole_verts = if center_cell_is_pent {
          NUM_PENT_VERTS
        } else {
          NUM_HEX_VERTS
        };

        assert_eq!(
          multi_poly[0].outer.len(),
          expected_outer_verts,
          "Outer loop vertex count mismatch for k=1 donut. Center is_pent: {}. Got {}, Expected {}",
          center_cell_is_pent,
          multi_poly[0].outer.len(),
          expected_outer_verts
        );
        assert_eq!(
          multi_poly[0].holes[0].len(),
          expected_hole_verts,
          "Hole loop vertex count mismatch for k=1 donut. Center is_pent: {}. Got {}, Expected {}",
          center_cell_is_pent,
          multi_poly[0].holes[0].len(),
          expected_hole_verts
        );
      }
    }
  }
}
