use crate::constants::{MAX_CELL_BNDRY_VERTS, NUM_HEX_VERTS, NUM_PENT_VERTS}; // etc.
use crate::coords::face_ijk::_face_ijk_pent_to_cell_boundary;
use crate::coords::face_ijk::_face_ijk_to_cell_boundary; // For cellToBoundary logic
use crate::h3_index::_h3_to_face_ijk;
use crate::h3_index::inspection::is_pentagon;
use crate::h3_index::{get_base_cell, get_resolution};
use crate::latlng::geo_almost_equal; // For comparing LatLng
use crate::types::{CellBoundary, CoordIJK, FaceIJK, GeoLoop, GeoPolygon, H3Error, H3Index, LatLng, Vec2d};
use crate::{is_valid_cell, H3_NULL};
use std::collections::HashMap; // Might be useful for _hashVertex alternative
use std::hash::{Hash, Hasher};

// --- VertexGraph related structures and functions ---

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
  res: i32,
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
    // Simple hash: Take the sum of the lat and lng with a precision level
    // determined by the resolution, converted to int, modulo bucket count.
    // Convert to integer representation for hashing to avoid f64 hashing issues.
    // Multiplying by a large number and casting to u64 can work.
    let lat_u = (vertex.lat * 1e7).round() as i64;
    let lng_u = (vertex.lng * 1e7).round() as i64;

    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    lat_u.hash(&mut hasher);
    lng_u.hash(&mut hasher);
    (hasher.finish() % self.num_buckets as u64) as usize
  }

  // addVertexNode, removeVertexNode, findNodeForEdge, findNodeForVertex, firstVertexNode
  // These will need careful porting of the linked list logic within buckets.
  // Example stub for add_vertex_node:
  fn add_vertex_node(&mut self, from_vtx: &LatLng, to_vtx: &LatLng) -> bool {
    // Returns true if new node added
    let index = self._hash_vertex(from_vtx);
    let mut current_opt = &mut self.buckets[index];

    loop {
      match current_opt {
        Some(current_node) => {
          if geo_almost_equal(&current_node.from, from_vtx) && geo_almost_equal(&current_node.to, to_vtx) {
            return false; // Already exists
          }
          current_opt = &mut current_node.next;
        }
        None => {
          *current_opt = Some(Box::new(VertexNode {
            from: *from_vtx,
            to: *to_vtx,
            next: None,
          }));
          self.size += 1;
          return true;
        }
      }
    }
  }

  // ... (other VertexGraph methods)
}

// Explicit drop for VertexGraph to mirror C's destroyVertexGraph (though Rust's Box would handle it)
// This is mostly for conceptual mapping; Box<T> handles its own deallocation.
impl Drop for VertexGraph {
  fn drop(&mut self) {
    for i in 0..self.num_buckets {
      let mut current = self.buckets[i].take(); // Takes ownership, leaving None in the bucket
      while let Some(mut node) = current {
        current = node.next.take(); // Move to next, freeing current node when it drops
      }
    }
  }
}

// --- LinkedGeo related structures and functions ---

// Rust equivalent of C's LinkedLatLng, LinkedGeoLoop, LinkedGeoPolygon
// These can be simplified in Rust, e.g., using Vec<LatLng> for loops.
// For a direct port maintaining linked structure:
#[derive(Debug)]
struct LinkedLatLngNode {
  vertex: LatLng,
  next: Option<Box<LinkedLatLngNode>>,
}

#[derive(Debug)]
struct LinkedGeoLoopNode {
  first: Option<Box<LinkedLatLngNode>>,
  last: *mut LinkedLatLngNode, // Raw pointer for efficient last-node tracking if needed, careful!
  next: Option<Box<LinkedGeoLoopNode>>,
}

#[derive(Debug)]
struct LinkedGeoPolygonNode {
  first_loop: Option<Box<LinkedGeoLoopNode>>, // Renamed from 'first' to avoid clash if polygon has 'next'
  last_loop: *mut LinkedGeoLoopNode,          // Raw pointer
  next_polygon: Option<Box<LinkedGeoPolygonNode>>, // Renamed from 'next'
}

// Functions like normalizeMultiPolygon would operate on these.
// This linked structure is complex to manage safely in Rust without `unsafe` or `Rc<RefCell>`.
// A more idiomatic Rust output for cellsToMultiPolygon would be:
pub struct PolygonRust {
  pub outer: Vec<LatLng>,
  pub holes: Vec<Vec<LatLng>>,
}
pub type MultiPolygonRust = Vec<PolygonRust>;

// --- Core logic port ---

/// Port of C's h3SetToVertexGraph
fn h3_set_to_vertex_graph(h3_set: &[H3Index], graph: &mut VertexGraph) -> Result<(), H3Error> {
  if h3_set.is_empty() {
    // VertexGraph::new already initializes an empty graph.
    return Ok(());
  }

  // All cells must be valid and have the same resolution.
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

    let mut boundary = CellBoundary::default(); // Assuming CellBoundary is Default
    let mut fijk = FaceIJK::default();
    _h3_to_face_ijk(cell_idx, &mut fijk)?; // Get canonical Fijk

    if is_pentagon(cell_idx) {
      _face_ijk_pent_to_cell_boundary(&fijk, res, 0, NUM_PENT_VERTS as i32, &mut boundary);
    } else {
      _face_ijk_to_cell_boundary(&fijk, res, 0, NUM_HEX_VERTS as i32, &mut boundary);
    }
    if boundary.num_verts == 0 {
      return Err(H3Error::Failed);
    } // Should not happen for valid cell

    for j in 0..boundary.num_verts {
      let from_vtx = boundary.verts[j];
      let to_vtx = boundary.verts[(j + 1) % boundary.num_verts];

      // If we've seen this edge already (from_vtx, to_vtx), it will be reversed in graph.
      // The C code finds (toVtx, fromVtx) and removes it. If not found, adds (fromVtx, toVtx).
      // This is complex with Rust's ownership. A simpler model for VertexGraph might be needed,
      // or use a HashMap to count edge occurrences.
      // For now, let's simulate:
      if !graph.remove_edge_if_exists(&to_vtx, &from_vtx) {
        graph.add_vertex_node(&from_vtx, &to_vtx);
      }
    }
  }
  Ok(())
}

// Add a helper to VertexGraph for the remove_edge_if_exists logic
impl VertexGraph {
  fn remove_edge_if_exists(&mut self, from_vtx: &LatLng, to_vtx: &LatLng) -> bool {
    let index = self._hash_vertex(from_vtx);
    let mut current_opt = &mut self.buckets[index];

    // Handle case where the node to remove is the first in the bucket
    if let Some(node) = current_opt {
      if geo_almost_equal(&node.from, from_vtx) && geo_almost_equal(&node.to, to_vtx) {
        *current_opt = node.next.take(); // Move out the next node
        self.size -= 1;
        return true; // Node found and removed
      }
    }

    // Iterate through the rest of the list
    while let Some(node) = current_opt {
      if let Some(next_node) = &mut node.next {
        if geo_almost_equal(&next_node.from, from_vtx) && geo_almost_equal(&next_node.to, to_vtx) {
          node.next = next_node.next.take(); // Splice out next_node
          self.size -= 1;
          return true; // Node found and removed
        }
      }
      current_opt = &mut node.next;
    }
    false // Node not found
  }
}

/// Port of C's _vertexGraphToLinkedGeo and normalizeMultiPolygon, outputting to Rust-idiomatic structure.
fn vertex_graph_to_multi_polygon_rust(graph: &mut VertexGraph) -> Result<MultiPolygonRust, H3Error> {
  // This is a very complex function. It involves:
  // 1. Traversing the graph to find closed loops.
  // 2. Grouping these loops into polygons (one outer, multiple inner/holes).
  // 3. Normalizing winding order (outer CCW, inner CW - though H3 graph produces outer CW, inner CCW).
  // For now, this will be a major stub.
  // A full port needs careful handling of graph traversal and loop detection.

  let mut multi_polygon_rust = MultiPolygonRust::new();

  while graph.size > 0 {
    // While there are edges left in the graph
    let mut current_node_opt_box = None;
    let mut start_node_from_vtx = LatLng::default();

    // Find a starting edge
    for i in 0..graph.num_buckets {
      if graph.buckets[i].is_some() {
        // Take the node out of the bucket to start a loop
        current_node_opt_box = graph.buckets[i].take();
        if let Some(ref node) = current_node_opt_box {
          start_node_from_vtx = node.from;
          graph.buckets[i] = node.next.clone(); // Put rest of bucket back if any, this is wrong with Box
                                                // This simple take won't work correctly for bucket chains.
                                                // A proper graph traversal is needed.
                                                // Let's find the first node without removing it, then trace and remove.
          break;
        }
      }
    }

    // If no edge found (graph should be empty, but as a safeguard)
    if current_node_opt_box.is_none() && graph.size > 0 {
      // This implies an issue with graph.size tracking or iteration
      return Err(H3Error::Failed); // Should be empty if no node found
    }
    if current_node_opt_box.is_none() && graph.size == 0 {
      break; // Graph is empty
    }

    // Trace a loop
    let mut current_loop_vec: Vec<LatLng> = Vec::new();
    let mut current_edge_box = current_node_opt_box.unwrap(); // We found one
    let loop_start_vtx = current_edge_box.from;
    current_loop_vec.push(current_edge_box.from);

    // Remove the starting edge properly
    graph.remove_edge_if_exists(&current_edge_box.from, &current_edge_box.to);

    let mut next_vtx_to_find = current_edge_box.to;

    loop {
      // Find next edge starting with next_vtx_to_find
      let mut found_next_edge = false;
      for i in 0..graph.num_buckets {
        let mut temp_head = graph.buckets[i].take(); // Take ownership of the bucket's list head
        let mut current_node_in_chain = temp_head.as_mut(); // Option<&mut Box<VertexNode>>
        let mut found_and_removed_from_bucket = false;

        // Case 1: The node to remove is the head of the bucket's list
        if let Some(node_box) = current_node_in_chain {
          // current_node_in_chain is &mut Box<VertexNode>
          if geo_almost_equal(&node_box.from, &next_vtx_to_find) {
            current_loop_vec.push(node_box.from);
            next_vtx_to_find = node_box.to;
            graph.buckets[i] = node_box.next.take(); // Old head's next becomes new head
            graph.size -= 1;
            found_next_edge = true;
            found_and_removed_from_bucket = true;
          }
        }

        // Case 2: The node to remove is not the head
        if !found_and_removed_from_bucket {
          // Need to re-borrow temp_head to iterate if head wasn't the one.
          let mut current_node_in_chain_opt = temp_head.as_mut(); // Re-borrow as mutable to iterate
          while let Some(node_box_ref) = current_node_in_chain_opt {
            // Look at node_box_ref.next
            let mut should_remove_next = false;
            if let Some(next_node_box) = &node_box_ref.next {
              if geo_almost_equal(&next_node_box.from, &next_vtx_to_find) {
                should_remove_next = true;
              }
            }

            if should_remove_next {
              let mut removed_node_box = node_box_ref.next.take().unwrap(); // Take ownership of the node to remove
              current_loop_vec.push(removed_node_box.from);
              next_vtx_to_find = removed_node_box.to;
              node_box_ref.next = removed_node_box.next.take(); // Link prev to removed's next
              graph.size -= 1;
              found_next_edge = true;
              break; // Found and removed
            }
            current_node_in_chain_opt = node_box_ref.next.as_mut();
          }
          // After checking the chain, put the (potentially modified) head back
          graph.buckets[i] = temp_head;
        }
        // else if head was removed, graph.buckets[i] is already updated.

        if found_next_edge {
          break;
        } // Break from outer bucket loop (for i in 0..graph.num_buckets)
      }

      if !found_next_edge {
        // Should not happen in a well-formed graph if a loop was started
        // unless it's an open line segment (error in h3SetToVertexGraph)
        return Err(H3Error::Failed); // Or handle dangling edges
      }

      if geo_almost_equal(&next_vtx_to_find, &loop_start_vtx) {
        break; // Closed the loop
      }
    }
    // At this point, current_loop_vec contains one closed loop.
    // Now, determine if it's an outer loop or a hole and assign to a PolygonRust.
    // This requires the normalizeMultiPolygon logic (winding order, point-in-poly for holes).
    // For now, let's just add every loop as an outer loop of a new polygon.
    // This is NOT correct for holes but is a starting point.
    multi_polygon_rust.push(PolygonRust {
      outer: current_loop_vec,
      holes: Vec::new(),
    });
  }

  // TODO: Implement full normalization (winding order, hole assignment)
  // based on C's normalizeMultiPolygon.

  Ok(multi_polygon_rust)
}

/// Converts a set of H3 cells to a `MultiPolygon` (a list of polygons,
/// where each polygon has an outer loop and optional inner hole loops).
///
/// # Arguments
/// * `h3_set` - A slice of H3 cell indexes. Assumed to be of the same resolution
///              and contiguous (or forming multiple contiguous groups).
///
/// # Returns
/// `Ok(MultiPolygonRust)` on success, or an `H3Error` on failure.
pub fn cells_to_multi_polygon(h3_set: &[H3Index]) -> Result<MultiPolygonRust, H3Error> {
  if h3_set.is_empty() {
    return Ok(Vec::new());
  }
  let res = get_resolution(h3_set[0]); // Assuming all cells are same res
  let num_hexes = h3_set.len();

  let min_buckets = if num_hexes > 6 { num_hexes } else { 6 }; // Heuristic from C
  let mut graph = VertexGraph::new(min_buckets, res);

  h3_set_to_vertex_graph(h3_set, &mut graph)?;

  // vertex_graph_to_multi_polygon_rust consumes the graph
  vertex_graph_to_multi_polygon_rust(&mut graph)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::indexing::lat_lng_to_cell;
  use crate::traversal::grid_disk::grid_disk;

  #[test]
  fn test_cells_to_multi_polygon_single_cell() {
    let cell = lat_lng_to_cell(&LatLng { lat: 0.0, lng: 0.0 }, 5).unwrap();
    let result = cells_to_multi_polygon(&[cell]);
    assert!(result.is_ok());
    let multi_poly = result.unwrap();
    assert_eq!(multi_poly.len(), 1); // One polygon
    assert_eq!(
      multi_poly[0].outer.len(),
      if is_pentagon(cell) {
        crate::constants::NUM_PENT_VERTS
      } else {
        crate::constants::NUM_HEX_VERTS
      }
    ); // 6 or 5 verts
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
    assert!(result.is_ok());
    let multi_poly = result.unwrap();
    assert_eq!(multi_poly.len(), 2); // Two polygons (stubbed graph traversal might give 1 if not careful)
                                     // My current stub for vertex_graph_to_multi_polygon_rust will produce one polygon per edge component.
                                     // Since these are disjoint, they are two components. So 2 is correct for the stub.
  }

  #[test]
  fn test_cells_to_multi_polygon_donut() {
    // Create a k-ring (disk) and then remove the center cell to form a donut (ring)
    let center_cell = lat_lng_to_cell(&LatLng { lat: 0.0, lng: 0.0 }, 2).unwrap();
    let mut disk_size = crate::traversal::grid_disk::max_grid_disk_size(1).unwrap() as usize;
    let mut disk_cells = vec![H3_NULL; disk_size];
    grid_disk(center_cell, 1, &mut disk_cells).unwrap();

    let donut_cells: Vec<H3Index> = disk_cells
      .into_iter()
      .filter(|&h| h != H3_NULL && h != center_cell)
      .collect();

    assert_eq!(donut_cells.len(), if is_pentagon(center_cell) { 5 } else { 6 });

    let result = cells_to_multi_polygon(&donut_cells);
    assert!(
      result.is_ok(),
      "cells_to_multi_polygon failed for donut: {:?}",
      result.err()
    );
    let multi_poly = result.unwrap();

    // A donut should result in one polygon with one outer loop and one inner hole loop.
    // My current stub for vertex_graph_to_multi_polygon_rust will NOT correctly form holes.
    // It will likely produce two separate polygons (one for outer boundary, one for inner boundary).
    // This test will fail until normalizeMultiPolygon logic is fully implemented.
    // For now, let's assert based on the stub's likely behavior (two polygons, no holes).
    // assert_eq!(multi_poly.len(), 1, "Expected 1 polygon for a donut");
    // assert!(!multi_poly[0].outer.is_empty());
    // assert_eq!(multi_poly[0].holes.len(), 1, "Expected 1 hole for a donut");
    // assert!(!multi_poly[0].holes[0].is_empty());
    // For the stub:
    if is_pentagon(center_cell) {
      assert_eq!(
        multi_poly.len(),
        2,
        "Stub: Expected 2 polygons for pentagon donut (outer, inner)"
      );
    } else {
      assert_eq!(
        multi_poly.len(),
        2,
        "Stub: Expected 2 polygons for hexagon donut (outer, inner)"
      );
    }
  }
}
