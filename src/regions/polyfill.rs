use crate::base_cells::{_is_base_cell_pentagon, baseCellNumToCell};
use crate::bbox::{bbox_contains_bbox, bbox_overlaps_bbox, bbox_to_cell_boundary, bboxes_from_geo_polygon};
use crate::constants::{
  CELL_SCALE_FACTOR, CHILD_SCALE_FACTOR, EPSILON, EPSILON_RAD, MAX_CELL_BNDRY_VERTS, MAX_EDGE_LENGTH_RADS, MAX_H3_RES,
  M_PI, M_PI_2, NORTH_POLE_CELLS, NUM_BASE_CELLS, RES0_BBOXES, SOUTH_POLE_CELLS, VALID_RANGE_BBOX,
};
use crate::coords::face_ijk::{Overage, _face_ijk_pent_to_cell_boundary, _face_ijk_to_cell_boundary}; // Overage might not be needed here directly
use crate::h3_index::inspection::is_valid_cell as h3_is_valid_cell; // Alias for clarity
use crate::h3_index::{
  _h3_to_face_ijk, get_base_cell, get_index_digit, get_mode, get_resolution, is_pentagon, set_index_digit,
  set_resolution,
};
use crate::hierarchy::parent_child::{cell_to_center_child, cell_to_children_size};
use crate::iterators::{iterInitParent, iterStepChild, IterCellsChildren};
use crate::polygon::{
  cell_boundary_crosses_polygon, cell_boundary_inside_polygon, flag_get_containment_mode, point_inside_polygon,
  validate_polygon_flags,
};
use crate::types::{BBox, CellBoundary, ContainmentMode, GeoPolygon, H3Error, H3Index, LatLng, H3_NULL};
use std::ptr;

/// Internal helper: Given a cell, find the next cell in the sequence of all cells
/// to check in the iteration (used by IterCellsPolygonCompact).
fn next_cell_for_polyfill(mut cell: H3Index) -> H3Index {
  let mut res = get_resolution(cell);
  loop {
    // If this is a base cell, set to next base cell (or H3_NULL if done)
    if res == 0 {
      let next_bc = get_base_cell(cell) + 1;
      return if next_bc < (NUM_BASE_CELLS as i32) {
        baseCellNumToCell(next_bc) // Assumes this returns H3_NULL if next_bc is invalid
      } else {
        H3_NULL
      };
    }

    // Faster cellToParent when we know the resolution is valid
    // and we're only moving up one level
    let mut parent = cell;
    set_resolution(&mut parent, res - 1);
    // In H3, digits beyond a cell's resolution are set to 7 (InvalidDigit).
    // When finding a parent, we effectively set the current res digit to 7.
    set_index_digit(&mut parent, res, crate::types::Direction::InvalidDigit);

    // If not the last sibling of parent, return next sibling
    let digit = get_index_digit(cell, res);
    if digit < crate::types::Direction::InvalidDigit {
      // Max valid child digit is IJ (6)
      let mut next_digit_val = digit as u8 + 1;
      // Skip K-axis for pentagon children if current parent is a pentagon
      if is_pentagon(parent) && next_digit_val == crate::types::Direction::KAxes as u8 {
        next_digit_val += 1;
      }

      if next_digit_val < crate::types::Direction::InvalidDigit as u8 {
        set_index_digit(
          &mut cell,
          res,
          crate::types::Direction::try_from(next_digit_val).unwrap_or(crate::types::Direction::InvalidDigit),
        );
        return cell;
      }
    }
    // Move up to the parent for the next loop iteration
    res -= 1;
    cell = parent;
  }
}

/// Iterator for traversing H3 cells within a polygon, producing a compact set.
#[derive(Debug)]
pub struct IterCellsPolygonCompact {
  /// Current H3 cell in the iteration. `H3_NULL` if exhausted or error.
  pub cell: H3Index,
  /// Error code, if any occurred during initialization or iteration.
  pub error: H3Error,
  _res: i32,                       // Target resolution
  _flags: u32,                     // Mode flags
  _polygon_ptr: *const GeoPolygon, // Raw pointer to avoid lifetime issues with _bboxes
  _bboxes_ptr: *mut BBox,          // Raw pointer to heap-allocated bboxes
  _num_bboxes: usize,              // Number of bboxes (1 outer + N holes)
  _started: bool,                  // Whether iteration has started
}

impl IterCellsPolygonCompact {
  /// Internal initializer without the first step.
  fn _new(polygon: &GeoPolygon, res: i32, flags: u32) -> Self {
    let mut iter = Self {
      cell: baseCellNumToCell(0), // Start with first base cell
      error: H3Error::Success,
      _res: res,
      _flags: flags,
      _polygon_ptr: polygon as *const GeoPolygon,
      _bboxes_ptr: ptr::null_mut(),
      _num_bboxes: 0,
      _started: false,
    };

    if res < 0 || res > MAX_H3_RES {
      iter.error = H3Error::ResDomain;
      iter.cell = H3_NULL;
      return iter;
    }

    if let Err(e) = validate_polygon_flags(flags) {
      iter.error = e;
      iter.cell = H3_NULL;
      return iter;
    }

    iter._num_bboxes = polygon.num_holes + 1;
    // Use H3_MEMORY from C for exact port, or Rust's global allocator.
    // For now, let's use Rust's allocator.
    let bboxes_layout = std::alloc::Layout::array::<BBox>(iter._num_bboxes).unwrap();
    // Safety: Layout is valid.
    iter._bboxes_ptr = unsafe { std::alloc::alloc(bboxes_layout) as *mut BBox };

    if iter._bboxes_ptr.is_null() {
      iter.error = H3Error::MemoryAlloc;
      iter.cell = H3_NULL;
      return iter;
    }
    // Safety: _bboxes_ptr is valid, _polygon_ptr is valid for lifetime of iter.
    unsafe {
      bboxes_from_geo_polygon(
        &*iter._polygon_ptr,
        std::slice::from_raw_parts_mut(iter._bboxes_ptr, iter._num_bboxes),
      );
    }
    iter
  }

  /// Steps the iterator to the next compact H3 cell.
  pub fn step(&mut self) {
    if self.cell == H3_NULL || self.error != H3Error::Success {
      self.destroy_internal_data(); // Ensure cleanup if error or already done
      self.cell = H3_NULL;
      return;
    }

    let mut current_cell_iter = self.cell;

    if self._started {
      current_cell_iter = next_cell_for_polyfill(current_cell_iter);
    } else {
      self._started = true;
    }

    // Safety: polygon_ptr and bboxes_ptr are assumed valid for the iter's lifetime.
    let polygon = unsafe { &*self._polygon_ptr };
    let bboxes = unsafe { std::slice::from_raw_parts(self._bboxes_ptr, self._num_bboxes) };

    if polygon.geoloop.num_verts == 0 && polygon.num_holes == 0 {
      self.destroy_internal_data();
      self.cell = H3_NULL;
      return;
    }

    let mode = match crate::polygon::flag_get_containment_mode(self._flags) {
      Ok(m) => m,
      Err(e) => {
        self.error = e;
        self.cell = H3_NULL;
        self.destroy_internal_data(); // Clean up
        return;
      }
    };

    while current_cell_iter != H3_NULL {
      let cell_res = get_resolution(current_cell_iter);

      if cell_res == self._res {
        // Target resolution
        let mut should_output = false;
        match mode {
          ContainmentMode::Center | ContainmentMode::Overlapping | ContainmentMode::OverlappingBbox => {
            let mut center = LatLng::default();
            match crate::indexing::cell_to_lat_lng(current_cell_iter) {
              Ok(center_val) => {
                // center_val is the LatLng
                if point_inside_polygon(polygon, bboxes, &center_val) {
                  should_output = true;
                }
              }
              Err(e) => {
                self.error = e; // cell_to_lat_lng failed
                self.cell = H3_NULL;
                return; // Exit step function on error
              }
            }
          }
          _ => {} // Handled below or not applicable for this first check
        }

        if !should_output
          && (mode == ContainmentMode::Full
            || mode == ContainmentMode::Overlapping
            || mode == ContainmentMode::OverlappingBbox)
        {
          match crate::indexing::cell_to_boundary(current_cell_iter) {
            Ok(boundary_val) => {
              // boundary_val is the CellBoundary
              let mut cell_bbox = BBox::default();
              match cell_to_bbox(current_cell_iter, &mut cell_bbox, false) {
                Ok(()) => {
                  if (mode == ContainmentMode::Full || mode == ContainmentMode::OverlappingBbox)
                    && cell_boundary_inside_polygon(polygon, bboxes, &boundary_val, &cell_bbox)
                  {
                    should_output = true;
                  } else if (mode == ContainmentMode::Overlapping || mode == ContainmentMode::OverlappingBbox)
                    && cell_boundary_crosses_polygon(polygon, bboxes, &boundary_val, &cell_bbox)
                  {
                    should_output = true;
                  }
                }
                Err(e) => {
                  self.error = e; // cell_to_bbox failed
                  self.cell = H3_NULL;
                  return; // Exit step function on error
                }
              }
            }
            Err(e) => {
              self.error = e; // cell_to_boundary failed
              self.cell = H3_NULL;
              return; // Exit step function on error
            }
          }
        }

        if !should_output && mode == ContainmentMode::OverlappingBbox {
          let mut cell_children_bbox = BBox::default(); // BBox to cover children
          if cell_to_bbox(current_cell_iter, &mut cell_children_bbox, true).is_ok() {
            if bbox_overlaps_bbox(&bboxes[0], &cell_children_bbox) {
              // Further checks if outer bbox overlaps...
              let bbox_boundary = bbox_to_cell_boundary(&cell_children_bbox);
              if bbox_contains_bbox(&cell_children_bbox, &bboxes[0]) || // cell bbox contains polygon
                               point_inside_polygon(polygon, bboxes, &bbox_boundary.verts[0]) || // polygon contains cell bbox
                               cell_boundary_crosses_polygon(polygon, bboxes, &bbox_boundary, &cell_children_bbox)
              {
                // polygon crosses cell bbox
                should_output = true;
              }
            }
          } else {
            self.error = H3Error::Failed;
            break;
          }
        }

        if should_output {
          self.cell = current_cell_iter;
          return;
        }
      } else if cell_res < self._res {
        // Coarser cell
        let mut children_bbox = BBox::default();
        if cell_to_bbox(current_cell_iter, &mut children_bbox, true).is_ok() {
          // coverChildren = true
          if bbox_overlaps_bbox(&bboxes[0], &children_bbox) {
            let children_bbox_boundary = bbox_to_cell_boundary(&children_bbox);
            if cell_boundary_inside_polygon(polygon, bboxes, &children_bbox_boundary, &children_bbox) {
              self.cell = current_cell_iter; // Output this coarse cell
              return;
            }
            // Intersects, so recurse: start with its first child.
            // cellToCenterChild then _zero_index_digits to get first child conceptually
            let center_child = cell_to_center_child(current_cell_iter, cell_res + 1).unwrap_or(H3_NULL); // Error case unlikely if parent valid
            if center_child == H3_NULL {
              self.error = H3Error::Failed;
              break;
            } // Should not happen
            current_cell_iter = center_child;
            continue; // Restart loop with the child cell
          }
        } else {
          self.error = H3Error::Failed;
          break;
        } /* cellToBBox error */
      }
      // If cell_res > self._res, or if bbox didn't overlap, or if fine-grained check failed:
      current_cell_iter = next_cell_for_polyfill(current_cell_iter);
    }

    // If loop finishes, we're done or an error occurred
    self.destroy_internal_data();
    self.cell = H3_NULL; // Mark as exhausted
  }

  /// Destroys internal data. Safe to call multiple times.
  fn destroy_internal_data(&mut self) {
    if !self._bboxes_ptr.is_null() {
      // Safety: _bboxes_ptr was allocated by us with corresponding layout.
      let layout = std::alloc::Layout::array::<BBox>(self._num_bboxes).unwrap();
      unsafe { std::alloc::dealloc(self._bboxes_ptr as *mut u8, layout) };
      self._bboxes_ptr = ptr::null_mut();
    }
  }
}

impl Drop for IterCellsPolygonCompact {
  fn drop(&mut self) {
    self.destroy_internal_data();
  }
}

/// Iterator for traversing all H3 cells within a polygon at a specific resolution.
#[derive(Debug)]
pub struct IterCellsPolygon {
  /// Current H3 cell in the iteration. `H3_NULL` if exhausted or error.
  pub cell: H3Index,
  /// Error code, if any occurred.
  pub error: H3Error,
  _cell_iter: IterCellsPolygonCompact, // Sub-iterator for compact cells
  _child_iter: IterCellsChildren,      // Sub-iterator for cell children
}

impl IterCellsPolygon {
  /// Initializes the iterator.
  pub fn new(polygon: &GeoPolygon, res: i32, flags: u32) -> Self {
    let cell_iter = IterCellsPolygonCompact::_new(polygon, res, flags); // Use internal init
    let mut child_iter = iterInitParent(cell_iter.cell, res); // iterInitParent handles H3_NULL for cell_iter.cell

    let initial_cell;
    if cell_iter.error == H3Error::Success {
      if cell_iter.cell == H3_NULL {
        // Compact iterator exhausted on init
        initial_cell = H3_NULL;
      } else if get_resolution(cell_iter.cell) == res {
        initial_cell = cell_iter.cell; // First compact cell is at target res
        child_iter.h = H3_NULL; // Mark child_iter as "done" for this compact cell
      } else {
        initial_cell = child_iter.h; // First child of first compact cell
      }
    } else {
      initial_cell = H3_NULL; // Error during compact_iter init
    }

    Self {
      cell: initial_cell,
      error: cell_iter.error,
      _cell_iter: cell_iter,
      _child_iter: child_iter,
    }
  }

  /// Steps the iterator to the next H3 cell.
  pub fn step(&mut self) {
    if self.cell == H3_NULL || self.error != H3Error::Success {
      self.destroy_internal_data();
      self.cell = H3_NULL;
      return;
    }

    // Try to step the child iterator first
    if self._child_iter.h != H3_NULL {
      iterStepChild(&mut self._child_iter);
      if self._child_iter.h != H3_NULL {
        self.cell = self._child_iter.h;
        return;
      }
    }

    // Child iterator exhausted, step the compact cell iterator
    self._cell_iter.step();
    if self._cell_iter.cell != H3_NULL {
      // Initialize new child iterator for the new compact cell
      self._child_iter = iterInitParent(self._cell_iter.cell, self._cell_iter._res);

      // If the compact cell itself is at target res, it's the next value.
      // Otherwise, the first child from the new child_iter is.
      if get_resolution(self._cell_iter.cell) == self._cell_iter._res {
        self.cell = self._cell_iter.cell;
        self._child_iter.h = H3_NULL; // Mark child_iter as "done" for this one
      } else {
        self.cell = self._child_iter.h;
      }
    } else {
      // Compact cell iterator exhausted or errored
      self.cell = H3_NULL;
      self.error = self._cell_iter.error; // Propagate error
      self.destroy_internal_data();
    }
  }

  /// Destroys internal data.
  fn destroy_internal_data(&mut self) {
    self._cell_iter.destroy_internal_data();
    // Child iterator doesn't own heap data directly in current Rust port
  }
}

impl Drop for IterCellsPolygon {
  fn drop(&mut self) {
    self.destroy_internal_data();
  }
}

// Public API functions
pub fn max_polygon_to_cells_size(polygon: &GeoPolygon, res: i32, flags: u32) -> Result<i64, H3Error> {
  // Similar logic to C's H3_EXPORT(maxPolygonToCellsSizeExperimental)
  if polygon.geoloop.num_verts == 0 && polygon.num_holes == 0 {
    return Ok(0);
  }

  let mut iter = IterCellsPolygonCompact::_new(polygon, res, flags);
  if iter.error != H3Error::Success {
    let err = iter.error;
    iter.destroy_internal_data(); // Explicitly destroy on error path
    return Err(err);
  }

  let mut count: i64 = 0;
  iter.step(); // Initial step
  while iter.cell != H3_NULL {
    if iter.error != H3Error::Success {
      let err = iter.error;
      iter.destroy_internal_data();
      return Err(err);
    }
    let children_count = cell_to_children_size(iter.cell, res)?;
    count = count.saturating_add(children_count);
    iter.step();
  }
  // Check for error one last time after loop finishes
  if iter.error != H3Error::Success {
    let err = iter.error;
    iter.destroy_internal_data();
    return Err(err);
  }
  iter.destroy_internal_data(); // Normal destruction
  Ok(count)
}

pub fn polygon_to_cells(polygon: &GeoPolygon, res: i32, flags: u32, out: &mut [H3Index]) -> Result<(), H3Error> {
  let mut iter = IterCellsPolygon::new(polygon, res, flags);
  if iter.error != H3Error::Success {
    return Err(iter.error);
  }

  let mut i: usize = 0;
  while iter.cell != H3_NULL {
    if iter.error != H3Error::Success {
      return Err(iter.error);
    }
    if i >= out.len() {
      iter.destroy_internal_data();
      return Err(H3Error::MemoryBounds);
    }
    out[i] = iter.cell;
    i += 1;
    iter.step();
  }
  // Check for error one last time after loop finishes
  if iter.error != H3Error::Success {
    return Err(iter.error);
  }

  // Fill remaining with H3_NULL if out buffer was larger than needed.
  // This depends on whether max_polygon_to_cells_size is exact or an overestimate.
  // If it's an overestimate, this is needed.
  // for j in i..out.len() {
  //     out[j] = H3_NULL;
  // }

  Ok(())
}

// Internal helper (should be in this module or a private utility module)
/// Calculates bounding box of a cell. If `cover_children` is true, expands the
/// bbox to guarantee coverage of all children at any finer resolution.
pub(crate) fn cell_to_bbox(cell: H3Index, out: &mut BBox, cover_children: bool) -> Result<(), H3Error> {
  let res = get_resolution(cell);

  if res == 0 {
    let base_cell = get_base_cell(cell);
    if base_cell < 0 || base_cell >= (NUM_BASE_CELLS as i32) {
      return Err(H3Error::CellInvalid);
    }
    *out = RES0_BBOXES[base_cell as usize];
  } else {
    let center = crate::indexing::cell_to_lat_lng(cell)?; // Use public API
    let edge_len_rad = MAX_EDGE_LENGTH_RADS[res as usize]; // From constants

    // A simple bbox approximation around the center
    // More accurate would be to get actual cell boundary and find its min/max lat/lng.
    // For performance, H3 often uses approximations for candidate searches.
    let lat_span = edge_len_rad; // Simplified
    let lng_span_at_equator = edge_len_rad;
    let lng_span = if center.lat.cos().abs() > EPSILON {
      lng_span_at_equator / center.lat.cos().abs()
    } else {
      // At or near a pole
      M_PI // Effectively covers all longitudes
    };

    out.north = center.lat + lat_span;
    out.south = center.lat - lat_span;
    out.east = center.lng + lng_span;
    out.west = center.lng - lng_span;
  }

  // Normalize and clamp
  out.north = out.north.min(M_PI_2);
  out.south = out.south.max(-M_PI_2);
  // Longitude will be handled by scaleBBox's internal normalization or bboxIsTransmeridian

  crate::bbox::scale_bbox(
    out,
    if cover_children {
      CHILD_SCALE_FACTOR
    } else {
      CELL_SCALE_FACTOR
    },
  );

  // Adjust for poles
  if cell == H3Index(NORTH_POLE_CELLS[res as usize]) {
    out.north = M_PI_2;
  }
  if cell == H3Index(SOUTH_POLE_CELLS[res as usize]) {
    out.south = -M_PI_2;
  }
  if out.north == M_PI_2 || out.south == -M_PI_2 {
    out.east = M_PI;
    out.west = -M_PI;
  }
  Ok(())
}

// Add tests for IterCellsPolygonCompact, IterCellsPolygon, max_polygon_to_cells_size, polygon_to_cells
// These tests would be similar to C's testPolygonToCells.c, testPolygonToCellsExperimental.c, etc.
// They require known polygons and expected H3 cell outputs for various resolutions and flags.
