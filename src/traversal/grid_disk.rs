// src/traversal/grid_disk.rs

use crate::constants::MAX_H3_RES;
use crate::coords::ijk::_rotate60_ccw;
use crate::h3_index::{get_resolution, is_pentagon}; // Assuming these are pub(crate) or pub
use crate::math::extensions::_ipow; // For K_ALL_CELLS_AT_RES_15 if needed, or direct const
use crate::traversal::neighbors::h3_neighbor_rotations; // From the module we just worked on
use crate::types::{Direction, H3Error, H3Index, H3_NULL}; // For allocating memory if fallback is needed

// K_ALL_CELLS_AT_RES_15 from C algos.c
// This constant defines the k-ring size that would encompass all cells at resolution 15.
// It's used to cap the size estimate from max_grid_disk_size to prevent overflow
// and to provide an upper bound related to the total number of cells.
const K_ALL_CELLS_AT_RES_15: i32 = 13_780_510;

/// Maximum number of H3 cells in a k-ring disk.
///
/// # Arguments
/// * `k` - The k-ring radius (k >= 0).
///
/// # Returns
/// `Ok(count)` of the maximum number of cells, or `H3Error::Domain` if k is negative.
pub fn max_grid_disk_size(k: i32) -> Result<i64, H3Error> {
  if k < 0 {
    return Err(H3Error::Domain);
  }
  // If k is very large, the formula 3*k*(k+1)+1 can overflow i64.
  // The C code caps this by K_ALL_CELLS_AT_RES_15.
  // If k is larger than or equal to a value that would contain all cells
  // at the finest resolution, we can return the total number of cells at res 15.
  // This check is slightly different from C's direct use of K_ALL_CELLS_AT_RES_15
  // in the formula, but aims for the same protective effect.
  // Maximum number of cells at res 15: 2 + 120 * 7^15
  // This value is 569_707_381_193_162.
  // The formula 3*k*(k+1)+1 for k=K_ALL_CELLS_AT_RES_15 is much larger than this.
  // So, if k could lead to an estimate larger than total cells at res 15, cap it.
  // A simpler approach matching C more directly:
  if k >= K_ALL_CELLS_AT_RES_15 {
    // crate::h3_index::get_num_cells(MAX_H3_RES) // This function doesn't exist yet, use const
    // NUM_CELLS_RES_15 needs to be defined in constants.rs
    // For now, let's use the C logic: 3*k*(k+1)+1.
    // If k is K_ALL_CELLS_AT_RES_15 or more, the C code still uses the formula,
    // but the actual gridDisk functions would likely hit cell limits or pentagons first.
    // The cap in C is more about preventing the formula itself from overflowing if k were enormous
    // and only relevant if getNumCells(MAX_H3_RES) was smaller than the formula's output.
    //
    // C code's maxGridDiskSize:
    // if (k >= K_ALL_CELLS_AT_RES_15) return H3_EXPORT(getNumCells)(MAX_H3_RES, out);
    // H3_EXPORT(getNumCells) needs to be ported. Let's make a const for now.
    // pub const NUM_CELLS_MAX_RES: i64 = 569_707_381_193_162; in constants.rs
    return Ok(crate::constants::NUM_CELLS_MAX_RES); // Assuming this constant exists
  }

  // Formula: 1 + sum(6*i) for i=1 to k = 1 + 6*k*(k+1)/2 = 1 + 3*k*(k+1)
  // Using i64 for k to prevent overflow during intermediate multiplication if k is large i32.
  let k_i64 = k as i64;
  Ok(3 * k_i64 * (k_i64 + 1) + 1)
}

/// Internal: BFS-based grid disk algorithm, handles pentagons correctly.
/// This is `_gridDiskDistancesInternal` from C algos.c.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `k` - The k-ring radius.
/// * `out` - Output array for H3 cells, pre-allocated and zeroed.
/// * `distances` - Output array for distances, pre-allocated and zeroed.
/// * `max_idx` - The allocated size of `out` and `distances`.
/// * `current_k` - The current distance from the origin in the recursion.
fn _grid_disk_distances_internal(
  origin: H3Index,
  k: i32,
  out: &mut [H3Index],
  distances: &mut [i32],
  max_idx: i64, // max_idx is usize from max_grid_disk_size's i64
  current_k: i32,
) -> Result<(), H3Error> {
  if origin == H3_NULL {
    // Should not happen with valid inputs
    return Ok(());
  }

  // Hash set like logic: find a slot for origin.
  // `max_idx` must be > 0 for this to be safe.
  // max_grid_disk_size ensures k=0 gives size 1.
  let mut offset = (origin.0 % max_idx as u64) as usize;
  loop {
    if out[offset] == H3_NULL {
      // Found an empty slot
      break;
    }
    if out[offset] == origin {
      // Already visited
      // If we found it via a shorter path earlier, we don't need to explore from here again.
      if distances[offset] <= current_k {
        return Ok(());
      }
      // Otherwise, we found a shorter path to this already processed cell, update distance.
      break;
    }
    // Collision, linear probe.
    offset = (offset + 1) % (max_idx as usize);
    // TODO: Add a cycle detection/limit here to prevent infinite loop on full table,
    // though max_grid_disk_size should be large enough.
  }

  out[offset] = origin;
  distances[offset] = current_k;

  if current_k >= k {
    return Ok(()); // Reached max depth
  }

  // Recurse to neighbors
  // Directions for iteration (K, J, JK, I, IK, IJ)
  // Using a fixed array to represent C's `DIRECTIONS` from algos.c if not already available
  const CELL_DIRECTIONS: [Direction; 6] = [
    Direction::KAxes,
    Direction::JAxes,
    Direction::JkAxes,
    Direction::IAxes,
    Direction::IkAxes,
    Direction::IjAxes,
  ];

  for dir_enum_val in 0..6 {
    let dir_to_neighbor = CELL_DIRECTIONS[dir_enum_val];

    // Skip K direction for pentagons if this logic matches C's use of DIRECTIONS
    // C's `DIRECTIONS` array starts with J_AXES (2) and goes up to IJ_AXES (6).
    // So C loop is for (d=0; d<6; d++) using `DIRECTIONS[d]`.
    // Our `CELL_DIRECTIONS` matches that order if K is 0, J is 1 etc.
    // No, `DIRECTIONS` in C's `_kRingInternal` (which is `_gridDiskDistancesInternal`) is:
    // `static const Direction DIRECTIONS[6] = {J_AXES_DIGIT, JK_AXES_DIGIT, K_AXES_DIGIT, IK_AXES_DIGIT, I_AXES_DIGIT, IJ_AXES_DIGIT};`
    // Let's use that exact order.
    const ACTUAL_C_DIRECTIONS: [Direction; 6] = [
      Direction::JAxes,
      Direction::JkAxes,
      Direction::KAxes,
      Direction::IkAxes,
      Direction::IAxes,
      Direction::IjAxes,
    ];

    let actual_dir = ACTUAL_C_DIRECTIONS[dir_enum_val];

    let mut rotations = 0;
    let mut next_neighbor = H3_NULL;

    let neighbor_result = h3_neighbor_rotations(origin, actual_dir, &mut rotations, &mut next_neighbor);

    // E_PENTAGON is an expected case when traversing off pentagons, meaning no neighbor in that dir.
    if neighbor_result.is_ok() {
      // Successfully found a neighbor
      _grid_disk_distances_internal(next_neighbor, k, out, distances, max_idx, current_k + 1)?;
    } else if matches!(neighbor_result, Err(H3Error::Pentagon)) {
      // It's a pentagon, and we tried to move into a deleted k-subsequence.
      // This is fine, just means no neighbor in that direction.
      continue;
    } else {
      // Any other error from h3_neighbor_rotations should be propagated.
      return neighbor_result;
    }
  }
  Ok(())
}

/// Produces H3 cells and their distances from the given origin cell, up to distance `k`.
/// Handles pentagons correctly.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `k` - The k-ring radius (k >= 0).
/// * `out_cells` - Output array for H3 cells. Must be pre-allocated to `max_grid_disk_size(k)` and zero-initialized.
/// * `out_distances` - Output array for distances. Must be pre-allocated like `out_cells` and zero-initialized.
///                     If `None`, distances will not be calculated or stored.
pub fn grid_disk_distances(
  origin: H3Index,
  k: i32,
  out_cells: &mut [H3Index],
  mut out_distances_opt: Option<&mut [i32]>,
) -> Result<(), H3Error> {
  let max_k_size = max_grid_disk_size(k)?;

  // Validate output buffer sizes
  if out_cells.len() < max_k_size as usize {
    return Err(H3Error::MemoryBounds);
  }
  if let Some(ref distances_slice) = out_distances_opt {
    if distances_slice.len() < max_k_size as usize {
      return Err(H3Error::MemoryBounds);
    }
  }

  // Zero out buffers (as C version expects)
  // Note: The caller is documented to provide zeroed arrays.
  // However, for safety within this Rust function, explicitly zeroing or careful handling is good.
  // If we rely on caller zeroing, these lines are not strictly needed here.
  // for item in out_cells.iter_mut() { *item = H3_NULL; }
  // if let Some(ref mut distances_slice) = out_distances_opt {
  //     for item in distances_slice.iter_mut() { *item = 0; }
  // }

  // The C version tries `gridDiskDistancesUnsafe` first. We can do that too.
  // For now, let's implement the safe version directly.
  // If `out_distances_opt` is None, we need a temporary buffer for internal use.
  if let Some(ref mut distances_slice) = out_distances_opt {
    _grid_disk_distances_internal(origin, k, out_cells, distances_slice, max_k_size, 0)
  } else {
    // Need to allocate a temporary distances array.
    // Requires an allocator. The C version does this using H3_MEMORY macros.
    // In Rust, if we want to avoid allocations here when distances are not requested,
    // the internal BFS needs to be adapted or `gridDiskUnsafe` needs to be callable
    // without distances.
    // For now, let's assume if distances are None, we use a temp vec.
    let mut temp_distances = vec![0i32; max_k_size as usize];
    _grid_disk_distances_internal(origin, k, out_cells, &mut temp_distances, max_k_size, 0)
  }
}

/// Produces H3 cells within `k` distance of the origin cell.
/// Handles pentagons correctly.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `k` - The k-ring radius (k >= 0).
/// * `out_cells` - Output array for H3 cells. Must be pre-allocated to `max_grid_disk_size(k)` and zero-initialized.
pub fn grid_disk(origin: H3Index, k: i32, out_cells: &mut [H3Index]) -> Result<(), H3Error> {
  grid_disk_distances(origin, k, out_cells, None)
}

/// "Unsafe" version of gridDiskDistances, outputting cells and their distances.
/// Output behavior is undefined or an error is returned when a pentagon
/// or its distortion area is encountered.
///
/// Cells are added to `out_cells` in order of increasing distance from the origin.
/// Within each k-ring, the order is spiral, starting from the I-axis direction
/// from the previous ring's I-axis cell.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `k` - The k-ring radius (k >= 0).
/// * `out_cells` - Output array for H3 cells.
/// * `out_distances_opt` - Optional output array for distances.
///
/// # Returns
/// `Ok(())` on success, `H3Error::Pentagon` if a pentagon is encountered,
/// or other `H3Error` on failure.
pub fn grid_disk_distances_unsafe(
  mut origin: H3Index,
  k: i32,
  out_cells: &mut [H3Index],
  mut out_distances_opt: Option<&mut [i32]>,
) -> Result<(), H3Error> {
  if k < 0 {
    return Err(H3Error::Domain);
  }

  let max_size = max_grid_disk_size(k)? as usize;
  if out_cells.len() < max_size {
    return Err(H3Error::MemoryBounds);
  }
  if let Some(ref dists) = out_distances_opt {
    if dists.len() < max_size {
      return Err(H3Error::MemoryBounds);
    }
  }

  let mut idx = 0;
  out_cells[idx] = origin; // Add original origin
  if let Some(ref mut dists_slice) = out_distances_opt {
    dists_slice[idx] = 0;
  }
  idx += 1;

  // C: if (H3_EXPORT(isPentagon)(origin)) { return E_PENTAGON; }
  // This check is after adding origin. If k=0, this is fine. If k>0, this means error.
  if is_pentagon(origin) {
    if k > 0 {
      // Only error if we intend to find neighbors
      for i in idx..out_cells.len() {
        out_cells[i] = H3_NULL;
      }
      if let Some(ref mut ds) = out_distances_opt {
        for i in idx..ds.len() {
          ds[i] = -1;
        }
      }
      return Err(H3Error::Pentagon);
    }
    // If k=0 and origin is pentagon, it's a valid output of 1 cell.
    // Fill rest if buffer is larger.
    if k == 0 {
      for i in idx..out_cells.len() {
        out_cells[i] = H3_NULL;
      }
      if let Some(ref mut ds) = out_distances_opt {
        for i in idx..ds.len() {
          ds[i] = -1;
        }
      }
      return Ok(());
    }
  }
  if k == 0 {
    // Handled origin, and if it was pentagon k>0 would have errored.
    for i in idx..out_cells.len() {
      out_cells[i] = H3_NULL;
    }
    if let Some(ref mut ds) = out_distances_opt {
      for i in idx..ds.len() {
        ds[i] = -1;
      }
    }
    return Ok(());
  }

  let mut ring_k_num = 1; // current ring (ring in C)
  let mut c_direction_idx = 0; // current side of the ring (direction in C)
  let mut pos_on_side = 0; // current position on the side of the ring (i in C)

  // This is `rotations` in C code, accumulates throughout.
  let mut path_accumulated_rotations: i32 = 0;

  const C_NEXT_RING_DIRECTION: Direction = Direction::IAxes;
  const C_SIDE_DIRECTIONS: [Direction; 6] = [
    Direction::JAxes,
    Direction::JkAxes,
    Direction::KAxes,
    Direction::IkAxes,
    Direction::IAxes,
    Direction::IjAxes,
  ];

  while ring_k_num <= k {
    if c_direction_idx == 0 && pos_on_side == 0 {
      // Start of a new ring (or first step of first ring)
      // Move `origin_param` (which is current_H) to the I-Axis start of this ring.
      let mut step_rot = 0; // Rotations for this specific step
      let neighbor_result = h3_neighbor_rotations(
        origin,
        C_NEXT_RING_DIRECTION,
        &mut step_rot,
        &mut origin, // origin_param is updated
      );
      if neighbor_result.is_err() {
        return neighbor_result;
      } // Propagate E_PENTAGON if h3nr fails
      path_accumulated_rotations = (path_accumulated_rotations + step_rot) % 6;
      if path_accumulated_rotations < 0 {
        path_accumulated_rotations += 6;
      }

      // C: if (H3_EXPORT(isPentagon)(origin)) { return E_PENTAGON; }
      // `origin_param` is now the start-of-ring cell.
      if is_pentagon(origin) {
        return Err(H3Error::Pentagon);
      }
    }

    // Now, `origin_param` is a cell on the current ring.
    // The C code adds it to output *before* the next move for the side.
    // But my loop adds it *after* the NEXT_RING_DIRECTION step.
    // The C code is:
    //   IF start_of_ring: move_IAXIS, update origin_param, check isPentagon(origin_param)
    //   move_SIDE_DIR, update origin_param
    //   out[idx++] = origin_param (this is the cell after move_SIDE_DIR)
    //   update i, direction
    //   check isPentagon(origin_param) (the cell just added)

    // My version (from previous successful C test alignment):
    // current_h_on_ring (start of ring) added.
    // Then loop: move current_h_on_ring, then add.
    // This seems to be where the difference is from the C snippet.

    // Let's strictly follow the C snippet's update and check order for `origin_param`
    let mut step_rot_side = 0;
    let mut dir_for_side_step = C_SIDE_DIRECTIONS[c_direction_idx as usize];
    for _ in 0..path_accumulated_rotations {
      // Pre-rotate like C kRing
      dir_for_side_step = crate::coords::ijk::_rotate60_ccw(dir_for_side_step);
    }

    let neighbor_result_side = h3_neighbor_rotations(
      origin, // current cell on ring
      dir_for_side_step,
      &mut step_rot_side,
      &mut origin, // origin_param becomes next cell on ring
    );
    if neighbor_result_side.is_err() {
      return neighbor_result_side;
    }
    path_accumulated_rotations = (path_accumulated_rotations + step_rot_side) % 6;
    if path_accumulated_rotations < 0 {
      path_accumulated_rotations += 6;
    }

    // Now origin_param is the new cell. Add it.
    if idx >= max_size {
      return Err(H3Error::MemoryBounds);
    }
    out_cells[idx] = origin;
    if let Some(ref mut dists_slice) = out_distances_opt {
      dists_slice[idx] = ring_k_num;
    }
    idx += 1;

    pos_on_side += 1;
    if pos_on_side == ring_k_num {
      // End of this side
      pos_on_side = 0;
      c_direction_idx += 1;
      if c_direction_idx == 6 {
        // End of this ring
        c_direction_idx = 0;
        ring_k_num += 1; // Move to next ring
      }
    }

    // C: if (H3_EXPORT(isPentagon)(origin)) { return E_PENTAGON; }
    // `origin_param` is the cell just added.
    if is_pentagon(origin) {
      return Err(H3Error::Pentagon);
    }
  }

  if idx != max_size {
    return Err(H3Error::Failed);
  }
  Ok(())
}

pub fn grid_disk_unsafe(origin: H3Index, k: i32, out: &mut [H3Index]) -> Result<(), H3Error> {
  grid_disk_distances_unsafe(origin, k, out, None)
}

/// Returns the "hollow" ring of H3 cells at exactly grid distance `k` from the origin cell.
/// Behavior is undefined if a pentagon is encountered.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `k` - The k-ring radius (k >= 0).
/// * `out_cells` - Output array for H3 cells.
///                 Must be pre-allocated to `1` if `k=0`, or `6*k` if `k > 0`.
///                 The array is filled contiguously.
///
/// # Returns
/// `Ok(())` on success, `H3Error::Pentagon` if a pentagon is encountered,
/// `H3Error::Domain` if `k < 0`, or `H3Error::MemoryBounds` if `out_cells` is too small.
pub fn grid_ring_unsafe(origin: H3Index, k: i32, out_cells: &mut [H3Index]) -> Result<(), H3Error> {
  if k < 0 {
    return Err(H3Error::Domain);
  }

  let expected_size = if k == 0 { 1 } else { (6 * k) as usize };
  if out_cells.len() < expected_size {
    return Err(H3Error::MemoryBounds);
  }
  // It's good practice to zero out if C relies on it for partially filled arrays,
  // but this function should fill contiguously up to expected_size.
  // for item in out_cells.iter_mut().take(expected_size) { *item = H3_NULL; }

  if k == 0 {
    out_cells[0] = origin;
    return Ok(());
  }

  let mut idx = 0;
  let mut rotations: i32 = 0;
  let mut current_h = origin;

  if is_pentagon(current_h) {
    return Err(H3Error::Pentagon);
  }

  // Move to the start of the target ring k
  for _ring_num in 0..k {
    let err = h3_neighbor_rotations(current_h, Direction::IAxes, &mut rotations, &mut current_h);
    if err.is_err() {
      return err;
    }
    if is_pentagon(current_h) {
      return Err(H3Error::Pentagon);
    }
  }

  // current_h is now the first cell in ring k
  let first_cell_in_ring = current_h;

  out_cells[idx] = current_h;
  idx += 1;

  const RING_TRAVERSAL_DIRECTIONS: [Direction; 6] = [
    Direction::JAxes,
    Direction::JkAxes,
    Direction::KAxes,
    Direction::IkAxes,
    Direction::IAxes,
    Direction::IjAxes,
  ];

  for side_idx in 0..6 {
    // For each of the 6 sides of the ring
    // For each cell along this side.
    // For k=1, inner loop runs 1 time per side. Total 6 steps. out[0] to out[5].
    // For k=2, inner loop runs 2 times per side. Total 12 steps. out[0] to out[11].
    for pos_on_side in 0..k {
      let err = h3_neighbor_rotations(
        current_h,
        RING_TRAVERSAL_DIRECTIONS[side_idx],
        &mut rotations,
        &mut current_h,
      );
      if err.is_err() {
        return err;
      }
      if is_pentagon(current_h) {
        return Err(H3Error::Pentagon);
      }

      // Add to output if it's not the very last cell that would duplicate the start
      if idx < expected_size {
        // This ensures we don't write past the allocated slice
        // and also handles the "don't add last duplicate" implicitly
        // if expected_size is exactly 6*k.
        out_cells[idx] = current_h;
        idx += 1;
      } else {
        // We've filled 6*k cells. If current_h is not the first_cell_in_ring,
        // it means distortion occurred.
        if current_h != first_cell_in_ring {
          return Err(H3Error::Pentagon);
        }
        // If it IS first_cell_in_ring, we're good, loop will terminate.
      }
    }
  }

  // Final check: after traversing all 6k cells, current_h should be back at the first_cell_in_ring.
  // This is implicitly handled by the logic above; if it wasn't, idx wouldn't be expected_size
  // or the else branch in the loop would have triggered.
  if idx != expected_size || (k > 0 && current_h != first_cell_in_ring) {
    // If idx didn't reach expected_size, something went wrong (e.g. an early pentagon error not caught)
    // If current_h isn't first_cell_in_ring, implies distortion
    return Err(H3Error::Pentagon); // Or Failed for idx mismatch
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::indexing::lat_lng_to_cell;
  use crate::latlng::_set_geo_degs;
  use crate::types::LatLng;
  use std::collections::HashSet;

  #[test]
  fn test_max_grid_disk_size() {
    assert_eq!(max_grid_disk_size(0), Ok(1));
    assert_eq!(max_grid_disk_size(1), Ok(7));
    assert_eq!(max_grid_disk_size(2), Ok(19));
    assert_eq!(max_grid_disk_size(-1), Err(H3Error::Domain));
    // Test the cap
    assert_eq!(
      max_grid_disk_size(K_ALL_CELLS_AT_RES_15),
      Ok(crate::constants::NUM_CELLS_MAX_RES)
    );
    assert_eq!(
      max_grid_disk_size(K_ALL_CELLS_AT_RES_15 + 1000),
      Ok(crate::constants::NUM_CELLS_MAX_RES)
    );
  }

  #[test]
  fn test_grid_disk_distances_k0() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let origin = lat_lng_to_cell(&geo, 5).unwrap();

    let k = 0;
    let max_size = max_grid_disk_size(k).unwrap() as usize;
    let mut cells = vec![H3_NULL; max_size];
    let mut distances = vec![0i32; max_size];

    assert!(grid_disk_distances(origin, k, &mut cells, Some(&mut distances)).is_ok());

    let mut count = 0;
    for i in 0..max_size {
      if cells[i] != H3_NULL {
        count += 1;
        assert_eq!(cells[i], origin);
        assert_eq!(distances[i], 0);
      }
    }
    assert_eq!(count, 1, "k=0 should return only origin");
  }

  #[test]
  fn test_grid_disk_k1() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let origin = lat_lng_to_cell(&geo, 5).unwrap();

    let k = 1;
    let max_size = max_grid_disk_size(k).unwrap() as usize;
    let mut cells = vec![H3_NULL; max_size];

    assert!(grid_disk(origin, k, &mut cells).is_ok());

    let mut found_cells = HashSet::new();
    let mut origin_found = false;
    for i in 0..max_size {
      if cells[i] != H3_NULL {
        assert!(crate::h3_index::inspection::is_valid_cell(cells[i]));
        found_cells.insert(cells[i]);
        if cells[i] == origin {
          origin_found = true;
        }
      }
    }
    assert_eq!(found_cells.len(), 7, "k=1 should return 7 unique cells for a hexagon");
    assert!(origin_found, "Origin cell should be in k=1 disk");
  }

  #[test]
  fn test_grid_disk_distances_pentagon_k1() {
    // Base cell 4 is a pentagon. Get its center child at res 1.
    let pent_origin_res0 = crate::base_cells::baseCellNumToCell(4);
    let pent_origin = crate::hierarchy::cell_to_center_child(pent_origin_res0, 1).unwrap();
    assert!(is_pentagon(pent_origin));

    let k = 1;
    let max_size = max_grid_disk_size(k).unwrap() as usize; // Should be 7
    let mut cells = vec![H3_NULL; max_size];
    let mut distances = vec![0i32; max_size];

    assert!(grid_disk_distances(pent_origin, k, &mut cells, Some(&mut distances)).is_ok());

    let mut found_cells_map = std::collections::HashMap::new();
    let mut cell_count = 0;
    for i in 0..max_size {
      if cells[i] != H3_NULL {
        cell_count += 1;
        assert!(crate::h3_index::inspection::is_valid_cell(cells[i]));
        found_cells_map.insert(cells[i], distances[i]);
      }
    }
    // Pentagon has 5 neighbors + itself = 6 cells in k=1 disk
    assert_eq!(cell_count, 6, "k=1 from pentagon should return 6 cells");
    assert_eq!(
      found_cells_map.len(),
      6,
      "k=1 from pentagon should return 6 unique cells"
    );
    assert_eq!(
      found_cells_map.get(&pent_origin),
      Some(&0),
      "Origin pentagon should be at distance 0"
    );

    for (cell, dist) in found_cells_map {
      if cell != pent_origin {
        assert_eq!(dist, 1, "Neighbor cells should be at distance 1");
      }
    }
  }

  #[test]
  fn test_grid_disk_unsafe_k1_hex() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let origin = lat_lng_to_cell(&geo, 5).unwrap();
    assert!(!is_pentagon(origin));

    let k = 1;
    let max_size = max_grid_disk_size(k).unwrap() as usize;
    let mut cells = vec![H3_NULL; max_size];
    let mut distances = vec![0i32; max_size];

    assert!(grid_disk_distances_unsafe(origin, k, &mut cells, Some(&mut distances)).is_ok());

    let mut found_cells = HashSet::new();
    let mut cell_count = 0;
    let mut dist_sum = 0;
    for i in 0..max_size {
      if cells[i] != H3_NULL {
        cell_count += 1;
        assert!(crate::h3_index::inspection::is_valid_cell(cells[i]));
        found_cells.insert(cells[i]);
        if cells[i] == origin {
          assert_eq!(distances[i], 0);
        } else {
          assert_eq!(distances[i], 1);
        }
        dist_sum += distances[i];
      }
    }
    assert_eq!(cell_count, 7);
    assert_eq!(found_cells.len(), 7);
    assert_eq!(dist_sum, 6); // 0 (origin) + 6 * 1 (neighbors)
  }

  #[test]
  fn test_grid_disk_unsafe_pentagon_origin() {
    let pent_origin_res0 = crate::base_cells::baseCellNumToCell(4);
    let pent_origin = crate::hierarchy::parent_child::cell_to_center_child(pent_origin_res0, 1).unwrap();
    assert!(is_pentagon(pent_origin));

    let k = 1;
    let max_size = max_grid_disk_size(k).unwrap() as usize;
    let mut cells = vec![H3_NULL; max_size];

    assert_eq!(grid_disk_unsafe(pent_origin, k, &mut cells), Err(H3Error::Pentagon));
  }

  #[test]
  fn test_grid_disk_unsafe_near_pentagon() {
    // Find a cell that is a direct neighbor of a pentagon
    let pent_origin_res0 = crate::base_cells::baseCellNumToCell(4); // BC 4 is pentagon
    let pent_origin_res1 = crate::hierarchy::parent_child::cell_to_center_child(pent_origin_res0, 1).unwrap();

    let mut pent_neighbors = vec![H3_NULL; 7];
    grid_disk(pent_origin_res1, 1, &mut pent_neighbors).unwrap(); // Use safe version to get neighbors

    let mut hex_near_pent = H3_NULL;
    for neighbor in pent_neighbors {
      if neighbor != H3_NULL && neighbor != pent_origin_res1 {
        hex_near_pent = neighbor;
        break;
      }
    }
    assert_ne!(hex_near_pent, H3_NULL, "Could not find a hex neighbor of the pentagon");
    assert!(!is_pentagon(hex_near_pent));

    // k=1 from this hex_near_pent should not hit the pentagon with gridDiskUnsafe
    // because gridDiskUnsafe's first ring traversal does not include the origin.
    // The unsafe k-ring only fails if one of the *output cells* is a pentagon.
    let k = 1;
    let max_size_k1 = max_grid_disk_size(k).unwrap() as usize;
    let mut cells_k1 = vec![H3_NULL; max_size_k1];
    assert_eq!(
      grid_disk_unsafe(hex_near_pent, k, &mut cells_k1),
      Err(H3Error::Pentagon),
      "k=1 from hex near pentagon should error in unsafe because the disk includes the pentagon"
    );

    // k=2 from this hex_near_pent *will* include the pentagon in its output, so unsafe should fail.
    let k_unsafe_fail = 2;
    let max_size_k_unsafe_fail = max_grid_disk_size(k_unsafe_fail).unwrap() as usize;
    let mut cells_k_unsafe_fail = vec![H3_NULL; max_size_k_unsafe_fail];
    assert_eq!(
      grid_disk_unsafe(hex_near_pent, k_unsafe_fail, &mut cells_k_unsafe_fail),
      Err(H3Error::Pentagon),
      "k=2 from hex near pentagon should error in unsafe as pentagon is in output"
    );
  }

  #[test]
  fn test_grid_ring_unsafe_k0() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let origin = lat_lng_to_cell(&geo, 9).unwrap();

    let mut ring_cells = [H3_NULL; 1];
    assert!(grid_ring_unsafe(origin, 0, &mut ring_cells).is_ok());
    assert_eq!(ring_cells[0], origin);
  }

  #[test]
  fn test_grid_ring_unsafe_k1_hex() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let origin = lat_lng_to_cell(&geo, 9).unwrap();
    assert!(!is_pentagon(origin));

    let k = 1;
    let expected_size = (6 * k) as usize;
    let mut ring_cells = vec![H3_NULL; expected_size];

    assert!(grid_ring_unsafe(origin, k, &mut ring_cells).is_ok());

    let mut found_cells = HashSet::new();
    let mut cell_count = 0;
    for i in 0..expected_size {
      if ring_cells[i] != H3_NULL {
        // Should all be non-null for successful hex case
        cell_count += 1;
        assert!(crate::h3_index::inspection::is_valid_cell(ring_cells[i]));
        assert_ne!(ring_cells[i], origin, "Ring cells should not be the origin for k>0");
        found_cells.insert(ring_cells[i]);
      }
    }
    assert_eq!(cell_count, expected_size);
    assert_eq!(found_cells.len(), expected_size, "All ring cells should be unique");

    // Compare with grid_disk to verify content (excluding origin)
    let max_disk_size = max_grid_disk_size(k).unwrap() as usize;
    let mut disk_cells = vec![H3_NULL; max_disk_size];
    let mut disk_distances = vec![0i32; max_disk_size];
    assert!(grid_disk_distances(origin, k, &mut disk_cells, Some(&mut disk_distances)).is_ok());

    for cell_in_disk in disk_cells.iter().zip(disk_distances.iter()) {
      if *cell_in_disk.0 != H3_NULL && *cell_in_disk.0 != origin {
        // Only neighbors
        assert_eq!(*cell_in_disk.1, k, "Cell from disk should be at distance k");
        assert!(
          found_cells.contains(cell_in_disk.0),
          "Cell from disk's k-ring found in grid_ring_unsafe output"
        );
      }
    }
  }

  #[test]
  fn test_grid_ring_unsafe_k2_hex() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419); // SF
    let origin = lat_lng_to_cell(&geo, 5).unwrap();
    assert!(!is_pentagon(origin));

    let k = 2;
    let expected_size = (6 * k) as usize; // 12 cells
    let mut ring_cells = vec![H3_NULL; expected_size];

    assert!(grid_ring_unsafe(origin, k, &mut ring_cells).is_ok());

    let mut found_cells = HashSet::new();
    for i in 0..expected_size {
      assert_ne!(ring_cells[i], H3_NULL, "Cell {} should be populated", i);
      assert!(crate::h3_index::inspection::is_valid_cell(ring_cells[i]));
      assert_ne!(ring_cells[i], origin);
      assert!(
        found_cells.insert(ring_cells[i]),
        "Cell {:x} was duplicate in ring",
        ring_cells[i].0
      );
    }
    assert_eq!(
      found_cells.len(),
      expected_size,
      "All ring cells should be unique for k=2 hex"
    );
  }

  #[test]
  fn test_grid_ring_unsafe_on_pentagon_origin() {
    let pent_origin_res0 = crate::base_cells::baseCellNumToCell(4);
    let pent_origin = crate::hierarchy::parent_child::cell_to_center_child(pent_origin_res0, 1).unwrap();
    assert!(is_pentagon(pent_origin));

    let k_val: i32 = 1; // Explicitly i32

    // Calculate size for the vec, ensuring k_val is treated as i32 for the multiplication
    // and then the result is cast to usize.
    let vec_size: usize = if k_val == 0 { 1 } else { (6 * k_val) as usize };
    let mut ring_cells = vec![H3_NULL; vec_size];

    assert_eq!(
      grid_ring_unsafe(pent_origin, k_val, &mut ring_cells),
      Err(H3Error::Pentagon)
    );
  }

  #[test]
  fn test_grid_ring_unsafe_near_pentagon_failure() {
    let pent_origin_res0 = crate::base_cells::baseCellNumToCell(4); // BC 4 is pentagon
    let pent_origin_res1 = crate::hierarchy::parent_child::cell_to_center_child(pent_origin_res0, 1).unwrap();

    let mut pent_neighbors_disk = vec![H3_NULL; max_grid_disk_size(1).unwrap() as usize];
    // Use the public grid_disk which uses the safe internal BFS, ensuring it works around pentagons
    grid_disk(pent_origin_res1, 1, &mut pent_neighbors_disk).unwrap();

    let mut hex_near_pent = H3_NULL;
    for neighbor_h3 in pent_neighbors_disk {
      if neighbor_h3 != H3_NULL && neighbor_h3 != pent_origin_res1 {
        hex_near_pent = neighbor_h3;
        break;
      }
    }
    assert_ne!(hex_near_pent, H3_NULL, "Could not find a hex neighbor of the pentagon");
    assert!(!is_pentagon(hex_near_pent));

    // k=1 ring from hex_near_pent should NOT hit the pentagon.
    let k1_val: i32 = 1; // Explicitly type k1_val as i32
    let ring_k1_size: usize = (6 * k1_val) as usize; // Calculate size, ensuring cast to usize
    let mut ring_k1 = vec![H3_NULL; ring_k1_size];
    assert_eq!(
      grid_ring_unsafe(hex_near_pent, k1_val, &mut ring_k1),
      Err(H3Error::Pentagon),
      "k=1 ring from hex near pentagon should error in unsafe because the ring includes the pentagon"
    );
    // k=2 ring from hex_near_pent WILL encounter the pentagon.
    let k2_val: i32 = 2; // Explicitly type k2_val as i32
    let ring_k2_size: usize = (6 * k2_val) as usize; // Calculate size, ensuring cast to usize
    let mut ring_k2 = vec![H3_NULL; ring_k2_size];
    assert_eq!(
      grid_ring_unsafe(hex_near_pent, k2_val, &mut ring_k2),
      Err(H3Error::Pentagon),
      "Line 825 error site"
    );
  }

  #[test]
  fn test_grid_ring_unsafe_invalid_k() {
    let origin = H3Index(0x85283473fffffff);
    let mut cells = [H3_NULL; 1];
    assert_eq!(grid_ring_unsafe(origin, -1, &mut cells), Err(H3Error::Domain));
  }

  #[test]
  fn test_grid_ring_unsafe_output_too_small() {
    let origin = H3Index(0x85283473fffffff);
    let mut cells = [H3_NULL; 5]; // k=1 needs 6
    assert_eq!(grid_ring_unsafe(origin, 1, &mut cells), Err(H3Error::MemoryBounds));
  }

  // TODO: Add tests for gridDiskUnsafe once implemented.
  // TODO: Add tests comparing unsafe and safe versions.
  // TODO: Add tests for larger k values and edge cases (poles, transmeridian).
}
