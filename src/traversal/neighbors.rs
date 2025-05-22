// src/traversal/neighbors.rs

use crate::base_cells::{
  _base_cell_is_cw_offset,
  _get_base_cell_direction, // Corrected name from C-style to snake_case
  _get_base_cell_neighbor,  // Corrected name
  _is_base_cell_pentagon,
  _is_base_cell_polar_pentagon,
  BASE_CELL_DATA,                // Corrected to uppercase
  BASE_CELL_NEIGHBORS,           // Corrected to uppercase
  BASE_CELL_NEIGHBOR_60CCW_ROTS, // Corrected to uppercase
  INVALID_BASE_CELL,
};
use crate::constants::{H3_CELL_MODE, MAX_H3_RES, NUM_BASE_CELLS};
use crate::coords::ijk::_rotate60_ccw;
use crate::types::Direction;
use crate::{is_pentagon, H3_NULL};

use crate::h3_index::{
  _h3_leading_non_zero_digit, _h3_rotate60_ccw, _h3_rotate60_cw, _h3_rotate_pent60_ccw, get_base_cell, get_index_digit,
  get_resolution, is_resolution_class_iii, set_base_cell, set_index_digit,
};
use crate::types::{H3Error, H3Index};

// Constants for new digit/adjustment tables (from C algos.c)
// These could also go into constants.rs or a dedicated algos_consts.rs

#[rustfmt::skip]
const NEW_DIGIT_II: [[Direction; 7]; 7] = [
    [Direction::Center, Direction::KAxes,  Direction::JAxes,  Direction::JkAxes,  Direction::IAxes,  Direction::IkAxes,  Direction::IjAxes],
    [Direction::KAxes,  Direction::IAxes,  Direction::JkAxes,  Direction::IjAxes,  Direction::IkAxes,  Direction::JAxes,   Direction::Center],
    [Direction::JAxes,  Direction::JkAxes, Direction::KAxes,  Direction::IAxes,   Direction::IjAxes,  Direction::Center,  Direction::IkAxes],
    [Direction::JkAxes, Direction::IjAxes, Direction::IAxes,  Direction::IkAxes,  Direction::Center,  Direction::KAxes,   Direction::JAxes],
    [Direction::IAxes,  Direction::IkAxes, Direction::IjAxes, Direction::Center,  Direction::JAxes,   Direction::JkAxes,  Direction::KAxes],
    [Direction::IkAxes, Direction::JAxes,  Direction::Center, Direction::KAxes,   Direction::JkAxes,  Direction::IjAxes,  Direction::IAxes],
    [Direction::IjAxes, Direction::Center, Direction::IkAxes, Direction::JAxes,   Direction::KAxes,   Direction::IAxes,   Direction::JkAxes]
];

#[rustfmt::skip]
const NEW_ADJUSTMENT_II: [[Direction; 7]; 7] = [
    [Direction::Center, Direction::Center, Direction::Center, Direction::Center,   Direction::Center, Direction::Center,   Direction::Center],
    [Direction::Center, Direction::KAxes,  Direction::Center, Direction::KAxes,    Direction::Center, Direction::IkAxes,   Direction::Center],
    [Direction::Center, Direction::Center, Direction::JAxes,  Direction::JkAxes,   Direction::Center, Direction::Center,   Direction::JAxes],
    [Direction::Center, Direction::KAxes,  Direction::JkAxes, Direction::JkAxes,   Direction::Center, Direction::Center,   Direction::Center],
    [Direction::Center, Direction::Center, Direction::Center, Direction::Center,   Direction::IAxes,  Direction::IAxes,    Direction::IjAxes],
    [Direction::Center, Direction::IkAxes, Direction::Center, Direction::Center,   Direction::IAxes,  Direction::IkAxes,   Direction::Center],
    [Direction::Center, Direction::Center, Direction::JAxes,  Direction::Center,   Direction::IjAxes, Direction::Center,   Direction::IjAxes]
];

#[rustfmt::skip]
const NEW_DIGIT_III: [[Direction; 7]; 7] = [
    [Direction::Center, Direction::KAxes,  Direction::JAxes,  Direction::JkAxes,  Direction::IAxes,  Direction::IkAxes,  Direction::IjAxes],
    [Direction::KAxes,  Direction::JAxes,  Direction::JkAxes, Direction::IAxes,   Direction::IkAxes,  Direction::IjAxes,  Direction::Center],
    [Direction::JAxes,  Direction::JkAxes, Direction::IAxes,  Direction::IkAxes,  Direction::IjAxes,  Direction::Center,  Direction::KAxes],
    [Direction::JkAxes, Direction::IAxes,  Direction::IkAxes, Direction::IjAxes,  Direction::Center,  Direction::KAxes,   Direction::JAxes],
    [Direction::IAxes,  Direction::IkAxes, Direction::IjAxes, Direction::Center,  Direction::KAxes,   Direction::JAxes,   Direction::JkAxes],
    [Direction::IkAxes, Direction::IjAxes, Direction::Center, Direction::KAxes,   Direction::JAxes,   Direction::JkAxes,  Direction::IAxes],
    [Direction::IjAxes, Direction::Center, Direction::KAxes,  Direction::JAxes,   Direction::JkAxes,  Direction::IAxes,   Direction::IkAxes]
];

#[rustfmt::skip]
const NEW_ADJUSTMENT_III: [[Direction; 7]; 7] = [
    [Direction::Center, Direction::Center, Direction::Center, Direction::Center,   Direction::Center, Direction::Center,   Direction::Center],
    [Direction::Center, Direction::KAxes,  Direction::Center, Direction::JkAxes,   Direction::Center, Direction::KAxes,    Direction::Center],
    [Direction::Center, Direction::Center, Direction::JAxes,  Direction::JAxes,    Direction::Center, Direction::Center,   Direction::IjAxes],
    [Direction::Center, Direction::JkAxes, Direction::JAxes,  Direction::JkAxes,   Direction::Center, Direction::Center,   Direction::Center],
    [Direction::Center, Direction::Center, Direction::Center, Direction::Center,   Direction::IAxes,  Direction::IkAxes,   Direction::IAxes],
    [Direction::Center, Direction::KAxes,  Direction::Center, Direction::Center,   Direction::IkAxes, Direction::IkAxes,   Direction::Center],
    [Direction::Center, Direction::Center, Direction::IjAxes, Direction::Center,   Direction::IAxes,  Direction::Center,   Direction::IjAxes]
];

/// Returns the hexagon index neighboring the origin, in the direction `dir`.
///
/// # Arguments
/// * `origin` - Origin H3Index.
/// * `dir` - Direction to move in (1-6).
/// * `rotations` - Input: current number of 60 degree CCW rotations relative to normal face IJK.
///                 Output: new number of rotations after traversal.
/// * `out` - H3Index of the specified neighbor.
///
/// # Returns
/// `Ok(())` on success, `H3Error` on failure.
pub(crate) fn h3_neighbor_rotations(
  origin: H3Index,
  mut dir: Direction, // mut to allow rotation
  rotations: &mut i32,
  out: &mut H3Index,
) -> Result<(), H3Error> {
  let mut current_h = origin; // Use a mutable copy for modifications

  if dir == Direction::Center || dir == Direction::InvalidDigit {
    // This indicates an invalid H3 an航空方向
    // The C code returns E_FAILED (1) for this.
    return Err(H3Error::Failed); // Or a more specific error like DomainError
  }

  // Ensure that rotations is modulo'd by 6 before any possible addition,
  // to protect against signed integer overflow.
  *rotations %= 6;
  if *rotations < 0 {
    *rotations += 6;
  } // Ensure positive modulo

  for _i in 0..*rotations {
    dir = _rotate60_ccw(dir);
  }

  let mut new_rotations = 0;
  let old_base_cell = get_base_cell(current_h);

  // Base cells less than zero or >= NUM_BASE_CELLS can not be represented in an index.
  // This should ideally be caught by isValidCell if origin was an external input.
  if old_base_cell < 0 || old_base_cell >= (NUM_BASE_CELLS as i32) {
    return Err(H3Error::CellInvalid);
  }
  let old_leading_digit = _h3_leading_non_zero_digit(current_h);

  // Adjust the indexing digits and, if needed, the base cell.
  let mut r = get_resolution(current_h) - 1; // Resolution of the parent cell.
  loop {
    if r == -1 {
      // `dir` here is the original input `dir` parameter possibly rotated by `*rotations`
      if _is_base_cell_pentagon(old_base_cell) && dir == Direction::KAxes {
        // Explicitly cannot move in K_AXES direction from a pentagon at base cell level.
        // This matches behavior of some H3 library versions.
        return Err(H3Error::Pentagon);
      }

      set_base_cell(
        &mut current_h,
        BASE_CELL_NEIGHBORS[old_base_cell as usize][dir as usize],
      );
      new_rotations = BASE_CELL_NEIGHBOR_60CCW_ROTS[old_base_cell as usize][dir as usize];

      if get_base_cell(current_h) == INVALID_BASE_CELL {
        // This implies origin was a hexagon, and 'dir' was K_AXES pointing to a pentagon.
        // Or, origin was pentagon, 'dir' was not K_AXES, but the neighbor lookup still failed (should not happen).
        // If origin was pentagon and dir was K, it's already handled above.
        // So, this path is for: Hexagon origin, dir=K, new_base_cell is a pentagon.

        // The C "fixup" logic:
        // We were trying to go in `dir` (which was K) from `old_base_cell` (hexagon).
        // This led to an invalid base cell, meaning the K-neighbor is a pentagon.
        // We must use an alternative path, IK_AXES_DIGIT, and adjust rotations.
        // Note: This adjustment is for the *digits* of current_h, which still holds origin's digits.

        // Set to IK neighbor instead
        set_base_cell(
          &mut current_h, // current_h at this point has old_base_cell, but its digits are from origin
          BASE_CELL_NEIGHBORS[old_base_cell as usize][Direction::IkAxes as usize],
        );
        new_rotations = BASE_CELL_NEIGHBOR_60CCW_ROTS[old_base_cell as usize][Direction::IkAxes as usize];

        // The digits of current_h (still origin's digits) need to be rotated as if we went via IK.
        // The C code also applies _h3Rotate60ccw(out) to the digits here.
        current_h = _h3_rotate60_ccw(current_h);
        // And the overall rotation count for subsequent levels is also adjusted.
        *rotations = (*rotations + 1) % 6;
      }
      break; // Exited loop for r
    } else {
      let old_digit = get_index_digit(current_h, r + 1);
      let next_dir: Direction;

      if old_digit == Direction::InvalidDigit {
        return Err(H3Error::CellInvalid); // Only possible on invalid input
      }

      if is_resolution_class_iii(r + 1) {
        set_index_digit(&mut current_h, r + 1, NEW_DIGIT_II[old_digit as usize][dir as usize]);
        next_dir = NEW_ADJUSTMENT_II[old_digit as usize][dir as usize];
      } else {
        set_index_digit(&mut current_h, r + 1, NEW_DIGIT_III[old_digit as usize][dir as usize]);
        next_dir = NEW_ADJUSTMENT_III[old_digit as usize][dir as usize];
      }

      if next_dir != Direction::Center {
        dir = next_dir;
        r -= 1;
      } else {
        // No more adjustment to perform
        break;
      }
    }
  }

  let new_base_cell = get_base_cell(current_h);
  if _is_base_cell_pentagon(new_base_cell) {
    let mut already_adjusted_k_subsequence = false;

    // Force rotation out of missing k-axes sub-sequence
    if _h3_leading_non_zero_digit(current_h) == Direction::KAxes {
      if old_base_cell != new_base_cell {
        // Traversed into the deleted k subsequence of a pentagon base cell.
        // Rotate out based on how we got here.
        if _base_cell_is_cw_offset(new_base_cell, BASE_CELL_DATA[old_base_cell as usize].home_fijk.face) {
          current_h = _h3_rotate60_cw(current_h);
        } else {
          current_h = _h3_rotate60_ccw(current_h);
        }
        already_adjusted_k_subsequence = true;
      } else {
        // Traversed into deleted k subsequence from within the same pentagon base cell.
        match old_leading_digit {
          Direction::Center => return Err(H3Error::Pentagon), // Undefined: K direction is deleted
          Direction::JkAxes => {
            current_h = _h3_rotate60_ccw(current_h);
            *rotations = (*rotations + 1) % 6;
          }
          Direction::IkAxes => {
            current_h = _h3_rotate60_cw(current_h);
            *rotations = (*rotations + 5) % 6; // +5 is equiv to -1 or 1 cw
          }
          _ => return Err(H3Error::Failed), // Should not occur, but fuzzer might find
        }
      }
    }

    for _i in 0..new_rotations {
      current_h = _h3_rotate_pent60_ccw(current_h);
    }

    // Account for differing orientation of the base cells
    if old_base_cell != new_base_cell {
      if _is_base_cell_polar_pentagon(new_base_cell) {
        if old_base_cell != 118 && old_base_cell != 8 && // Not from the two "aligned" neighbors of polar pentagons
                   _h3_leading_non_zero_digit(current_h) != Direction::JkAxes
        {
          *rotations = (*rotations + 1) % 6;
        }
      } else if _h3_leading_non_zero_digit(current_h) == Direction::IkAxes && !already_adjusted_k_subsequence {
        // Account for distortion by deleted k subsequence for non-polar pentagons
        *rotations = (*rotations + 1) % 6;
      }
    }
  } else {
    // Hexagon base cell
    for _i in 0..new_rotations {
      current_h = _h3_rotate60_ccw(current_h);
    }
  }

  *rotations = (*rotations + new_rotations) % 6;
  *out = current_h;

  Ok(())
}

/// Get the H3 `Direction` from an origin H3 cell to a neighboring H3 cell.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `destination` - The H3 cell for which the direction is desired.
///
/// # Returns
/// The `Direction` from origin to destination, or `Direction::InvalidDigit`
/// if the cells are not neighbors.
// This is `directionForNeighbor` in C.
pub(crate) fn direction_for_neighbor(origin: H3Index, destination: H3Index) -> Direction {
  // eprintln!(
  //   "  DI_FOR_NB: direction_for_neighbor called with origin=0x{:x}, destination=0x{:x}",
  //   origin.0, destination.0
  // );

  if origin == destination {
    // eprintln!("    DI_FOR_NB: Origin and destination are the same. Returning Center.");
    return Direction::Center;
  }

  let is_pent = is_pentagon(origin); // is_pentagon from h3_index::inspection
  // eprintln!("    DI_FOR_NB: is_pentagon(origin) = {}", is_pent);

  let start_dir_u8: u8 = if is_pent {
    Direction::JAxes as u8
  } else {
    Direction::KAxes as u8
  };
  // Valid H3 directions for neighbors are K, J, JK, I, IK, IJ (1 through 6)
  // Loop from start_dir_u8 up to IJ_AXES_DIGIT (6)
  for dir_val in start_dir_u8..=(Direction::IjAxes as u8) {
    let direction = match Direction::try_from(dir_val) {
      Ok(d) => d,
      Err(_) => continue, // Should not happen if loop is 1..=6 or 2..=6
    };
    // eprintln!("    DI_FOR_NB: Testing direction {:?}", direction);

    let mut neighbor: H3Index = H3_NULL; // Initialize
    let mut rotations = 0; // Does not affect which neighbor is found, only its final orientation relative to origin's system

    let neighbor_result = h3_neighbor_rotations(origin, direction, &mut rotations, &mut neighbor);

    // eprintln!(
    //   "    DI_FOR_NB: h3_neighbor_rotations(0x{:x}, {:?}, ...) yielded H3 = 0x{:x}, rotations_out = {}, err = {:?}",
    //   origin.0, direction, neighbor.0, rotations, neighbor_result
    // );

    if neighbor_result.is_ok() && neighbor == destination {
      // eprintln!("    DI_FOR_NB: FOUND neighbor in direction {:?}", direction);
      return direction;
    } else {
      if neighbor_result.is_err() && neighbor_result != Err(H3Error::Pentagon) {
        // Log unexpected errors from h3_neighbor_rotations
        // eprintln!(
        //   "    DI_FOR_NB: h3_neighbor_rotations errored unexpectedly for direction {:?}: {:?}",
        //   direction, neighbor_result
        // );
      }
    }
  }

  // eprintln!("    DI_FOR_NB: Did not find direct neighbor after checking all valid directions.");
  Direction::InvalidDigit
}

/// Returns whether or not the provided H3 cells are neighbors.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `destination` - The destination H3 cell.
///
/// # Returns
/// `Ok(true)` if the cells are neighbors, `Ok(false)` otherwise.
/// Returns an `H3Error` if the input cells are invalid or not comparable.
// This is H3_EXPORT(areNeighborCells) in C.
pub fn are_neighbor_cells(origin: H3Index, destination: H3Index) -> Result<bool, H3Error> {
  // Make sure they're H3 cell indexes
  if crate::h3_index::get_mode(origin) != H3_CELL_MODE || crate::h3_index::get_mode(destination) != H3_CELL_MODE {
    return Err(H3Error::CellInvalid);
  }

  // Optimization: The C code has a fast path for cells sharing a common parent.
  // We can add that here too for performance parity.
  // For now, the direct check using direction_for_neighbor or gridDisk is functionally correct.

  // Cells cannot be neighbors with themselves
  if origin == destination {
    return Ok(false);
  }

  // Only cells in the same resolution can be neighbors
  if get_resolution(origin) != get_resolution(destination) {
    return Err(H3Error::ResMismatch);
  }

  // Check if `is_valid_cell` is necessary here.
  // The functions called by direction_for_neighbor (like get_base_cell) might
  // behave unexpectedly with truly garbage input bits even if mode/res match.
  // The C API version of areNeighborCells doesn't explicitly call isValidCell internally first,
  // but relies on downstream functions to handle or error out.
  // For a public Rust API, validating inputs is good.
  if !crate::h3_index::inspection::is_valid_cell(origin) || !crate::h3_index::inspection::is_valid_cell(destination) {
    return Err(H3Error::CellInvalid);
  }
  let dir = direction_for_neighbor(origin, destination);
  let result = dir != Direction::InvalidDigit;
  // eprintln!(
  //   "ARE_NEIGHBOR_CELLS_DEBUG: origin={:x}, dest={:x}, dir_found={:?}, result={}",
  //   origin.0, destination.0, dir, result
  // );
  // This relies heavily on the base cells
  Ok(direction_for_neighbor(origin, destination) != Direction::InvalidDigit)
  // From C tests, H3 library's _areNeighborCells function *always* calculates
  // the neighbor direction AND returns E_SUCCESS regardless of direction's validity

  // Original code
  // let mut out = 0;
  // neighborResult.map(|_ok_val| out = 1).map_err(|err| {
  //     if err == H3Error::NotNeighbors && !isPentagon(origin) {
  //         return err;
  //     }
  //     H3Error::Failed
  // })?;
  // Ok(out == 1)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::h3_index; // For direct access to setters/getters if needed for test setup
  use crate::indexing::lat_lng_to_cell; // Assuming this is pub from indexing
  use crate::latlng::_set_geo_degs;
  use crate::traversal::grid_disk::grid_disk;
  use crate::types::LatLng;
  use crate::types::H3_NULL; // Assuming grid_disk is ported and available

  #[test]
  fn test_direction_for_neighbor() {
    let mut sf_geo = LatLng::default();
    _set_geo_degs(&mut sf_geo, 37.779265, -122.419277);
    let sf_h3 = lat_lng_to_cell(&sf_geo, 9).unwrap();

    let mut ring = [H3_NULL; 7];
    assert!(grid_disk(sf_h3, 1, &mut ring).is_ok()); // Assuming grid_disk is ported

    let mut found_count = 0;
    for i in 0..7 {
      if ring[i] == H3_NULL || ring[i] == sf_h3 {
        continue;
      }
      found_count += 1;
      let dir = direction_for_neighbor(sf_h3, ring[i]);
      assert_ne!(
        dir,
        Direction::InvalidDigit,
        "Should find a direction for direct neighbor"
      );
      assert_ne!(dir, Direction::Center, "Direction should not be Center");

      // Round trip check: move from origin in `dir` should yield `ring[i]`
      let mut rotations = 0;
      let mut recovered_neighbor = H3_NULL;
      assert!(h3_neighbor_rotations(sf_h3, dir, &mut rotations, &mut recovered_neighbor).is_ok());
      assert_eq!(recovered_neighbor, ring[i], "Direction round trip mismatch");
    }
    assert_eq!(found_count, 6, "Found 6 distinct neighbors");

    // Test non-neighbor (e.g., from k-ring 2)
    let mut ring2 = [H3_NULL; 19];
    assert!(grid_disk(sf_h3, 2, &mut ring2).is_ok());
    let non_neighbor = ring2
      .iter()
      .find(|&&h| h != H3_NULL && h != sf_h3 && !ring.contains(&h))
      .unwrap();
    assert_eq!(
      direction_for_neighbor(sf_h3, *non_neighbor),
      Direction::InvalidDigit,
      "Non-neighbor should yield InvalidDigit"
    );

    // Test with self
    assert_eq!(
      direction_for_neighbor(sf_h3, sf_h3),
      Direction::Center,
      "Direction to self is Center"
    );
  }

  #[test]
  fn test_direction_for_neighbor_pentagon() {
    let pentagon_res2 = H3Index(0x820807fffffffff); // BC 4, all D=0 up to res 2
    assert!(h3_index::inspection::is_pentagon(pentagon_res2));

    let mut ring = [H3_NULL; 7];
    assert!(grid_disk(pentagon_res2, 1, &mut ring).is_ok());

    let mut found_count = 0;
    for i in 0..7 {
      if ring[i] == H3_NULL || ring[i] == pentagon_res2 {
        continue;
      }
      found_count += 1;
      let dir = direction_for_neighbor(pentagon_res2, ring[i]);
      assert_ne!(
        dir,
        Direction::InvalidDigit,
        "Should find direction for pentagon neighbor"
      );
      assert_ne!(dir, Direction::KAxes, "Direction from pentagon should not be KAxes");

      let mut rotations = 0;
      let mut recovered_neighbor = H3_NULL;
      assert!(h3_neighbor_rotations(pentagon_res2, dir, &mut rotations, &mut recovered_neighbor).is_ok());
      assert_eq!(recovered_neighbor, ring[i], "Pentagon direction round trip mismatch");
    }
    assert_eq!(found_count, 5, "Pentagon found 5 distinct neighbors");
  }

  #[test]
  fn test_are_neighbor_cells() {
    let mut sf_geo = LatLng::default();
    _set_geo_degs(&mut sf_geo, 37.779, -122.419);
    let sf_h3_res9 = lat_lng_to_cell(&sf_geo, 9).unwrap();

    let mut ring1 = [H3_NULL; 7];
    assert!(grid_disk(sf_h3_res9, 1, &mut ring1).is_ok());

    for i in 0..7 {
      if ring1[i] == H3_NULL {
        continue;
      }
      if ring1[i] == sf_h3_res9 {
        assert_eq!(
          are_neighbor_cells(sf_h3_res9, ring1[i]),
          Ok(false),
          "Cell is not neighbor with itself"
        );
      } else {
        assert_eq!(
          are_neighbor_cells(sf_h3_res9, ring1[i]),
          Ok(true),
          "Cell is neighbor with ring1 cell"
        );
      }
    }

    let mut ring2 = [H3_NULL; 19];
    assert!(grid_disk(sf_h3_res9, 2, &mut ring2).is_ok());
    for i in 0..19 {
      if ring2[i] == H3_NULL || ring2[i] == sf_h3_res9 || ring1.contains(&ring2[i]) {
        continue; // Skip self, nulls, and direct neighbors
      }
      assert_eq!(
        are_neighbor_cells(sf_h3_res9, ring2[i]),
        Ok(false),
        "Cell is not neighbor with ring2 cell"
      );
    }

    // Different resolutions
    let sf_h3_res8 = lat_lng_to_cell(&sf_geo, 8).unwrap();
    assert_eq!(are_neighbor_cells(sf_h3_res9, sf_h3_res8), Err(H3Error::ResMismatch));

    // Invalid input
    assert_eq!(are_neighbor_cells(sf_h3_res9, H3_NULL), Err(H3Error::CellInvalid));
  }
}
