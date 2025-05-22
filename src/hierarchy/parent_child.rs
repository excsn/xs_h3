// src/hierarchy/parent_child.rs

use crate::constants::MAX_H3_RES;
use crate::h3_index::inspection::is_valid_cell as h3_is_valid_cell;
use crate::h3_index::{
  _h3_leading_non_zero_digit, get_index_digit, get_resolution, is_pentagon, set_index_digit, set_resolution,
}; // Alias to avoid conflict if used locally

use crate::iterators::{iterInitParent, iterStepChild, IterCellsChildren}; // For cellToChildren
use crate::math::extensions::_ipow;
use crate::types::{Direction, H3Error, H3Index};

/// Helper: Zero out index digits from `start_res_inclusive` to `end_res_inclusive`.
/// No-op if `start_res_inclusive > end_res_inclusive`.
/// `start_res_inclusive` and `end_res_inclusive` are 1-based H3 resolution digit numbers.
pub(crate) fn _zero_index_digits(mut h: H3Index, start_res_inclusive: i32, end_res_inclusive: i32) -> H3Index {
  if start_res_inclusive > end_res_inclusive {
    return h;
  }
  // Iterate from the finest digit to be zeroed up to the coarsest.
  for r in start_res_inclusive..=end_res_inclusive {
    set_index_digit(&mut h, r, Direction::Center);
  }
  h
}

/// Returns whether one resolution is a valid child resolution for a cell.
/// Each resolution is considered a valid child resolution of itself.
fn _has_child_at_res(parent_h: H3Index, child_res: i32) -> bool {
  let parent_res = get_resolution(parent_h);
  child_res >= parent_res && child_res <= MAX_H3_RES
}

/// Produces the parent H3 index of `h` at `parent_res`.
pub fn cell_to_parent(h: H3Index, parent_res: i32) -> Result<H3Index, H3Error> {
  // println!(
  //   "---- cell_to_parent ---- START Input H3: {:x} (Res {}), Target Parent Res: {}",
  //   h.0,
  //   get_resolution(h), // Get resolution before any modification for logging
  //   parent_res
  // );

  if !h3_is_valid_cell(h) {
    // println!("  cell_to_parent: Input H3 {:x} is invalid.", h.0);
    return Err(H3Error::CellInvalid);
  }

  let child_res = get_resolution(h);
  if parent_res < 0 || parent_res > child_res || parent_res > MAX_H3_RES {
    // child_res check also implicitly done by is_valid_cell
    // println!(
    //   "  cell_to_parent: parent_res {} out of domain for child_res {}.",
    //   parent_res, child_res
    // );
    return Err(H3Error::ResDomain); // More appropriate for parent_res > child_res
  }

  if parent_res == child_res {
    // println!(
    //   "  cell_to_parent: Parent res same as child res. Returning original H3: {:x}",
    //   h.0
    // );
    return Ok(h);
  }

  let mut parent_h = h;
  // println!("  cell_to_parent: Initial parent_h (copy of input): {:x}", parent_h.0);

  set_resolution(&mut parent_h, parent_res);
  
  // println!(
  //   "  cell_to_parent: After set_resolution({}), parent_h: {:x}",
  //   parent_res, parent_h.0
  // );

  // Set digits from parent_res + 1 through child_res to InvalidDigit (7).
  // Digits from child_res + 1 to MAX_H3_RES were already 7 in the valid child `h`.
  for r_digit_to_invalidate in (parent_res + 1)..=child_res {
    // Get the digit *before* invalidating it for accurate logging
    let original_digit_at_r = get_index_digit(parent_h, r_digit_to_invalidate);
    // println!(
    //   "    cell_to_parent: Invalidating digit for resolution {} (original digit was {:?})",
    //   r_digit_to_invalidate, original_digit_at_r
    // );
    
    set_index_digit(&mut parent_h, r_digit_to_invalidate, Direction::InvalidDigit);
    
    // println!(
    //   "    cell_to_parent: After invalidating digit for res {}, parent_h is now: {:x}",
    //   r_digit_to_invalidate, parent_h.0
    // );
  }

  // println!("---- cell_to_parent ---- END Final parent_h: {:x}", parent_h.0);
  Ok(parent_h)
}

/// Determines the exact number of children (or grandchildren, etc.) for a given cell.
pub fn cell_to_children_size(h: H3Index, child_res: i32) -> Result<i64, H3Error> {
  if !h3_is_valid_cell(h) {
    return Err(H3Error::CellInvalid);
  }
  if !_has_child_at_res(h, child_res) {
    return Err(H3Error::ResDomain);
  }

  let n = child_res - get_resolution(h);

  if is_pentagon(h) {
    // Assuming is_pentagon is from h3_index::inspection
    // Formula for pentagon children: 1 (self) + 5 * (sum_{i=0}^{n-1} 7^i)
    // Sum is (7^n - 1) / (7 - 1) = (7^n - 1) / 6
    Ok(1 + 5 * (_ipow(7, n as i64) - 1) / 6)
  } else {
    Ok(_ipow(7, n as i64))
  }
}

/// Returns the center child of the given cell at the specified child resolution.
pub fn cell_to_center_child(h: H3Index, child_res: i32) -> Result<H3Index, H3Error> {
  if !h3_is_valid_cell(h) {
    return Err(H3Error::CellInvalid);
  }
  if !_has_child_at_res(h, child_res) {
    return Err(H3Error::ResDomain);
  }

  let parent_res = get_resolution(h);
  let mut child_h = h; // Start with a copy of the parent

  // Set to the child resolution
  set_resolution(&mut child_h, child_res);
  // Zero out all digit places from the parent's res + 1 up to the child's res.
  // This makes all finer digits "Center".
  child_h = _zero_index_digits(child_h, parent_res + 1, child_res);

  Ok(child_h)
}

/// Fills `children` with all H3 cells that are children of `h` at `child_res`.
pub fn cell_to_children(h: H3Index, child_res: i32, children: &mut [H3Index]) -> Result<(), H3Error> {
  let expected_size = cell_to_children_size(h, child_res)?; // Also validates h and child_res
  if children.len() < expected_size as usize {
    return Err(H3Error::MemoryBounds);
  }

  let mut iter = iterInitParent(h, child_res); // Assumes iterators.rs is ported
  let mut i: usize = 0;
  while iter.h != crate::types::H3_NULL {
    if i >= children.len() {
      // Should not happen if expected_size is correct
      return Err(H3Error::MemoryBounds);
    }
    children[i] = iter.h;
    i += 1;
    iterStepChild(&mut iter);
  }
  // Zero out remaining part of buffer if iter produced fewer than expected_size
  // This can happen if expected_size was for a hexagon but h was a pentagon.
  // `cellToChildrenSize` handles this, so this loop might not be strictly needed if size is exact.
  for k in i..(expected_size as usize) {
    // Fill up to actual expected_size
    if k < children.len() {
      children[k] = crate::types::H3_NULL;
    }
  }

  Ok(())
}

/// Validates a child position in the context of a given parent.
fn validate_child_pos(child_pos: i64, parent: H3Index, child_res: i32) -> Result<(), H3Error> {
  let max_child_count = cell_to_children_size(parent, child_res)?;
  if child_pos < 0 || child_pos >= max_child_count {
    return Err(H3Error::Domain); // childPos out of range
  }
  Ok(())
}

/// Returns the position of the child cell within an ordered list of all children
/// of the cell's parent at the specified `parent_res`.
///
/// # Arguments
/// * `child` - The H3 cell index of the child.
/// * `parent_res` - The resolution of the parent. Must be coarser than or equal to child's resolution.
///
/// # Returns
/// `Ok(position)` on success, or an `H3Error` if inputs are invalid.
pub fn cell_to_child_pos(child: H3Index, parent_res: i32) -> Result<i64, H3Error> {
  if !h3_is_valid_cell(child) {
    return Err(H3Error::CellInvalid);
  }
  let child_res = get_resolution(child);
  if parent_res < 0 || parent_res > child_res || parent_res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }

  let mut current_pos: i64 = 0;

  // This loop iterates from the coarsest digit level defining the child (parent_res + 1)
  // down to the finest digit level (child_res).
  for r_level_of_digit in (parent_res + 1)..=child_res {
    let digit_val = get_index_digit(child, r_level_of_digit);

    // Determine the nature of the parent of this specific digit.
    // This parent is at resolution (r_level_of_digit - 1).
    let parent_of_this_digit = cell_to_parent(child, r_level_of_digit - 1)?;
    let is_immediate_parent_pentagon = is_pentagon(parent_of_this_digit);

    // How many children at the final `child_res` are represented by a single slot at this `r_level_of_digit`?
    let children_per_slot_at_final_res = _ipow(7, (child_res - r_level_of_digit) as i64);

    let mut adjusted_digit_for_offset_calc = digit_val as i32;
    if is_immediate_parent_pentagon {
      if digit_val == Direction::KAxes {
        return Err(H3Error::CellInvalid);
      } // Should be caught by child's is_valid_cell
      if digit_val > Direction::KAxes {
        adjusted_digit_for_offset_calc -= 1;
      }
    }

    if adjusted_digit_for_offset_calc == 0 { // This digit is Center (or effective Center for pentagon path)
       // No contribution from *preceding* slots at this r_level_of_digit.
       // The position is within the "center slot" of its immediate parent.
    } else {
      // This digit is not Center. Add the count of all children from preceding slots.
      // Count of children in the center slot of the immediate parent:
      current_pos += if is_immediate_parent_pentagon {
        1 + 5 * (children_per_slot_at_final_res - 1) / 6
      } else {
        children_per_slot_at_final_res // A hex parent's center slot is just one hex-equivalent block
      };
      // Count of children in the (adjusted_digit_for_offset_calc - 1) preceding hex-equivalent slots:
      current_pos += (adjusted_digit_for_offset_calc as i64 - 1) * children_per_slot_at_final_res;
    }
  }
  Ok(current_pos)
}

/// Returns the child cell at a given `child_pos` within an ordered list of all
/// children of the `parent` cell, at the specified `child_res`.
///
/// # Arguments
/// * `child_pos` - The position in the ordered list (0-indexed).
/// * `parent` - The H3 cell index of the parent.
/// * `child_res` - The resolution of the desired child cell. Must be finer than or equal to parent's resolution.
///
/// # Returns
/// `Ok(H3Index)` of the child cell, or an `H3Error` if inputs are invalid.
pub fn child_pos_to_cell(child_pos: i64, parent: H3Index, child_res: i32) -> Result<H3Index, H3Error> {
  if !h3_is_valid_cell(parent) {
    return Err(H3Error::CellInvalid);
  }
  if child_res < 0 || child_res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }

  let parent_res = get_resolution(parent);
  if child_res < parent_res {
    return Err(H3Error::ResMismatch);
  }

  validate_child_pos(child_pos, parent, child_res)?; // Validates child_pos range

  let mut current_child_h = parent;
  set_resolution(&mut current_child_h, child_res); // Start with parent, set to child res (finer digits are 7s)

  let mut current_idx_remaining = child_pos;
  let is_original_parent_pentagon = is_pentagon(parent);

  if is_original_parent_pentagon {
    // Pentagon logic
    let mut current_parent_for_pent_logic = parent; // Tracks the parent at each sub-resolution step
    for r in (parent_res + 1)..=child_res {
      let is_current_parent_pentagon = is_pentagon(current_parent_for_pent_logic);
      let children_at_finer_res_count = _ipow(7, (child_res - r) as i64);

      if is_current_parent_pentagon {
        let pentagon_center_child_descendants = 1 + 5 * (children_at_finer_res_count - 1) / 6;
        if current_idx_remaining < pentagon_center_child_descendants {
          set_index_digit(&mut current_child_h, r, Direction::Center);
          // current_parent_for_pent_logic remains a pentagon for next iter
        } else {
          current_idx_remaining -= pentagon_center_child_descendants;
          let digit_val = (current_idx_remaining / children_at_finer_res_count) + 1; // +1 because digit 0 used by center
                                                                                     // This digit needs to skip K-axis (1)
          let actual_digit = if digit_val >= Direction::KAxes as i64 {
            Direction::try_from((digit_val + 1) as u8).unwrap_or(Direction::InvalidDigit)
          } else {
            Direction::try_from(digit_val as u8).unwrap_or(Direction::InvalidDigit)
          };
          set_index_digit(&mut current_child_h, r, actual_digit);
          // current_parent_for_pent_logic becomes a hexagon for next iter
          // We need to update it to reflect the child chosen at this step `r`.
          // Simplest: construct it based on current_child_h up to resolution `r`.
          let mut temp_next_parent = current_child_h;
          set_resolution(&mut temp_next_parent, r);
          temp_next_parent = _zero_index_digits(temp_next_parent, r + 1, child_res);
          current_parent_for_pent_logic = temp_next_parent;
        }
      } else {
        // Hexagon parent logic for this step
        let digit_val = current_idx_remaining / children_at_finer_res_count;
        set_index_digit(
          &mut current_child_h,
          r,
          Direction::try_from(digit_val as u8).unwrap_or(Direction::InvalidDigit),
        );
        // current_parent_for_pent_logic remains a hexagon
        // (or rather, its pentagon-ness doesn't change based on this step)
        if r < child_res {
          // Update for next iter if needed
          let mut temp_next_parent = current_child_h;
          set_resolution(&mut temp_next_parent, r);
          temp_next_parent = _zero_index_digits(temp_next_parent, r + 1, child_res);
          current_parent_for_pent_logic = temp_next_parent;
        }
      }
      current_idx_remaining %= children_at_finer_res_count;
    }
  } else {
    // Hexagon parent logic (simpler)
    for r in (parent_res + 1)..=child_res {
      let children_at_finer_res_count = _ipow(7, (child_res - r) as i64);
      let digit_val = current_idx_remaining / children_at_finer_res_count;
      set_index_digit(
        &mut current_child_h,
        r,
        Direction::try_from(digit_val as u8).unwrap_or(Direction::InvalidDigit),
      );
      current_idx_remaining %= children_at_finer_res_count;
    }
  }

  // The above logic should fill digits up to child_res. Digits beyond child_res should be InvalidDigit (7).
  // _zero_index_digits clears from parent_res+1 to child_res to Center initially if needed by center_child.
  // Here, we start with `parent` (which has 7s beyond its res) and set its res to `child_res` (still 7s).
  // Then we overwrite digits from parent_res+1 to child_res.
  // So the digits beyond child_res should remain 7.

  Ok(current_child_h)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::h3_index::_set_h3_index as set_h3_index_internal; // Renamed for clarity
  use crate::h3_index::get_resolution;
  use crate::types::{Direction, H3Index, H3_NULL};

  #[test]
  fn test_cell_to_parent() {
    let mut sf_geo = crate::types::LatLng::default();
    crate::latlng::_set_geo_degs(&mut sf_geo, 37.779, -122.419);
    let child_h = crate::indexing::lat_lng_to_cell(&sf_geo, 10).unwrap();

    let parent_h_res9 = cell_to_parent(child_h, 9).unwrap();
    assert_eq!(get_resolution(parent_h_res9), 9);
    assert_eq!(parent_h_res9.0, 0x089283082877ffff);

    let parent_h_res5 = cell_to_parent(child_h, 5).unwrap();
    assert_eq!(get_resolution(parent_h_res5), 5);
    assert_eq!(parent_h_res5.0, 0x085283083fffffff);

    assert_eq!(cell_to_parent(child_h, 10), Ok(child_h)); // Parent at same res is self
    assert_eq!(cell_to_parent(child_h, 11), Err(H3Error::ResDomain));
    assert_eq!(cell_to_parent(child_h, -1), Err(H3Error::ResDomain));
    assert_eq!(cell_to_parent(H3_NULL, 5), Err(H3Error::CellInvalid));
  }

  #[test]
  fn test_cell_to_children_size() {
    let mut h_hex = H3Index::default();
    set_h3_index_internal(&mut h_hex, 5, 10, Direction::Center); // Hexagon

    assert_eq!(cell_to_children_size(h_hex, 5), Ok(1));
    assert_eq!(cell_to_children_size(h_hex, 6), Ok(7));
    assert_eq!(cell_to_children_size(h_hex, 7), Ok(49));
    assert_eq!(cell_to_children_size(h_hex, 4), Err(H3Error::ResDomain));

    let mut h_pent = H3Index::default();
    set_h3_index_internal(&mut h_pent, 5, 4, Direction::Center); // Pentagon (BC 4)
    assert!(crate::h3_index::inspection::is_pentagon(h_pent));

    assert_eq!(cell_to_children_size(h_pent, 5), Ok(1));
    assert_eq!(cell_to_children_size(h_pent, 6), Ok(1 + 5 * (7 - 1) / 6)); // 1 + 5 = 6
    assert_eq!(cell_to_children_size(h_pent, 7), Ok(1 + 5 * (49 - 1) / 6)); // 1 + 5*8 = 41
  }

  #[test]
  fn test_cell_to_center_child() {
    let mut h_hex_res5 = H3Index::default();
    set_h3_index_internal(&mut h_hex_res5, 5, 10, Direction::IjAxes); // D5=IJ

    let center_child_res5 = cell_to_center_child(h_hex_res5, 5).unwrap();
    assert_eq!(center_child_res5, h_hex_res5, "Center child at same res is self");

    let center_child_res6 = cell_to_center_child(h_hex_res5, 6).unwrap();
    assert_eq!(get_resolution(center_child_res6), 6);
    assert_eq!(
      crate::h3_index::get_index_digit(center_child_res6, 6),
      Direction::Center
    );
    // Check that coarser digits match parent
    for r in 1..=5 {
      assert_eq!(
        crate::h3_index::get_index_digit(center_child_res6, r),
        crate::h3_index::get_index_digit(h_hex_res5, r)
      );
    }

    // Test with pentagon
    let mut h_pent_res2 = H3Index::default();
    set_h3_index_internal(&mut h_pent_res2, 2, 4, Direction::Center); // BC4, all center digits
    let center_child_pent_res4 = cell_to_center_child(h_pent_res2, 4).unwrap();
    assert_eq!(get_resolution(center_child_pent_res4), 4);
    assert_eq!(crate::h3_index::get_base_cell(center_child_pent_res4), 4);
    for r in 1..=4 {
      assert_eq!(
        crate::h3_index::get_index_digit(center_child_pent_res4, r),
        Direction::Center
      );
    }
    assert!(crate::h3_index::inspection::is_pentagon(center_child_pent_res4));
  }
  #[test]
  fn test_cell_to_child_pos_and_back() {
    let mut parent_hex = H3Index::default();
    crate::h3_index::_set_h3_index(&mut parent_hex, 2, 10, Direction::Center); // Res 2 hex

    let child_res = 4;
    let children_count = cell_to_children_size(parent_hex, child_res).unwrap();
    let mut children = vec![H3Index::default(); children_count as usize];
    cell_to_children(parent_hex, child_res, &mut children).unwrap();

    for (expected_pos, &child_h) in children.iter().enumerate() {
      if child_h == crate::types::H3_NULL {
        continue;
      }
      let pos = cell_to_child_pos(child_h, get_resolution(parent_hex)).unwrap();
      assert_eq!(
        pos, expected_pos as i64,
        "cellToChildPos mismatch for child {:x}",
        child_h.0
      );

      let recovered_child = child_pos_to_cell(pos, parent_hex, child_res).unwrap();
      assert_eq!(recovered_child, child_h, "childPosToCell mismatch for pos {}", pos);
    }
  }

  #[test]
  fn test_cell_to_child_pos_and_back_pentagon() {
    let mut parent_pent = H3Index::default();
    crate::h3_index::_set_h3_index(&mut parent_pent, 1, 4, Direction::Center); // Res 1 pentagon (BC 4)
    assert!(is_pentagon(parent_pent));

    let child_res = 3;
    let children_count = cell_to_children_size(parent_pent, child_res).unwrap();
    let mut children = vec![H3Index::default(); children_count as usize];
    cell_to_children(parent_pent, child_res, &mut children).unwrap();

    for (expected_pos, &child_h) in children.iter().enumerate() {
      if child_h == crate::types::H3_NULL {
        continue;
      } // Should not happen if list is compact
      assert!(h3_is_valid_cell(child_h));

      let pos = cell_to_child_pos(child_h, get_resolution(parent_pent)).unwrap();
      assert_eq!(
        pos, expected_pos as i64,
        "cellToChildPos mismatch for pentagon child {:x} (parent {:x}, childRes {}). Expected pos {}, got {}",
        child_h.0, parent_pent.0, child_res, expected_pos, pos
      );

      let recovered_child = child_pos_to_cell(pos, parent_pent, child_res).unwrap();
      assert_eq!(
        recovered_child, child_h,
        "childPosToCell mismatch for pentagon parent {:x} at pos {}. Expected {:x}, got {:x}",
        parent_pent.0, pos, child_h.0, recovered_child.0
      );
    }
  }

  #[test]
  fn test_cell_to_child_pos_errors() {
    let child = crate::indexing::lat_lng_to_cell(&crate::types::LatLng { lat: 0.0, lng: 0.0 }, 8).unwrap();
    assert_eq!(cell_to_child_pos(child, -1), Err(H3Error::ResDomain));
    assert_eq!(cell_to_child_pos(child, 16), Err(H3Error::ResDomain));
    assert_eq!(cell_to_child_pos(child, 9), Err(H3Error::ResDomain)); // Parent res finer
  }

  #[test]
  fn test_child_pos_to_cell_errors() {
    let parent = crate::indexing::lat_lng_to_cell(&crate::types::LatLng { lat: 0.0, lng: 0.0 }, 5).unwrap();
    assert_eq!(child_pos_to_cell(0, parent, 4), Err(H3Error::ResMismatch)); // Child res coarser
    assert_eq!(child_pos_to_cell(0, parent, 16), Err(H3Error::ResDomain));
    assert_eq!(child_pos_to_cell(-1, parent, 6), Err(H3Error::Domain)); // Negative pos
    assert_eq!(child_pos_to_cell(100, parent, 6), Err(H3Error::Domain)); // Pos too large (max 7 for k=1)
  }
}
