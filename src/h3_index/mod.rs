#![allow(clippy::cast_possible_truncation)] // For H3_GET_INDEX_DIGIT style operations

pub mod inspection;
pub mod string_conv;

use crate::base_cells::{
  _base_cell_is_cw_offset, _face_ijk_to_base_cell, _face_ijk_to_base_cell_ccwrot60, _is_base_cell_pentagon,
  BASE_CELL_DATA, INVALID_BASE_CELL, INVALID_ROTATIONS,
};
use crate::coords::face_ijk::{Overage, _adjust_overage_class_ii};
use crate::coords::ijk::{
  _down_ap7, _down_ap7r, _ijk_normalize, _ijk_sub, _neighbor, _rotate60_ccw, _rotate60_cw, _unit_ijk_to_digit, _up_ap7,
  _up_ap7r, UNIT_VECS,
};
use crate::types::{CoordIJK, Direction, FaceIJK, H3Index};
use crate::{constants::*, H3_NULL}; // From base_cells.rs

pub use inspection::{is_pentagon, get_num_cells};
pub use string_conv::{string_to_h3, h3_to_string, h3_to_string_alloc};

// H3Index bit layout accessors/mutators (macros from C translated to inline functions)

/// Gets the mode of the H3 index.
#[inline(always)]
#[must_use]
pub const fn get_mode(h: H3Index) -> u8 {
  ((h.0 & H3_MODE_MASK) >> H3_MODE_OFFSET) as u8
}

/// Sets the mode of the H3 index.
#[inline(always)]
pub fn set_mode(h: &mut H3Index, mode: u8) {
  h.0 = (h.0 & H3_MODE_MASK_NEGATIVE) | ((mode as u64) << H3_MODE_OFFSET);
}

/// Gets the resolution of the H3 index.
#[inline(always)]
#[must_use]
pub const fn get_resolution(h: H3Index) -> i32 {
  ((h.0 & H3_RES_MASK) >> H3_RES_OFFSET) as i32
}

/// Sets the resolution of the H3 index.
#[inline(always)]
pub fn set_resolution(h: &mut H3Index, res: i32) {
  h.0 = (h.0 & H3_RES_MASK_NEGATIVE) | ((res as u64) << H3_RES_OFFSET);
}

/// Gets the base cell of the H3 index.
#[inline(always)]
#[must_use]
pub const fn get_base_cell(h: H3Index) -> i32 {
  ((h.0 & H3_BC_MASK) >> H3_BC_OFFSET) as i32
}

/// Sets the base cell of the H3 index.
#[inline(always)]
pub fn set_base_cell(h: &mut H3Index, bc: i32) {
  h.0 = (h.0 & H3_BC_MASK_NEGATIVE) | ((bc as u64) << H3_BC_OFFSET);
}

/// Gets the H3 digit at the given resolution `res` from the H3 index.
/// `res` must be between 1 and `get_resolution(h)`.
#[inline(always)]
#[must_use]
pub fn get_index_digit(h: H3Index, res: i32) -> Direction {
  // res is 1-based for digits
  let val = ((h.0 >> ((MAX_H3_RES - res) * H3_PER_DIGIT_OFFSET as i32)) & H3_DIGIT_MASK) as u8;
  // This conversion is safe because H3_DIGIT_MASK ensures val is 0-7.
  unsafe { std::mem::transmute(val) }
}

/// Sets the H3 digit at the given resolution `res` in the H3 index.
/// `res` must be between 1 and `get_resolution(h)`.
#[inline(always)]
pub fn set_index_digit(h: &mut H3Index, res: i32, digit: Direction) {
  // res is 1-based for digits
  let offset = (MAX_H3_RES - res) * H3_PER_DIGIT_OFFSET as i32;
  h.0 = (h.0 & !(H3_DIGIT_MASK << offset)) | ((digit as u64) << offset);
}

/// Gets the reserved bits of the H3 index. Should be 0 for valid cell indexes.
#[inline(always)]
#[must_use]
pub const fn get_reserved_bits(h: H3Index) -> u8 {
  ((h.0 & H3_RESERVED_MASK) >> H3_RESERVED_OFFSET) as u8
}

/// Sets the reserved bits of the H3 index.
#[inline(always)]
pub fn set_reserved_bits(h: &mut H3Index, v: u8) {
  h.0 = (h.0 & H3_RESERVED_MASK_NEGATIVE) | ((v as u64) << H3_RESERVED_OFFSET);
}

/// Gets the high bit of the H3 index (should be 0).
#[inline(always)]
#[must_use]
pub const fn get_high_bit(h: H3Index) -> u8 {
  ((h.0 & H3_HIGH_BIT_MASK) >> 63) as u8 // H3_MAX_OFFSET is 63 in C for the shift amount
}

/// Sets the high bit of the H3 index.
#[inline(always)]
pub fn set_high_bit(h: &mut H3Index, v: u8) {
  h.0 = (h.0 & H3_HIGH_BIT_MASK_NEGATIVE) | ((v as u64) << 63);
}

/// Initializes an H3 index with a given resolution, base cell, and initial digit for all resolution levels.
/// The mode is set to H3_CELL_MODE.
///
/// # Arguments
/// * `h` - Output: The H3 index to initialize.
/// * `res` - The H3 resolution.
/// * `base_cell` - The H3 base cell number.
/// * `init_digit` - The H3 `Direction` to initialize all resolution digits to.
pub(crate) fn _set_h3_index(h: &mut H3Index, res: i32, base_cell: i32, init_digit: Direction) {
  h.0 = H3_INIT; // Start with the base pattern (all digits 7, mode 0, res 0, bc 0)
  set_mode(h, H3_CELL_MODE);
  set_resolution(h, res);
  set_base_cell(h, base_cell);
  for r in 1..=res {
    set_index_digit(h, r, init_digit);
  }
}

/// Returns whether or not a resolution is a Class III grid.
/// Odd resolutions are Class III, even are Class II.
#[inline]
#[must_use]
pub(crate) fn is_resolution_class_iii(r: i32) -> bool {
  // r % 2 != 0 would also work
  (r % 2) == 1
}

/// Returns the highest resolution non-zero digit in an H3Index.
#[inline]
#[must_use]
pub(crate) fn _h3_leading_non_zero_digit(h: H3Index) -> Direction {
  let res = get_resolution(h);
  for r in 1..=res {
    let digit = get_index_digit(h, r);
    if digit != Direction::Center {
      return digit;
    }
  }
  Direction::Center // If all are 0 (Center), then Center is the leading non-zero (or only) digit.
}

/// Rotate an H3Index 60 degrees counter-clockwise.
#[inline]
pub(crate) fn _h3_rotate60_ccw(mut h: H3Index) -> H3Index {
  let res = get_resolution(h);
  for r in 1..=res {
    let old_digit = get_index_digit(h, r); // Read before mutable borrow
    let new_digit = crate::coords::ijk::_rotate60_ccw(old_digit); // Assuming _rotate60_ccw for Direction is in ijk module
    set_index_digit(&mut h, r, new_digit);
  }
  h
}

/// Rotate an H3Index 60 degrees clockwise.
#[inline]
pub(crate) fn _h3_rotate60_cw(mut h: H3Index) -> H3Index {
  let res = get_resolution(h);
  for r in 1..=res {
    let old_digit = get_index_digit(h, r); // Read before mutable borrow
    let new_digit = crate::coords::ijk::_rotate60_cw(old_digit); // Assuming _rotate60_cw for Direction is in ijk module
    set_index_digit(&mut h, r, new_digit);
  }
  h
}
/// Rotate an H3Index 60 degrees counter-clockwise about a pentagonal center.
#[inline]
pub(crate) fn _h3_rotate_pent60_ccw(mut h: H3Index) -> H3Index {
  let res = get_resolution(h);
  let mut found_first_non_zero_digit = false;
  for r in 1..=res {
    let old_digit = get_index_digit(h, r); // Read before mutable borrow
    let rotated_digit = crate::coords::ijk::_rotate60_ccw(old_digit);
    set_index_digit(&mut h, r, rotated_digit);

    // After set_index_digit, h has been modified.
    // We need to re-read the digit if _h3_leading_non_zero_digit depends on the current h.
    if !found_first_non_zero_digit && get_index_digit(h, r) != Direction::Center {
      found_first_non_zero_digit = true;
      // _h3_leading_non_zero_digit reads the current state of h
      if _h3_leading_non_zero_digit(h) == Direction::KAxes {
        h = _h3_rotate60_ccw(h); // Recursive call, takes a new mutable copy
      }
    }
  }
  h
}

/// Rotate an H3Index 60 degrees clockwise about a pentagonal center.
#[inline]
pub(crate) fn _h3_rotate_pent60_cw(mut h: H3Index) -> H3Index {
  let res = get_resolution(h);
  let mut found_first_non_zero_digit = false;
  for r in 1..=res {
    let old_digit = get_index_digit(h, r); // Read before mutable borrow
    let rotated_digit = crate::coords::ijk::_rotate60_cw(old_digit);
    set_index_digit(&mut h, r, rotated_digit);

    if !found_first_non_zero_digit && get_index_digit(h, r) != Direction::Center {
      found_first_non_zero_digit = true;
      if _h3_leading_non_zero_digit(h) == Direction::KAxes {
        h = _h3_rotate60_cw(h); // Recursive call
      }
    }
  }
  h
}

/// Convert a `FaceIJK` address to the corresponding `H3Index`.
///
/// # Arguments
/// * `fijk` - The `FaceIJK` address.
/// * `res` - The cell resolution.
///
/// # Returns
/// The encoded `H3Index`, or `H3_NULL` on failure (e.g., invalid input).
// This is `_faceIjkToH3` in C
pub(crate) fn _face_ijk_to_h3(fijk: &FaceIJK, res: i32) -> H3Index {
  // initialize the index
  let mut h = H3Index(H3_INIT); // Start with the base pattern
  set_mode(&mut h, H3_CELL_MODE as u8); // Set mode to Cell
  set_resolution(&mut h, res);

  // check for res 0/base cell
  if res == 0 {
    // MAX_FACE_COORD is 2 in C, implies max valid IJK component for res 0 lookup
    if fijk.coord.i > 2 || fijk.coord.j > 2 || fijk.coord.k > 2 {
      return H3_NULL; // out of range input
    }
    let base_cell = _face_ijk_to_base_cell(fijk);
    if base_cell == INVALID_BASE_CELL {
      // Check if lookup failed
      return H3_NULL;
    }
    set_base_cell(&mut h, base_cell);
    return h;
  }

  // We need to find the correct base cell FaceIJK for this H3 index.
  // Start with the passed in face and resolution `res` IJK coordinates
  // in that face's coordinate system.
  let mut fijk_bc = *fijk; // Make a mutable copy to find the base cell Fijk

  // Build the H3Index from finest res up to coarsest (res 1).
  // The loop iterates from `res-1` down to `0`.
  // `r_plus_1` is the "current digit resolution" being determined (from 1 to `res`).
  let ijk_mut_for_bc_finding = &mut fijk_bc.coord; // Mutable borrow for finding BC's Fijk

  for r_iter in (0..res).rev() {
    let r_plus_1 = r_iter + 1; // Current digit's resolution (1 to res)
    let last_ijk_for_digit = *ijk_mut_for_bc_finding; // Store current IJK before up-aperture

    let mut last_center_of_parent = CoordIJK::default();

    if is_resolution_class_iii(r_plus_1) {
      // Class III parent == up_ap7 (counter-clockwise)
      _up_ap7(ijk_mut_for_bc_finding); // Modifies ijk_mut_for_bc_finding to be parent in IJK
      last_center_of_parent = *ijk_mut_for_bc_finding; // This is now the parent's IJK
      _down_ap7(&mut last_center_of_parent); // Find center of this parent's children set
    } else {
      // Class II parent == up_ap7r (clockwise)
      _up_ap7r(ijk_mut_for_bc_finding);
      last_center_of_parent = *ijk_mut_for_bc_finding;
      _down_ap7r(&mut last_center_of_parent);
    }

    let mut diff = CoordIJK::default();
    _ijk_sub(&last_ijk_for_digit, &last_center_of_parent, &mut diff);
    _ijk_normalize(&mut diff); // diff is now the unit vector for the digit

    set_index_digit(&mut h, r_plus_1, _unit_ijk_to_digit(&diff));
  }

  // fijk_bc.coord should now hold the IJK of the base cell within fijk_bc.face's system.
  // MAX_FACE_COORD is 2
  if fijk_bc.coord.i > 2 || fijk_bc.coord.j > 2 || fijk_bc.coord.k > 2 {
    return H3_NULL; // Base cell IJK out of range for face-based lookup
  }

  // Lookup the correct base cell ID and its canonical rotation on this face
  let base_cell = _face_ijk_to_base_cell(&fijk_bc);
  if base_cell == INVALID_BASE_CELL {
    // Check if lookup failed
    return H3_NULL;
  }
  set_base_cell(&mut h, base_cell);

  let num_rots = _face_ijk_to_base_cell_ccwrot60(&fijk_bc);
  if num_rots == INVALID_ROTATIONS {
    // Check if rotation lookup failed
    return H3_NULL;
  }

  // Rotate if necessary to get canonical base cell orientation.
  if _is_base_cell_pentagon(base_cell) {
    // Force rotation out of missing k-axes sub-sequence if necessary
    if _h3_leading_non_zero_digit(h) == Direction::KAxes {
      // This K_AXES_DIGIT check is for the H3 index `h` itself, not fijk_bc.coord
      if _base_cell_is_cw_offset(base_cell, fijk_bc.face) {
        h = _h3_rotate_pent60_cw(h);
      } else {
        h = _h3_rotate_pent60_ccw(h);
      }
    }

    for _i in 0..num_rots {
      h = _h3_rotate_pent60_ccw(h);
    }
  } else {
    // Hexagon base cell
    for _i in 0..num_rots {
      h = _h3_rotate60_ccw(h);
    }
  }
  h
}

/// Convert an `H3Index` to a `FaceIJK` address.
///
/// # Arguments
/// * `h` - The `H3Index`.
/// * `fijk` - Output: The corresponding `FaceIJK` address.
///
/// # Returns
/// `Ok(())` on success, or an `H3Error` if `h` is invalid.
/// Convert an `H3Index` to its canonical `FaceIJK` address.
///
/// # Arguments
/// * `h` - The `H3Index`.
/// * `fijk` - Output: The corresponding `FaceIJK` address.
///
/// # Returns
/// `Ok(())` on success, or an `H3Error` if `h` is invalid.
pub(crate) fn _h3_to_face_ijk(h: H3Index, fijk: &mut FaceIJK) -> Result<(), crate::types::H3Error> {
  let base_cell = get_base_cell(h);
  if base_cell < 0 || base_cell >= (NUM_BASE_CELLS as i32) {
    *fijk = FaceIJK::default(); // Zero out for safety on error
    return Err(crate::types::H3Error::CellInvalid);
  }

  // Start with the "home" face and IJK+ coordinates for the base cell of h.
  *fijk = BASE_CELL_DATA[base_cell as usize].home_fijk;

  // If h's resolution is 0, we're done (fijk is already correct).
  let res_h = get_resolution(h);
  if res_h == 0 {
    return Ok(());
  }

  // The H3 digits are in the canonical orientation. The base cell's home Fijk might
  // be in a different orientation on the face. We need to adjust `h`'s digits
  // to be relative to the base cell's home Fijk orientation *before* applying them.
  // The C code achieves this by rotating `h` based on `numRots` from `_faceIjkToBaseCellCCWrot60`.
  // However, `h` is *already* canonical. The `fijk.coord` from `baseCellData` is the one
  // that needs to align with how `h`'s digits will be applied.
  // The C logic for `_h3ToFaceIjk` first gets the base cell data.
  // Then it applies pentagon leading digit 5 rotation to `h`.
  // Then it calls `_h3ToFaceIjkWithInitializedFijk`.
  // Then it does overage adjustment.

  let mut h_for_digits = h; // Use a copy for potential digit rotation
  if _is_base_cell_pentagon(base_cell) && _h3_leading_non_zero_digit(h_for_digits) == Direction::IkAxes {
    // IK_AXES_DIGIT is 5
    h_for_digits = _h3_rotate60_cw(h_for_digits);
  }

  if !_h3_to_face_ijk_with_initialized_fijk(h_for_digits, fijk) {
    // No overage was possible from the initial face, so we're done.
    return Ok(());
  }

  // If we're here, there was a possibility of overage.
  // `fijk` now contains the coordinates on the initial home face.
  // We need to check and adjust for overage.
  let orig_ijk_on_home_face = fijk.coord; // Save before Class III adjustment for overage

  let mut res_for_overage = res_h;
  if is_resolution_class_iii(res_h) {
    // Overage is always Class II logic. Temporarily use Class II coords.
    _down_ap7r(&mut fijk.coord);
    res_for_overage += 1;
  }

  let pent_leading_4 =
    _is_base_cell_pentagon(base_cell) && _h3_leading_non_zero_digit(h_for_digits) == Direction::IAxes; // I_AXES_DIGIT is 4

  if _adjust_overage_class_ii(fijk, res_for_overage, pent_leading_4, false) != Overage::NoOverage {
    // Overage occurred.
    // If the base cell is a pentagon, we have the potential for secondary overages.
    if _is_base_cell_pentagon(base_cell) {
      while _adjust_overage_class_ii(fijk, res_for_overage, false, false) != Overage::NoOverage {
        // Loop until no more overages on new faces for this pentagonal case
      }
    }
    // If we had adjusted resolution for Class III overage check, convert IJK back to original res.
    if res_for_overage != res_h {
      _up_ap7r(&mut fijk.coord); // Convert Class II IJK back to Class III parent
    }
  } else {
    // No overage occurred.
    // If we had adjusted resolution for Class III, revert to original IJK on home face.
    if res_for_overage != res_h {
      fijk.coord = orig_ijk_on_home_face;
    }
  }
  Ok(())
}

/// Convert an H3Index to a `FaceIJK` address on a specified icosahedral face.
///
/// # Arguments
/// * `h` - The `H3Index`.
/// * `fijk` - Output: The `FaceIJK` address, initialized with the desired face
///   and normalized base cell coordinates from `h`.
///
/// # Returns
/// `true` if the possibility of overage exists (meaning the cell might not be
/// entirely on the initial face), `false` otherwise.
/// This function itself doesn't return an error, but an invalid `h` might lead to
/// unexpected `fijk` output. It's assumed `h` is valid.
pub(crate) fn _h3_to_face_ijk_with_initialized_fijk(h: H3Index, fijk: &mut FaceIJK) -> bool {
  let ijk_mut = &mut fijk.coord; // mutable borrow for in-place modifications
  let res_h = get_resolution(h);

  // Start with the fijk specific to the base cell of h.
  // The fijk argument is already initialized with the "home" face and IJK of the base cell.
  // So, fijk.coord here is effectively the base cell's center in its home face's system.
  // We now need to apply the digits of h to traverse down to the cell's resolution.

  // Center base cell hierarchy is entirely on this face,
  // unless the base cell is a pentagon, or resolution is 0.
  let mut possible_overage = true; // Assume overage is possible
  let base_cell = get_base_cell(h);

  if !_is_base_cell_pentagon(base_cell) && (res_h == 0 || (ijk_mut.i == 0 && ijk_mut.j == 0 && ijk_mut.k == 0)) {
    possible_overage = false; // No overage possible if on center of non-pent face, or res 0.
  }

  for r in 1..=res_h {
    if is_resolution_class_iii(r) {
      // Class III == rotate ccw
      _down_ap7(ijk_mut);
    } else {
      // Class II == rotate cw
      _down_ap7r(ijk_mut);
    }
    _neighbor(ijk_mut, get_index_digit(h, r));
  }

  possible_overage
}

// We need to add definitions for `BASE_CELL_DATA` and other base cell helpers
// in `src/base_cells.rs` and import them.
// For now, to allow this to compile, I'll put dummy definitions here or assume they exist.
// You should create `src/base_cells.rs` next.

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{coords::ijk::_ijk_matches, types::H3_NULL};

  #[test]
  fn test_get_set_mode() {
    let mut h = H3Index(0);
    for mode_val in 0..=15u8 {
      set_mode(&mut h, mode_val);
      assert_eq!(get_mode(h), mode_val, "Mode set/get mismatch for {}", mode_val);
    }
  }

  #[test]
  fn test_get_set_resolution() {
    let mut h = H3Index(0);
    for res_val in 0..=MAX_H3_RES {
      set_resolution(&mut h, res_val);
      assert_eq!(
        get_resolution(h),
        res_val,
        "Resolution set/get mismatch for {}",
        res_val
      );
    }
  }

  #[test]
  fn test_get_set_base_cell() {
    let mut h = H3Index(0);
    // NUM_BASE_CELLS is 122. Max value for 7 bits is 127.
    for bc_val in 0..122i32 {
      // Only valid base cells
      set_base_cell(&mut h, bc_val);
      assert_eq!(get_base_cell(h), bc_val, "Base cell set/get mismatch for {}", bc_val);
    }
  }

  #[test]
  fn test_get_set_index_digit() {
    let mut h = H3Index(0);
    set_resolution(&mut h, MAX_H3_RES); // Set max resolution to use all digit slots

    for res_level in 1..=MAX_H3_RES {
      for digit_val_u8 in 0..=6u8 {
        // Valid digits are 0-6 for this test
        let digit_val: Direction = unsafe { std::mem::transmute(digit_val_u8) };
        set_index_digit(&mut h, res_level, digit_val);
        assert_eq!(
          get_index_digit(h, res_level),
          digit_val,
          "Digit set/get mismatch for res {}, digit {}",
          res_level,
          digit_val_u8
        );
      }
    }
  }

  #[test]
  fn test_get_set_reserved_bits() {
    let mut h = H3Index(0);
    for val in 0..=0b111u8 {
      // 3 reserved bits
      set_reserved_bits(&mut h, val);
      assert_eq!(get_reserved_bits(h), val, "Reserved bits set/get mismatch for {}", val);
    }
  }

  #[test]
  fn test_get_set_high_bit() {
    let mut h = H3Index(0);
    set_high_bit(&mut h, 1);
    assert_eq!(get_high_bit(h), 1, "High bit set to 1");
    set_high_bit(&mut h, 0);
    assert_eq!(get_high_bit(h), 0, "High bit set to 0");
  }

  #[test]
  fn test_set_h3_index_function() {
    // Renamed from _set_h3_index to avoid conflict
    let mut h = H3Index::default();
    _set_h3_index(&mut h, 5, 12, Direction::KAxes); // KAxes is 1
    assert_eq!(get_resolution(h), 5, "resolution for _set_h3_index");
    assert_eq!(get_base_cell(h), 12, "base cell for _set_h3_index");
    assert_eq!(get_mode(h), H3_CELL_MODE, "mode for _set_h3_index");
    for r in 1..=5 {
      assert_eq!(get_index_digit(h, r), Direction::KAxes, "digit for _set_h3_index");
    }
    for r in 6..=MAX_H3_RES {
      // Check that higher res digits are 7 (Invalid)
      assert_eq!(
        get_index_digit(h, r),
        Direction::InvalidDigit,
        "blanked digit for _set_h3_index"
      );
    }
    // C literal 0x85184927fffffffL corresponds to:
    // Mode 1, Res 5, BC 12 (0x0C), Digits 1,1,1,1,1 then 7s
    // 0b_0_0001_000_0101_0001100_001_001_001_001_001_111_111_111_111_111_111_111_111_111_111
    //    H M    Rsv Res  BaseCell D1  D2  D3  D4  D5  D6  D7  D8  D9 D10 D11 D12 D13 D14 D15
    // Expected: Mode=1, Res=5, BC=12. Digits 1-5 are 1 (KAxes). Digits 6-15 are 7 (Invalid).
    // H3_INIT part (digits 7): 0x00001FFFFFFFFFFF
    // Digits 1-5 as 1 (001):
    // D1: 001 << (3*14)
    // D2: 001 << (3*13)
    // D3: 001 << (3*12)
    // D4: 001 << (3*11)
    // D5: 001 << (3*10)
    // This sums to 0x0000049249249240 (if only these digits were set and others 0)
    // With H3_INIT as base, we clear digits 1-5, then OR these.
    // Mask to clear D1-D5 from H3_INIT (which has all 7s):
    // 0x00000000000FFFFF (digits 6-15 as 7s)
    //   | ((0b001001001001001u64) << 30) // D1-D5 as 1s
    // Base value for mode, res, bc:
    // (1<<59) | (5<<52) | (12<<45) = 0x8518000000000000
    // Digits part (1-5 are 1, 6-15 are 7):
    // ( ( (1u64<<15)-1 ) << (3*10) ) ^ ( ( (1u64<<(3*5))-1 ) << (3*10) ) // set 6-15 to 7 (0x1FFFFC0000000)
    //   | ( ( (1u64<<(3*5)) / 7 ) * 1 ) // set 1-5 to 1 (0x0000049249240) -> this is too complex
    // A simpler construction:
    // Start with 0. Set mode, res, bc. Then loop set_index_digit.
    let mut expected_h_val = 0u64;
    expected_h_val |= (H3_CELL_MODE as u64) << H3_MODE_OFFSET;
    expected_h_val |= (5u64) << H3_RES_OFFSET;
    expected_h_val |= (12u64) << H3_BC_OFFSET;
    for r in 1..=5 {
      expected_h_val |= (Direction::KAxes as u64) << ((MAX_H3_RES - r) * H3_PER_DIGIT_OFFSET as i32);
    }
    for r in 6..=MAX_H3_RES {
      expected_h_val |= (Direction::InvalidDigit as u64) << ((MAX_H3_RES - r) * H3_PER_DIGIT_OFFSET as i32);
    }
    assert_eq!(
      h.0, expected_h_val,
      "index matches reconstructed literal. Got {:x}, Expected {:x}",
      h.0, expected_h_val
    );
    // The C literal 0x85184927fffffffL is correct for this setup.
    assert_eq!(h.0, 0x85184927fffffff_u64, "index matches C literal");
  }

  #[test]
  fn test_is_resolution_class_iii_test() {
    // Renamed from is_resolution_class_iii
    assert!(!is_resolution_class_iii(0), "Res 0 is Class II");
    assert!(is_resolution_class_iii(1), "Res 1 is Class III");
    assert!(!is_resolution_class_iii(2), "Res 2 is Class II");
  }

  #[test]
  fn test_h3_leading_non_zero_digit() {
    let mut h = H3Index::default();
    _set_h3_index(&mut h, 5, 0, Direction::Center); // All center digits
    assert_eq!(_h3_leading_non_zero_digit(h), Direction::Center);

    set_index_digit(&mut h, 3, Direction::JAxes); // Set a non-zero digit
    assert_eq!(_h3_leading_non_zero_digit(h), Direction::JAxes);

    set_index_digit(&mut h, 1, Direction::KAxes); // Set an earlier non-zero digit
    assert_eq!(_h3_leading_non_zero_digit(h), Direction::KAxes);
  }

  #[test]
  fn test_h3_rotations() {
    let mut h_i = H3Index::default(); // Represents I-axis direction at some res
    _set_h3_index(&mut h_i, 1, 0, Direction::IAxes);

    let mut h_j = H3Index::default();
    _set_h3_index(&mut h_j, 1, 0, Direction::JAxes);

    let mut h_k = H3Index::default();
    _set_h3_index(&mut h_k, 1, 0, Direction::KAxes);

    let mut h_ij = H3Index::default();
    _set_h3_index(&mut h_ij, 1, 0, Direction::IjAxes);

    let mut h_jk = H3Index::default();
    _set_h3_index(&mut h_jk, 1, 0, Direction::JkAxes);

    let mut h_ik = H3Index::default();
    _set_h3_index(&mut h_ik, 1, 0, Direction::IkAxes);

    assert_eq!(_h3_rotate60_ccw(h_i), h_ij, "I ccw -> IJ");
    assert_eq!(_h3_rotate60_cw(h_i), h_ik, "I cw -> IK");

    // Pentagonal rotations (should behave same as hexagonal if not on actual pentagon bc)
    assert_eq!(_h3_rotate_pent60_ccw(h_i), h_ij, "Pent I ccw -> IJ");
    assert_eq!(_h3_rotate_pent60_cw(h_i), h_ik, "Pent I cw -> IK");

    // Test a pentagon base cell that's not polar, e.g. BC 14
    // Its K digit (1) is deleted. So leading non-zero won't be K.
    let mut h_pent_center = H3Index::default();
    _set_h3_index(&mut h_pent_center, 1, 14, Direction::Center); // Res 1 pentagon center child

    let mut h_pent_j = H3Index::default();
    _set_h3_index(&mut h_pent_j, 1, 14, Direction::JAxes); // Res 1 pentagon J child

    // Rotate pentagon J child CCW. _h3_leading_non_zero_digit is J (2). Not K (1). No special rotation.
    // J (2) ccw -> JK (3).
    let mut h_pent_j_ccw_expected = H3Index::default();
    _set_h3_index(&mut h_pent_j_ccw_expected, 1, 14, Direction::JkAxes);
    assert_eq!(
      _h3_rotate_pent60_ccw(h_pent_j),
      h_pent_j_ccw_expected,
      "Pent J child ccw"
    );
  }

  #[test]
  fn test_face_ijk_h3_roundtrip_res0_simple() {
    // Test a few specific FaceIJKs at res 0
    let fijk_orig_on_face0_center = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 0, j: 0, k: 0 },
    };
    let h3_on_face0_center = _face_ijk_to_h3(&fijk_orig_on_face0_center, 0);
    assert_ne!(h3_on_face0_center, H3_NULL, "face 0 center should be valid");
    let mut fijk_rt_face0_center = FaceIJK::default();
    assert!(_h3_to_face_ijk(h3_on_face0_center, &mut fijk_rt_face0_center).is_ok());
    // Base cell 16's home is Face 0, IJK {0,0,0}
    assert_eq!(fijk_rt_face0_center.face, 0);
    assert!(_ijk_matches(
      &fijk_rt_face0_center.coord,
      &CoordIJK { i: 0, j: 0, k: 0 }
    ));

    // Base cell 4 (pentagon) home is Face 0, IJK {2,0,0}
    let fijk_orig_bc4_home = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 2, j: 0, k: 0 },
    };
    let h3_bc4 = _face_ijk_to_h3(&fijk_orig_bc4_home, 0);
    assert_ne!(h3_bc4, H3_NULL, "BC4 home Fijk should be valid");
    assert_eq!(get_base_cell(h3_bc4), 4, "BC4 H3Index has base cell 4");
    let mut fijk_rt_bc4 = FaceIJK::default();
    assert!(_h3_to_face_ijk(h3_bc4, &mut fijk_rt_bc4).is_ok());
    assert_eq!(fijk_rt_bc4.face, fijk_orig_bc4_home.face); // Should be on home face
    assert!(_ijk_matches(&fijk_rt_bc4.coord, &fijk_orig_bc4_home.coord));
  }

  #[test]
  fn test_is_resolution_class_iii() {
    assert!(!is_resolution_class_iii(0), "Res 0 is Class II");
    assert!(is_resolution_class_iii(1), "Res 1 is Class III");
    assert!(!is_resolution_class_iii(2), "Res 2 is Class II");
  }

  // Tests for _face_ijk_to_h3 and _h3_to_face_ijk roundtrip
  #[test]
  fn test_face_ijk_h3_roundtrip_res0() {
    for face_idx in 0..NUM_ICOSA_FACES {
      for i in 0..=2 {
        // Max face coord for res 0
        for j in 0..=2 {
          for k in 0..=2 {
            let mut fijk_orig = FaceIJK {
              face: face_idx,
              coord: CoordIJK { i, j, k },
            };
            // Skip invalid combinations for res 0 based on C's _faceIjkToBaseCell logic
            if i + j + k > 2 && (i > 0 && j > 0 && k > 0) {
              continue;
            } // Simplified check, more precise one might be needed

            _ijk_normalize(&mut fijk_orig.coord); // Ensure coord is H3-normalized for fair comparison later

            let h3_index = _face_ijk_to_h3(&fijk_orig, 0);
            if h3_index == H3_NULL {
              continue;
            } // Skip if input fijk was invalid for res 0

            let mut fijk_roundtrip = FaceIJK::default();
            let result = _h3_to_face_ijk(h3_index, &mut fijk_roundtrip);
            assert!(
              result.is_ok(),
              "Roundtrip _h3ToFaceIjk failed for {:?}->{:x}",
              fijk_orig,
              h3_index.0
            );

            assert_eq!(
              fijk_roundtrip.face, fijk_orig.face,
              "Face mismatch for res 0. Orig: {:?}, H3: {:x}, Roundtrip: {:?}",
              fijk_orig, h3_index.0, fijk_roundtrip
            );
            assert!(
              _ijk_matches(&fijk_roundtrip.coord, &fijk_orig.coord),
              "IJK mismatch for res 0. Orig: {:?}, H3: {:x}, Roundtrip: {:?}",
              fijk_orig,
              h3_index.0,
              fijk_roundtrip
            );
          }
        }
      }
    }
  }

  #[test]
  fn test_face_ijk_h3_roundtrip_res0_all_valid_inputs() {
    for face_idx in 0..NUM_ICOSA_FACES {
      for i_val in 0..=2 {
        // MAX_FACE_COORD for res 0
        for j_val in 0..=2 {
          for k_val in 0..=2 {
            let mut fijk_orig = FaceIJK {
              face: face_idx,
              coord: CoordIJK {
                i: i_val,
                j: j_val,
                k: k_val,
              },
            };

            // Check if this Fijk would map to a valid base cell and rotation
            // This pre-validates our input to _face_ijk_to_h3
            let bc_check = _face_ijk_to_base_cell(&fijk_orig);
            let rot_check = _face_ijk_to_base_cell_ccwrot60(&fijk_orig);

            if bc_check == INVALID_BASE_CELL || rot_check == INVALID_ROTATIONS {
              continue; // This Fijk doesn't represent a base cell center on this face
            }

            // Normalize the original coord *after* checking its validity for lookup,
            // because the lookup uses raw i,j,k for table indexing.
            // The H3 index produced will correspond to the *normalized* version.
            let mut normalized_ijk_orig = fijk_orig.coord;
            _ijk_normalize(&mut normalized_ijk_orig);
            // The face remains the same for the canonical representation used by _h3_to_face_ijk
            let canonical_fijk_orig_for_comparison = FaceIJK {
              face: fijk_orig.face,
              coord: normalized_ijk_orig,
            };

            let h3_index = _face_ijk_to_h3(&fijk_orig, 0);
            assert_ne!(
              h3_index, H3_NULL,
              "Expected valid H3Index from valid FaceIJK input: {:?}",
              fijk_orig
            );

            let mut fijk_roundtrip = FaceIJK::default();
            let result = _h3_to_face_ijk(h3_index, &mut fijk_roundtrip);
            assert!(
              result.is_ok(),
              "Roundtrip _h3ToFaceIjk failed for {:?} -> {:x}",
              fijk_orig,
              h3_index.0
            );

            // The round-tripped fijk should match the *canonical home fijk* of the base cell.
            // The original fijk_orig might have been an alias on a different face.
            let base_cell_of_h3 = get_base_cell(h3_index);
            let expected_fijk_rt = BASE_CELL_DATA[base_cell_of_h3 as usize].home_fijk;

            assert_eq!(fijk_roundtrip.face, expected_fijk_rt.face,
                        "Face mismatch for res 0. OrigInputFijk: {:?}, H3: {:x} (BC: {}), ExpectedRT_Fijk: {:?}, ActualRT_Fijk: {:?}",
                        fijk_orig, h3_index.0, base_cell_of_h3, expected_fijk_rt, fijk_roundtrip);
            assert!(_ijk_matches(&fijk_roundtrip.coord, &expected_fijk_rt.coord),
                        "IJK mismatch for res 0. OrigInputFijk: {:?}, H3: {:x} (BC: {}), ExpectedRT_Fijk: {:?}, ActualRT_Fijk: {:?}",
                        fijk_orig, h3_index.0, base_cell_of_h3, expected_fijk_rt, fijk_roundtrip);
          }
        }
      }
    }
  }

  fn roundtrip_h3_through_fijk_at_res(h_coarse: H3Index, target_res: i32) {
    let mut iter = crate::iterators::iterInitParent(h_coarse, target_res);
    while iter.h != H3_NULL {
      let fijk_producing_h3_child = iter.h;

      let mut fijk_canonical = FaceIJK::default();
      let h3_to_fijk_res = _h3_to_face_ijk(fijk_producing_h3_child, &mut fijk_canonical);
      assert!(
        h3_to_fijk_res.is_ok(),
        "H3 to FaceIJK failed for {:x}",
        fijk_producing_h3_child.0
      );

      let h3_roundtrip = _face_ijk_to_h3(&fijk_canonical, target_res);
      assert_eq!(
        h3_roundtrip, fijk_producing_h3_child,
        "H3->FaceIJK->H3 roundtrip failed. Orig H3: {:x}, Canonical Fijk: {:?}, Roundtrip H3: {:x}",
        fijk_producing_h3_child.0, fijk_canonical, h3_roundtrip.0
      );

      crate::iterators::iterStepChild(&mut iter);
    }
  }

  #[test]
  fn test_h3_fijk_h3_roundtrip_finer_res() {
    // Test for a few base cells to keep test time reasonable.
    // Base cell 0 (hexagon)
    let bc0 = _face_ijk_to_h3(
      &FaceIJK {
        face: 1,
        coord: CoordIJK { i: 1, j: 0, k: 0 },
      },
      0,
    );
    assert_ne!(bc0, H3_NULL);
    roundtrip_h3_through_fijk_at_res(bc0, 1);
    roundtrip_h3_through_fijk_at_res(bc0, 2);

    // Base cell 4 (pentagon)
    let bc4 = _face_ijk_to_h3(
      &FaceIJK {
        face: 0,
        coord: CoordIJK { i: 2, j: 0, k: 0 },
      },
      0,
    );
    assert_ne!(bc4, H3_NULL);
    roundtrip_h3_through_fijk_at_res(bc4, 1);
    roundtrip_h3_through_fijk_at_res(bc4, 2);

    // A non-central base cell
    let bc15 = _face_ijk_to_h3(
      &FaceIJK {
        face: 4,
        coord: CoordIJK { i: 1, j: 0, k: 0 },
      },
      0,
    );
    assert_ne!(bc15, H3_NULL);
    roundtrip_h3_through_fijk_at_res(bc15, 1);
  }
  #[test]
  fn test_overage_cases() {
    // Example: H3Index 820c7ffffffffff (Res 2, BC 2, D1=1(K), D2=0(Center))
    // Base Cell 2 home: Face 1, IJK {0,0,0}
    // This cell is on face 6.
    let h_overage = H3Index(0x820c7ffffffffff); // Use from_str_radix if H3_NULL is not 0
    let mut fijk_overage = FaceIJK::default();
    assert!(_h3_to_face_ijk(h_overage, &mut fijk_overage).is_ok());

    assert_eq!(fijk_overage.face, 6, "Overage cell expected on face 6");
    // The exact IJK coord here would need to be derived from C H3 output for this specific index.
    // For example, let's assume it's {1,0,0} on face 6 for this cell.
    let expected_overage_ijk = CoordIJK { i: 1, j: 0, k: 0 }; // Placeholder, get actual value
                                                              // assert!(_ijk_matches(&fijk_overage.coord, &expected_overage_ijk), "Overage IJK coord mismatch");

    let h_roundtrip = _face_ijk_to_h3(&fijk_overage, get_resolution(h_overage));
    assert_eq!(h_roundtrip, h_overage, "Overage case H3->Fijk->H3 roundtrip failed");

    // Another case: A cell whose FaceIJK from _geoToFaceIjk would be near an edge,
    // then _faceIjkToH3 needs to apply rotations correctly.
    // E.g., Fijk {face:0, coord:{2,2,1}} at res 0 (sum=5, clearly overage as max_dim=2)
    // _faceIjkToH3 should correctly map this to a base cell (likely on an adjacent face)
    // and then _h3ToFaceIjk for that base cell should yield its canonical home FaceIJK.
    let fijk_input_overage = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 2, j: 2, k: 1 },
    };
    let h3_from_overage_input = _face_ijk_to_h3(&fijk_input_overage, 0);
    assert_ne!(
      h3_from_overage_input, H3_NULL,
      "Expected valid H3 from overage Fijk input"
    );

    let mut fijk_rt_from_overage = FaceIJK::default();
    assert!(_h3_to_face_ijk(h3_from_overage_input, &mut fijk_rt_from_overage).is_ok());

    // The fijk_rt_from_overage should be the *canonical home Fijk* of the resulting base cell.
    let bc_of_h3 = get_base_cell(h3_from_overage_input);
    let canonical_home_fijk = BASE_CELL_DATA[bc_of_h3 as usize].home_fijk;

    assert_eq!(
      fijk_rt_from_overage.face, canonical_home_fijk.face,
      "Face for H3 from overage Fijk"
    );
    assert!(
      _ijk_matches(&fijk_rt_from_overage.coord, &canonical_home_fijk.coord),
      "IJK for H3 from overage Fijk"
    );
  }

  #[test]
  fn test_pentagon_k_axis_rotation_in_face_ijk_to_h3() {
    // We need an fijk that resolves to a pentagon base cell (e.g., BC4)
    // and whose H3 digits determined *before* the final rotation step
    // would have a leading K.

    // BC4's home fijk is {face:0, coord:{2,0,0}}
    // Let's construct an fijk at res=1 that should map to BC4, and try to
    // make its H3 digit effectively KAxes before the final rotation.
    // This means its "diff" vector in _face_ijk_to_h3's loop should be KAxes's unit vector.
    // last_center for res=1, parent_ijk={2,0,0} (BC4 home):
    //   _down_ap7({2,0,0}) -> {4,0,0} (normalized)
    // We want diff = K_AXES_UNIT_VECTOR = {0,0,1} (normalized)
    // So, last_ijk_for_digit = last_center + diff = {4,0,0} + {0,0,1} = {4,0,1}
    // This {4,0,1} is the fijk.coord at res=1 on face 0 that should produce KAxes digit for BC4.
    let fijk_input_k_oriented = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 4, j: 0, k: 1 },
    };
    let res = 1;

    let mut h_result = _face_ijk_to_h3(&fijk_input_k_oriented, res);
    // Expected: H3Index for BC4, Res 1.
    // Digit for res 1 was KAxes (1).
    // _h3_leading_non_zero_digit(h_with_K_digit) would be KAxes.
    // BC4 on face 0 is NOT CW offset. So _h3_rotate_pent60_ccw is applied.
    // KAxes (1) rotated CCW is IkAxes (5).
    // So the final H3Index should have D1 = IkAxes.

    assert_ne!(h_result, H3_NULL, "Expected valid H3Index");
    assert_eq!(get_mode(h_result), H3_CELL_MODE as u8, "Mode should be cell");
    assert_eq!(get_resolution(h_result), res, "Resolution should be 1");
    assert_eq!(get_base_cell(h_result), 4, "Base cell should be 4 (pentagon)");
    assert_eq!(
      get_index_digit(h_result, 1),
      Direction::IkAxes,
      "Digit 1 should be IkAxes after pentagonal K-axis rotation"
    );

    // Now, test the reverse with the canonical H3Index produced.
    let mut fijk_roundtrip = FaceIJK::default();
    assert!(_h3_to_face_ijk(h_result, &mut fijk_roundtrip).is_ok());

    // The canonical fijk for this H3Index (BC4, Res1, D1=IkAxes) needs to be known.
    // BC4 home: {0,{2,0,0}}.
    // _h3_to_face_ijk_with_initialized_fijk(h_result, &mut fijk_home_of_bc4):
    //   h_result (D1=IkAxes=5). BC4 is pentagon. LeadingNonZero is IkAxes (5). NOT KAxes. So NO _h3Rotate60cw on h_result.
    //   fijk.coord starts as {2,0,0}.
    //   r=1 (ClassIII): _down_ap7 on {2,0,0} -> {4,0,0}.
    //   _neighbor on {4,0,0} by digit IkAxes (5): {4,0,0} + {1,0,1} = {5,0,1}. Normalized: {4,0,0}.
    //   No, {5,0,1} normalized is {4,0,0}. What?
    //   UNIT_VECS[5] = {1,0,1}. {4,0,0} + {1,0,1} = {5,0,1}.
    //   _ijk_normalize of {5,0,1}: i=5,j=0,k=1. All non-negative. Min is 0. Stays {5,0,1}.
    // So, canonical fijk should be {face:0, coord:{5,0,1}}.
    let fijk_expected_canonical = FaceIJK {
      face: 0,
      coord: CoordIJK { i: 5, j: 0, k: 1 },
    };
    assert_eq!(
      fijk_roundtrip.face, fijk_expected_canonical.face,
      "Roundtrip face mismatch"
    );
    assert!(
      _ijk_matches(&fijk_roundtrip.coord, &fijk_expected_canonical.coord),
      "Roundtrip IJK mismatch. Expected {:?}, Got {:?}",
      fijk_expected_canonical.coord,
      fijk_roundtrip.coord
    );
  }
}
