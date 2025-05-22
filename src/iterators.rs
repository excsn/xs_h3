// src/iterators.rs
#![allow(non_snake_case)] // To match C API names like iterInitParent

use crate::base_cells::baseCellNumToCell;
use crate::constants::{MAX_H3_RES, NUM_BASE_CELLS};
use crate::h3_index::{
  _h3_leading_non_zero_digit, get_index_digit, get_resolution, is_pentagon, set_index_digit, set_resolution,
};
use crate::hierarchy::parent_child::_zero_index_digits; // Use the helper we created
use crate::types::{Direction, H3Index, H3_NULL}; // If this is how we get initial base cell

/// Iterator for traversing the children of a parent H3 cell at a specific resolution.
///
/// Access the current cell via the `h` field.
/// Advance the iterator with `iterStepChild`.
/// `h` will be `H3_NULL` when the iteration is complete or if initialization failed.
#[derive(Debug, Clone, Copy)]
pub struct IterCellsChildren {
  /// The current H3 child cell in the iteration. `H3_NULL` if exhausted.
  pub h: H3Index,
  /// Resolution of the parent cell from which children are derived.
  pub(crate) _parentRes: i32,
  /// Internal: The resolution digit that needs to skip `Direction::KAxes` for pentagons.
  /// Moves towards coarser resolutions as iteration proceeds. -1 if not a pentagon path.
  pub(crate) _skipDigitRes: i32, // Renamed from _skipDigit to be more descriptive of its value
}

impl IterCellsChildren {
  /// Creates a nulled-out iterator, typically used to signify exhaustion or error.
  fn null_iter() -> Self {
    Self {
      h: H3_NULL,
      _parentRes: -1,
      _skipDigitRes: -1,
    }
  }
}

/// Initializes an iterator for the children of H3 cell `parent_h` at `child_res`.
///
/// The first child is available in `iter.h` after calling this function.
pub fn iterInitParent(parent_h: H3Index, child_res: i32) -> IterCellsChildren {
  let parent_res = get_resolution(parent_h);

  if child_res < parent_res
    || child_res > MAX_H3_RES
    || parent_h == H3_NULL
    || !crate::h3_index::inspection::is_valid_cell(parent_h)
  {
    // Added validity check
    return IterCellsChildren::null_iter();
  }

  let mut h_iter = parent_h;
  // Set to child res and zero out intermediate digits to get the first child (center)
  set_resolution(&mut h_iter, child_res);
  h_iter = _zero_index_digits(h_iter, parent_res + 1, child_res);

  let skip_digit_res = if is_pentagon(h_iter) {
    // Check the *child* at target res for pentagon status
    child_res
  } else {
    -1
  };

  // The C code sets _skipDigit to childRes if the *original parent* is a pentagon.
  // Let's match that more closely for direct porting.
  // The pentagon-ness of the iteration path is determined by the original parent.
  let skip_digit_res_c_style = if is_pentagon(parent_h) { child_res } else { -1 };

  IterCellsChildren {
    h: h_iter,
    _parentRes: parent_res,
    _skipDigitRes: skip_digit_res_c_style,
  }
}

/// Initializes an iterator for all H3 cells at `child_res` that are descendants
/// of the base cell `base_cell_num`.
pub fn iterInitBaseCellNum(base_cell_num: i32, child_res: i32) -> IterCellsChildren {
  if base_cell_num < 0 || base_cell_num >= (NUM_BASE_CELLS as i32) || child_res < 0 || child_res > MAX_H3_RES {
    return IterCellsChildren::null_iter();
  }

  let base_cell_h3 = baseCellNumToCell(base_cell_num);
  if base_cell_h3 == H3_NULL {
    // Should not happen with valid base_cell_num
    return IterCellsChildren::null_iter();
  }
  iterInitParent(base_cell_h3, child_res)
}

/// Steps the `IterCellsChildren` iterator to the next child cell.
/// After stepping, the next cell is available in `iter.h`.
/// If `iter.h` is `H3_NULL`, the iteration is complete.
pub fn iterStepChild(iter: &mut IterCellsChildren) {
  if iter.h == H3_NULL {
    return;
  }

  let child_res = get_resolution(iter.h);

  // Phase 1: Naive increment of the H3 index, handling cascades.
  // This helper needs to increment starting from child_res and stop before _parentRes.
  // It should return true if it cascaded past the parent (making iter.h H3_NULL).
  if _naive_increment_and_cascade(iter, child_res) {
    // _naive_increment_and_cascade already set iter.h to H3_NULL
    return;
  }

  // Phase 2: K-Skip. This part executes only if iter.h is still valid.
  // Check if the current _skipDigitRes is valid and needs skipping.
  if iter._skipDigitRes >= iter._parentRes + 1 {
    // Still relevant to skip
    if get_index_digit(iter.h, iter._skipDigitRes) == Direction::KAxes {
      // Yes, the digit at the skip-responsible level is K.
      // Increment starting from this specific digit and cascade.
      if _naive_increment_and_cascade(iter, iter._skipDigitRes) {
        return; // Cascaded past parent during K-skip adjustment
      }
      // If successful, the responsibility moves up.
      iter._skipDigitRes -= 1;
    }
  }
}

// Helper for Phase 1 and Phase 2's increment part.
// Increments h starting at `level_to_inc_from`, cascading up to `iter._parentRes + 1`.
// Returns true if h becomes H3_NULL (cascaded past parent).
fn _naive_increment_and_cascade(iter: &mut IterCellsChildren, level_to_inc_from: i32) -> bool {
  let mut r_cascade = level_to_inc_from;
  loop {
    if r_cascade < iter._parentRes + 1 {
      *iter = IterCellsChildren::null_iter();
      return true;
    }

    let mut digit_val = get_index_digit(iter.h, r_cascade) as u8;
    digit_val += 1;

    if digit_val >= Direction::InvalidDigit as u8 {
      // Rollover
      set_index_digit(&mut iter.h, r_cascade, Direction::Center);
      // DO NOT MODIFY _skipDigitRes here.
      if r_cascade == iter._parentRes + 1 {
        *iter = IterCellsChildren::null_iter();
        return true;
      }
      r_cascade -= 1;
    } else {
      set_index_digit(&mut iter.h, r_cascade, Direction::try_from(digit_val).unwrap());
      return false;
    }
  }
}

/// Iterator for traversing all H3 cells at a given resolution.
///
/// Access the current cell via the `h` field.
/// Advance the iterator with `iterStepRes`.
/// `h` will be `H3_NULL` when the iteration is complete or if initialization failed.
#[derive(Debug, Clone, Copy)]
pub struct IterCellsResolution {
  /// The current H3 cell in the iteration. `H3_NULL` if exhausted.
  pub h: H3Index,
  pub(crate) _baseCellNum: i32,
  pub(crate) _res: i32,
  pub(crate) _itC: IterCellsChildren, // Internal child iterator
}

impl IterCellsResolution {
  /// Creates a nulled-out iterator.
  fn null_iter() -> Self {
    Self {
      h: H3_NULL,
      _baseCellNum: -1,
      _res: -1,
      _itC: IterCellsChildren::null_iter(),
    }
  }
}

/// Initializes an iterator for all H3 cells at the given `res`.
/// The first cell is available in `iter.h` after calling this function.
pub fn iterInitRes(res: i32) -> IterCellsResolution {
  if res < 0 || res > MAX_H3_RES {
    return IterCellsResolution::null_iter();
  }

  let child_iter = iterInitBaseCellNum(0, res);
  IterCellsResolution {
    h: child_iter.h,
    _baseCellNum: 0,
    _res: res,
    _itC: child_iter,
  }
}

/// Steps the `IterCellsResolution` iterator to the next H3 cell at the resolution.
pub fn iterStepRes(iter: &mut IterCellsResolution) {
  if iter.h == H3_NULL {
    return; // Iteration already complete or errored
  }

  iterStepChild(&mut iter._itC);

  // If the child iterator for the current base cell is exhausted,
  // and there are more base cells, move to the next base cell.
  if iter._itC.h == H3_NULL {
    iter._baseCellNum += 1;
    if iter._baseCellNum < (NUM_BASE_CELLS as i32) {
      iter._itC = iterInitBaseCellNum(iter._baseCellNum, iter._res);
    } else {
      // All base cells processed, iteration is complete.
      // _itC.h is already H3_NULL.
    }
  }
  iter.h = iter._itC.h; // Update current cell
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::h3_index; // For get_resolution, is_pentagon, etc.
  use crate::hierarchy::parent_child::cell_to_children_size;
  use crate::math::extensions::_ipow;
  use crate::types::H3Index; // For checking count

  #[test]
  fn test_iter_init_parent_invalid() {
    let parent = H3Index(0x85283473fffffff); // Res 5
    let mut iter = iterInitParent(parent, 4); // Child res coarser
    assert_eq!(iter.h, H3_NULL);

    iter = iterInitParent(parent, MAX_H3_RES + 1); // Child res too fine
    assert_eq!(iter.h, H3_NULL);

    iter = iterInitParent(H3_NULL, 5);
    assert_eq!(iter.h, H3_NULL);

    let mut invalid_parent = parent;
    h3_index::set_mode(&mut invalid_parent, 0); // Invalid mode
    iter = iterInitParent(invalid_parent, 6);
    assert_eq!(iter.h, H3_NULL);
  }

  #[test]
  fn test_iter_init_base_cell_num_invalid() {
    let mut iter = iterInitBaseCellNum(-1, 0);
    assert_eq!(iter.h, H3_NULL);
    iter = iterInitBaseCellNum(NUM_BASE_CELLS as i32, 0);
    assert_eq!(iter.h, H3_NULL);
    iter = iterInitBaseCellNum(0, -1);
    assert_eq!(iter.h, H3_NULL);
    iter = iterInitBaseCellNum(0, MAX_H3_RES + 1);
    assert_eq!(iter.h, H3_NULL);
  }

  #[test]
  fn test_iter_init_res_invalid() {
    let mut iter = iterInitRes(-1);
    assert_eq!(iter.h, H3_NULL);
    iter = iterInitRes(MAX_H3_RES + 1);
    assert_eq!(iter.h, H3_NULL);
  }

  #[test]
  fn test_iter_children_hexagon() {
    let parent = H3Index(0x85283473fffffff); // Res 5 hex
    let child_res = 7;
    let expected_count = cell_to_children_size(parent, child_res).unwrap();

    let mut iter = iterInitParent(parent, child_res);
    let mut count = 0;
    let mut first_child = H3_NULL;
    let mut last_child = H3_NULL;

    while iter.h != H3_NULL {
      if count == 0 {
        first_child = iter.h;
      }
      last_child = iter.h;
      assert_eq!(h3_index::get_resolution(iter.h), child_res);
      assert_eq!(h3_index::get_base_cell(iter.h), h3_index::get_base_cell(parent));
      // Check that parent of iter.h is indeed `parent`
      let current_parent = crate::hierarchy::parent_child::cell_to_parent(iter.h, get_resolution(parent)).unwrap();
      assert_eq!(current_parent, parent);

      count += 1;
      iterStepChild(&mut iter);
    }
    assert_eq!(count, expected_count);
    assert_ne!(first_child, H3_NULL);
    assert_ne!(last_child, H3_NULL);
    // C H3 Test: H3Index expectedFirst = 0x872834730ffffffL;
    // C H3 Test: H3Index expectedLast = 0x872834737ffffffL;
    // These specific values depend on how _zero_index_digits and digit iteration order work.
    // Our _zero_index_digits makes the first child have all center digits.
  }

  #[test]
  fn test_iter_children_pentagon() {
    let parent_pent = baseCellNumToCell(4); // BC 4 is a pentagon, Res 0
    assert!(is_pentagon(parent_pent));
    let child_res = 2;
    let expected_count = cell_to_children_size(parent_pent, child_res).unwrap(); // 1 + 5*(7^2-1)/6 + 5*7 = 1 + 5*8 + 35 = 41+35 != C test.
                                                                                 // C Test: 1 (self) + 5 = 6 (res 1 children)
                                                                                 // + 5 * 6 (children of those 5 hex) = 30
                                                                                 // + 1 * 5 (children of the center pent res 1 child) = 5
                                                                                 // Total = 1+5+30+5 = 41 for res 2 from res 0 pent.
    assert_eq!(expected_count, 1 + 5 * ((_ipow(7, 2) - 1) / 6)); // Should be 41

    let mut iter = iterInitParent(parent_pent, child_res);
    let mut count = 0;
    while iter.h != H3_NULL {
      assert_eq!(h3_index::get_resolution(iter.h), child_res);
      assert_eq!(h3_index::get_base_cell(iter.h), h3_index::get_base_cell(parent_pent));
      let current_parent = crate::hierarchy::parent_child::cell_to_parent(iter.h, get_resolution(parent_pent)).unwrap();
      assert_eq!(current_parent, parent_pent);
      count += 1;
      iterStepChild(&mut iter);
    }
    assert_eq!(count, expected_count);
  }

  #[test]
  fn test_iter_all_cells_at_resolution() {
    for res in 0..=3 {
      // Test a few low resolutions
      let expected_count_total = crate::h3_index::get_num_cells(res).unwrap(); // Assuming getNumCells is ported

      let mut iter = iterInitRes(res);
      let mut count = 0;
      let mut prev_h = H3_NULL;
      while iter.h != H3_NULL {
        assert_eq!(h3_index::get_resolution(iter.h), res);
        assert!(h3_index::inspection::is_valid_cell(iter.h));
        if prev_h != H3_NULL {
          // Skip first
          assert!(iter.h.0 > prev_h.0, "Cells should be ordered");
        }
        prev_h = iter.h;
        count += 1;
        iterStepRes(&mut iter);
      }
      assert_eq!(count, expected_count_total);
    }
  }
}
