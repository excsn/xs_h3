use crate::constants::{H3_RESERVED_MASK_NEGATIVE, MAX_H3_RES};
use crate::h3_index::inspection::is_valid_cell as h3_is_valid_cell;
use crate::h3_index::{get_reserved_bits, get_resolution, is_pentagon, set_reserved_bits};
use crate::hierarchy::parent_child::{cell_to_children, cell_to_children_size, cell_to_parent};
use crate::iterators::{iterInitParent}; // If cell_to_children uses this
use crate::types::{H3Error, H3Index, H3_NULL};

use std::collections::HashMap; // For a safer compaction approach

/// Calculates the exact number of H3 cells that will result from uncompacting
/// the given set of H3 cells to the target resolution `res`.
///
/// # Arguments
/// * `compacted_set` - A slice of H3 cell indexes (potentially compacted).
/// * `res` - The target resolution for uncompaction. Must be finer than or
///           equal to the resolution of every cell in `compacted_set`.
///
/// # Returns
/// `Ok(count)` of uncompacted cells, or an `H3Error` if input is invalid.
pub fn uncompact_cells_size(compacted_set: &[H3Index], res: i32) -> Result<i64, H3Error> {
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }

  let mut count: i64 = 0;
  for &cell in compacted_set {
    if cell == H3_NULL {
      continue; // Skip null/empty entries
    }
    if !h3_is_valid_cell(cell) {
      return Err(H3Error::CellInvalid);
    } // Validate each cell

    let cell_res = get_resolution(cell);
    if cell_res > res {
      return Err(H3Error::ResMismatch); // Cannot uncompact to a coarser resolution
    }

    let children_count = cell_to_children_size(cell, res)?;
    count = count.saturating_add(children_count);
  }
  Ok(count)
}

/// Uncompacts the given set of H3 cells to the target resolution `res`.
///
/// # Arguments
/// * `compacted_set` - A slice of H3 cell indexes.
/// * `res` - The target resolution for uncompaction.
/// * `out_set` - The output slice to fill with uncompacted cells. Must be
///               correctly sized using `uncompact_cells_size`.
///
/// # Returns
/// `Ok(())` on success, or an `H3Error` if input is invalid or `out_set` is too small.
pub fn uncompact_cells(compacted_set: &[H3Index], res: i32, out_set: &mut [H3Index]) -> Result<(), H3Error> {
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }

  let mut current_idx: usize = 0;
  let out_len = out_set.len();

  for &cell in compacted_set {
    if cell == H3_NULL {
      continue;
    }
    if !h3_is_valid_cell(cell) {
      return Err(H3Error::CellInvalid);
    }

    let cell_res = get_resolution(cell);
    if cell_res > res {
      return Err(H3Error::ResMismatch);
    }

    if cell_res == res {
      // Cell is already at the target resolution
      if current_idx >= out_len {
        return Err(H3Error::MemoryBounds);
      }
      out_set[current_idx] = cell;
      current_idx += 1;
    } else {
      // Uncompact further
      // cell_to_children uses an iterator internally which is efficient.
      // Or, we could reimplement the recursive/iterative uncompaction logic here.
      // For now, let's rely on cell_to_children.
      let children_count = cell_to_children_size(cell, res)? as usize;
      if current_idx + children_count > out_len {
        return Err(H3Error::MemoryBounds);
      }

      // Call cell_to_children to fill the segment of out_set
      // cell_to_children needs a mutable slice.
      cell_to_children(cell, res, &mut out_set[current_idx..(current_idx + children_count)])?;
      current_idx += children_count;
    }
  }

  // If the actual number of cells written is less than claimed by uncompact_cells_size,
  // fill the rest of out_set with H3_NULL (as per C API expectation if numOut was larger).
  // However, uncompact_cells_size should be exact.
  // This loop ensures any remaining slots (if out_set was overallocated by caller) are H3_NULL.
  for i in current_idx..out_len {
    out_set[i] = H3_NULL;
  }

  Ok(())
}

/// Compacts a set of H3 cells of the same resolution, returning the smallest
/// set of H3 cells that cover the same area.
///
/// # Arguments
/// * `h3_set` - A slice of H3 cell indexes, all at the same resolution.
///              The input set may be modified (sorted) by this function.
///              Consider passing a mutable copy if the original order is important.
/// * `out_compacted_set` - Output slice to fill with compacted cells.
///                         Its length should generally be at least `h3_set.len()`
///                         in the worst case (no compaction possible).
///
/// # Returns
/// `Ok(num_compacted_cells)` which is the number of cells written to `out_compacted_set`,
/// or an `H3Error` on failure.
pub fn compact_cells(
  mut h3_set: &mut [H3Index], // Takes mutable slice for in-place sort & modification
  out_compacted_set: &mut [H3Index],
) -> Result<usize, H3Error> {
  if h3_set.is_empty() {
    return Ok(0);
  }

  // Sort the input set to group cells by parent. H3Index sort is numeric.
  h3_set.sort_unstable();

  // Basic validation: all cells must be valid and have the same resolution.
  // Also check for duplicates early if they are problematic.
  // The C version errors on duplicates.
  let initial_res = get_resolution(h3_set[0]);
  if !h3_is_valid_cell(h3_set[0]) {
    return Err(H3Error::CellInvalid);
  }

  for i in 0..h3_set.len() {
    if h3_set[i] == H3_NULL {
      continue;
    } // Allow H3_NULL to be skipped
    if !h3_is_valid_cell(h3_set[i]) {
      return Err(H3Error::CellInvalid);
    }
    if get_resolution(h3_set[i]) != initial_res {
      return Err(H3Error::ResMismatch);
    }
    if i > 0 && h3_set[i] == h3_set[i - 1] && h3_set[i] != H3_NULL {
      return Err(H3Error::DuplicateInput);
    }
  }

  if initial_res == 0 {
    // Nothing to compact at res 0
    let mut count = 0;
    for (i, &cell) in h3_set.iter().enumerate() {
      if cell != H3_NULL {
        if count >= out_compacted_set.len() {
          return Err(H3Error::MemoryBounds);
        }
        out_compacted_set[count] = cell;
        count += 1;
      }
    }
    // Fill rest with H3_NULL if out_compacted_set is larger
    for i in count..out_compacted_set.len() {
      out_compacted_set[i] = H3_NULL;
    }
    return Ok(count);
  }

  // Iterative compaction:
  // `current_res_set` holds cells at the current resolution being processed.
  // `next_res_set_buffer` holds parents from `current_res_set` that might be fully covered.
  // `uncompactable_this_round` holds cells from `current_res_set` whose parents were not fully covered.

  let mut current_res_set: Vec<H3Index> = h3_set.iter().cloned().filter(|&h| h != H3_NULL).collect();
  let mut compacted_count: usize = 0; // Count of cells written to out_compacted_set

  while !current_res_set.is_empty() {
    let current_res = get_resolution(current_res_set[0]);
    if current_res == 0 {
      // Base cells cannot be compacted further
      for cell in current_res_set {
        if compacted_count >= out_compacted_set.len() {
          return Err(H3Error::MemoryBounds);
        }
        out_compacted_set[compacted_count] = cell;
        compacted_count += 1;
      }
      break; // Done
    }

    let parent_res = current_res - 1;
    let mut parent_counts: HashMap<H3Index, u8> = HashMap::new();
    let mut uncompactable_this_round: Vec<H3Index> = Vec::new();

    // Count children per parent
    for &cell in &current_res_set {
      let parent = cell_to_parent(cell, parent_res)?; // Should not error if logic is right
      *parent_counts.entry(parent).or_insert(0) += 1;
    }

    let mut next_res_set_buffer: Vec<H3Index> = Vec::new();

    for (parent, count) in parent_counts {
      let children_needed = if is_pentagon(parent) { 6 } else { 7 }; // Center (is a child) + 5/6 neighbors
      if count == children_needed {
        // All children present, parent is compactable
        next_res_set_buffer.push(parent);
      } else {
        // Not all children present, add original children of this parent to uncompactable
        for &cell in &current_res_set {
          if cell_to_parent(cell, parent_res)? == parent {
            // Should not error
            uncompactable_this_round.push(cell);
          }
        }
      }
    }

    // Add uncompactable cells from this round to the final output
    for cell in uncompactable_this_round {
      if compacted_count >= out_compacted_set.len() {
        return Err(H3Error::MemoryBounds);
      }
      out_compacted_set[compacted_count] = cell;
      compacted_count += 1;
    }

    current_res_set = next_res_set_buffer; // Prepare for next round of compaction
    if !current_res_set.is_empty() {
      current_res_set.sort_unstable(); // Sort for next round grouping
                                       // Deduplicate parents before next round (important if multiple groups compacted to same grandparent)
      current_res_set.dedup();
    }
  }

  // Fill rest of out_compacted_set with H3_NULL
  for i in compacted_count..out_compacted_set.len() {
    out_compacted_set[i] = H3_NULL;
  }

  Ok(compacted_count)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::hierarchy::parent_child::cell_to_children; // For setting up tests
  use crate::types::H3Index; // Explicitly use our H3Index type

  #[test]
  fn test_uncompact_cells_size() {
    let mut compacted = [H3Index(0x85283473fffffff)]; // Res 5 cell
    assert_eq!(uncompact_cells_size(&compacted, 5), Ok(1));
    assert_eq!(uncompact_cells_size(&compacted, 6), Ok(7));
    assert_eq!(uncompact_cells_size(&compacted, 7), Ok(49));
    assert_eq!(uncompact_cells_size(&compacted, 4), Err(H3Error::ResMismatch));
    assert_eq!(uncompact_cells_size(&[H3_NULL], 5), Ok(0)); // Skips H3_NULL

    let mut pent_compacted = [crate::h3_index::_face_ijk_to_h3(
      &crate::types::FaceIJK {
        face: 0,
        coord: crate::types::CoordIJK { i: 2, j: 0, k: 0 },
      },
      0,
    )]; // BC 4 (pentagon)
    assert_eq!(get_resolution(pent_compacted[0]), 0);
    assert!(is_pentagon(pent_compacted[0]));
    assert_eq!(uncompact_cells_size(&pent_compacted, 0), Ok(1));
    assert_eq!(uncompact_cells_size(&pent_compacted, 1), Ok(1 + 5 * (7 - 1) / 6)); // 6
    assert_eq!(uncompact_cells_size(&pent_compacted, 2), Ok(1 + 5 * (49 - 1) / 6));
    // 41
  }

  #[test]
  fn test_uncompact_cells() {
    let compacted = [H3Index(0x85283473fffffff)]; // Res 5 cell

    // Uncompact to same res
    let mut out_res5 = [H3_NULL; 1];
    assert!(uncompact_cells(&compacted, 5, &mut out_res5).is_ok());
    assert_eq!(out_res5[0], compacted[0]);

    // Uncompact to finer res
    let mut out_res6 = vec![H3_NULL; 7]; // Size for children of a hex
    assert!(uncompact_cells(&compacted, 6, &mut out_res6).is_ok());
    let mut children_check = vec![H3_NULL; 7];
    cell_to_children(compacted[0], 6, &mut children_check).unwrap();
    out_res6.sort_unstable(); // cell_to_children might not sort
    children_check.sort_unstable();
    assert_eq!(out_res6, children_check);

    // Test MemoryBounds error
    let mut too_small_out = [H3_NULL; 6];
    assert_eq!(
      uncompact_cells(&compacted, 6, &mut too_small_out),
      Err(H3Error::MemoryBounds)
    );
  }

  #[test]
  fn test_compact_cells_simple_parent() {
    let parent = H3Index(0x85283473fffffff); // Res 5
    let children_count = cell_to_children_size(parent, 6).unwrap() as usize;
    let mut children = vec![H3_NULL; children_count];
    cell_to_children(parent, 6, &mut children).unwrap();

    let mut compacted_out = vec![H3_NULL; children_count]; // Max possible size
    let compacted_len = compact_cells(&mut children, &mut compacted_out).unwrap();

    assert_eq!(compacted_len, 1);
    assert_eq!(compacted_out[0], parent);
  }

  #[test]
  fn test_compact_cells_no_compaction() {
    let mut cells = [
      H3Index(0x86283470fffffff), // Res 6, child 0 of parent above
      H3Index(0x86283472fffffff), // Res 6, child 1 of parent above
      H3Index(0x86283474fffffff), // Res 6, child 2 of parent above
    ];
    let mut cells_clone = cells; // Clone for assertion, as compact_cells sorts

    let mut compacted_out = vec![H3_NULL; cells.len()];
    let compacted_len = compact_cells(&mut cells, &mut compacted_out).unwrap();

    cells_clone.sort_unstable(); // Sort original to compare content
    compacted_out[0..compacted_len].sort_unstable();

    assert_eq!(compacted_len, 3);
    assert_eq!(&compacted_out[0..compacted_len], &cells_clone[..]);
  }

  #[test]
  fn test_compact_cells_duplicate_input() {
    let mut cells = [
      H3Index(0x86283470fffffff),
      H3Index(0x86283470fffffff), // Duplicate
    ];
    let mut compacted_out = vec![H3_NULL; cells.len()];
    assert_eq!(
      compact_cells(&mut cells, &mut compacted_out),
      Err(H3Error::DuplicateInput)
    );
  }
  #[test]
  fn test_compact_cells_mixed_res() {
    let mut cells = [
      H3Index(0x85283473fffffff), // Res 5
      H3Index(0x86283470fffffff), // Res 6
    ];
    let mut compacted_out = vec![H3_NULL; cells.len()];
    assert_eq!(compact_cells(&mut cells, &mut compacted_out), Err(H3Error::ResMismatch));
  }

  #[test]
  fn test_compact_cells_pentagon_children() {
    let pent_parent = crate::h3_index::_face_ijk_to_h3(
      // BC 4, Res 0
      &crate::types::FaceIJK {
        face: 0,
        coord: crate::types::CoordIJK { i: 2, j: 0, k: 0 },
      },
      0,
    );
    let children_count = cell_to_children_size(pent_parent, 1).unwrap() as usize; // Should be 6
    assert_eq!(children_count, 6);
    let mut children = vec![H3_NULL; children_count];
    cell_to_children(pent_parent, 1, &mut children).unwrap();

    let mut compacted_out = vec![H3_NULL; children_count];
    let compacted_len = compact_cells(&mut children, &mut compacted_out).unwrap();
    assert_eq!(compacted_len, 1);
    assert_eq!(compacted_out[0], pent_parent);
  }
}
