// examples/uncompact_compact.rs
use xs_h3::{
  cell_to_children, cell_to_children_size, compact_cells, h3_to_string_alloc, uncompact_cells, uncompact_cells_size,
  H3Error, H3Index, H3_NULL,
};

fn main() -> Result<(), H3Error> {
  // println!("--- H3 Compaction and Uncompaction Example ---");

  // Start with a cell at a moderate resolution, e.g., res 4
  let parent_cell_str = "8428309ffffffff"; // A res 4 cell
  let parent_cell = xs_h3::string_to_h3(parent_cell_str)?;

  // println!(
  //   "Parent cell: {} (Hex: {:x}) at res {}",
  //   h3_to_string_alloc(parent_cell),
  //   parent_cell.0,
  //   xs_h3::get_resolution(parent_cell)
  // );

  // 1. Get all its children at resolution 6
  let child_res = 6;
  let num_children = cell_to_children_size(parent_cell, child_res)? as usize;
  // println!("Number of children at res {}: {}", child_res, num_children);

  let mut children_cells = vec![H3_NULL; num_children];
  cell_to_children(parent_cell, child_res, &mut children_cells)?;

  // println!("First few children:");
  for i in 0..std::cmp::min(5, num_children) {
    // println!(
    //   "  Child {}: {} ({:x})",
    //   i,
    //   h3_to_string_alloc(children_cells[i]),
    //   children_cells[i].0
    // );
  }

  // 2. Compact this set of children
  // The result of compacting all children of a single parent should be the parent itself.
  let mut compacted_set_output = vec![H3_NULL; num_children]; // Max possible size for output
  let num_compacted = compact_cells(&mut children_cells, &mut compacted_set_output)?;
  // Note: compact_cells sorts its input `children_cells`

  // println!("\nCompacted set (should be 1 cell, the parent):");
  for i in 0..num_compacted {
    // println!(
    //   "  Compacted {}: {} ({:x})",
    //   i,
    //   h3_to_string_alloc(compacted_set_output[i]),
    //   compacted_set_output[i].0
    // );
  }
  assert_eq!(num_compacted, 1);
  assert_eq!(compacted_set_output[0], parent_cell);
  // println!("Compaction successful: children compacted back to the parent.");

  // 3. Uncompact the parent cell back to resolution 6
  let uncompact_target_res = 6;
  let compacted_input_slice = &compacted_set_output[0..num_compacted]; // Slice of valid compacted cells

  let uncompacted_size = uncompact_cells_size(compacted_input_slice, uncompact_target_res)? as usize;
  // println!(
  //   "\nSize of uncompacting the compacted set to res {}: {}",
  //   uncompact_target_res, uncompacted_size
  // );
  assert_eq!(uncompacted_size, num_children); // Should be the original number of children

  let mut uncompacted_cells_output = vec![H3_NULL; uncompacted_size];
  uncompact_cells(
    compacted_input_slice,
    uncompact_target_res,
    &mut uncompacted_cells_output,
  )?;

  // println!("Uncompacted cells (first few):");
  // Sort both original children and uncompacted children to compare content
  children_cells.sort_unstable(); // Was sorted by compact_cells, but let's be sure if we didn't clone
  uncompacted_cells_output.sort_unstable();

  let mut match_count = 0;
  for i in 0..std::cmp::min(5, uncompacted_size) {
    // println!(
    //   "  Uncompacted {}: {} ({:x})",
    //   i,
    //   h3_to_string_alloc(uncompacted_cells_output[i]),
    //   uncompacted_cells_output[i].0
    // );
  }
  // Verify all uncompacted cells match the original children
  let original_children_valid: Vec<H3Index> = children_cells.into_iter().filter(|&h| h != H3_NULL).collect();
  let uncompacted_valid: Vec<H3Index> = uncompacted_cells_output.into_iter().filter(|&h| h != H3_NULL).collect();

  assert_eq!(
    original_children_valid, uncompacted_valid,
    "Uncompacted cells do not match original children set."
  );
  // println!("Uncompaction successful: parent uncompacted back to original children set.");

  Ok(())
}
