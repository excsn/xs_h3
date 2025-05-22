// tests/traversal_tests.rs

use std::collections::HashSet;
use xs_h3::*;

#[test]
fn test_cli_are_neighbor_cells() {
  // "areNeighborCells -o 85283473fffffff -d 85283477fffffff" "true"
  let origin1 = H3Index(0x85283473fffffff);
  let dest1 = H3Index(0x85283477fffffff);
  assert_eq!(are_neighbor_cells(origin1, dest1), Ok(true));

  // "areNeighborCells -o 85283473fffffff -d 85283472fffffff" "false"
  // This H3 index (0x85283472fffffff) is actually invalid because its digits
  // D3, D4, D5 are 7 (InvalidDigit) for res 5.
  // Rust's is_valid_cell correctly catches this.
  let dest2 = H3Index(0x85283472fffffff);
  assert_eq!(are_neighbor_cells(origin1, dest2), Err(H3Error::CellInvalid)); // Corrected expectation

  // "areNeighborCells -o 85283473fffffff -d 852834727fffffff 2>&1" "Error 5: Cell argument was not valid"
  let dest_invalid_by_cli = H3Index(0x852834727fffffff); // Also invalid due to digits (D4 is 7)
  assert_eq!(
    are_neighbor_cells(origin1, dest_invalid_by_cli),
    Err(H3Error::CellInvalid)
  );
}

#[test]
fn test_cli_grid_distance() {
  // "gridDistance -o 85283473fffffff -d 8528342bfffffff" "2"
  let origin = H3Index(0x85283473fffffff);
  let destination = H3Index(0x8528342bfffffff);
  assert_eq!(grid_distance(origin, destination), Ok(2));
}

fn parse_h3_list_from_cli_output(s: &str) -> Vec<H3Index> {
  s.trim_matches(|c| c == '[' || c == ']' || c == ' ')
    .split(',')
    .filter_map(|item| {
      let cleaned = item.trim().trim_matches('"');
      if cleaned.is_empty() {
        None
      } else {
        string_to_h3(cleaned).ok()
      }
    })
    .collect()
}

#[test]
fn test_cli_grid_disk() {
  // "gridDisk -c 85283473fffffff -k 1"
  // "[ \"85283473fffffff\", \"85283447fffffff\", \"8528347bfffffff\", \"85283463fffffff\", \"85283477fffffff\", \"8528340ffffffff\", \"8528340bfffffff\" ]"
  let origin = H3Index(0x85283473fffffff);
  let k = 1;
  let expected_str = "[ \"85283473fffffff\", \"85283447fffffff\", \"8528347bfffffff\", \"85283463fffffff\", \"85283477fffffff\", \"8528340ffffffff\", \"8528340bfffffff\" ]";
  let expected_cells = parse_h3_list_from_cli_output(expected_str);

  let max_size = max_grid_disk_size(k).unwrap() as usize;
  let mut out_cells = vec![H3_NULL; max_size];
  assert!(grid_disk(origin, k, &mut out_cells).is_ok());

  let result_set: HashSet<H3Index> = out_cells.into_iter().filter(|&h| h != H3_NULL).collect();
  let expected_set: HashSet<H3Index> = expected_cells.into_iter().collect();

  assert_eq!(
    result_set.len(),
    expected_set.len(),
    "Number of cells mismatch for gridDisk"
  );
  assert_eq!(result_set, expected_set, "Cell content mismatch for gridDisk");
}

#[test]
fn test_cli_grid_ring() {
  // "gridRing -c 85283473fffffff -k 1"
  // "[ \"8528340bfffffff\", \"85283447fffffff\", \"8528347bfffffff\", \"85283463fffffff\", \"85283477fffffff\", \"8528340ffffffff\" ]"
  let origin = H3Index(0x85283473fffffff);
  let k = 1;
  let expected_str = "[ \"8528340bfffffff\", \"85283447fffffff\", \"8528347bfffffff\", \"85283463fffffff\", \"85283477fffffff\", \"8528340ffffffff\" ]";
  let expected_cells = parse_h3_list_from_cli_output(expected_str);

  let ring_size = if k == 0 { 1 } else { (6 * k) as usize };
  let mut out_cells = vec![H3_NULL; ring_size];
  assert!(
    grid_ring_unsafe(origin, k, &mut out_cells).is_ok(),
    "grid_ring_unsafe failed"
  );

  let result_set: HashSet<H3Index> = out_cells.into_iter().filter(|&h| h != H3_NULL).collect();
  let expected_set: HashSet<H3Index> = expected_cells.into_iter().collect();

  assert_eq!(
    result_set.len(),
    expected_set.len(),
    "Number of cells mismatch for gridRingUnsafe"
  );
  assert_eq!(result_set, expected_set, "Cell content mismatch for gridRingUnsafe");
}

// TODO: Tests for gridDiskDistances, gridPathCells using CLI data (parsing will be more complex for nested arrays).
