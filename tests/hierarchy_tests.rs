// tests/hierarchy_tests.rs

use xs_h3::*;

#[test]
fn test_cli_cell_to_parent() {
  // "cellToParent -c 8928342e20fffff --resolution 3 -f newline" "832834fffffffff"
  let child = H3Index(0x8928342e20fffff); // Res 9
  let parent_res = 3;
  let expected_parent = H3Index(0x832834fffffffff);

  match cell_to_parent(child, parent_res) {
    Ok(p) => assert_eq!(p, expected_parent),
    Err(e) => panic!("cellToParent failed: {:?}", e),
  }
}

#[test]
fn test_cli_cell_to_children_size() {
  // "cellToChildrenSize -c 85283473fffffff -r 6" "7"
  let parent = H3Index(0x85283473fffffff); // Res 5
  let child_res = 6;
  assert_eq!(cell_to_children_size(parent, child_res), Ok(7));
}

#[test]
fn test_cli_cell_to_center_child() {
  // "cellToCenterChild -c 85283473fffffff --resolution 7 -f newline" "872834700ffffff"
  let parent = H3Index(0x85283473fffffff); // Res 5
  let child_res = 7;
  let expected_center_child = H3Index(0x872834700ffffff);
  match cell_to_center_child(parent, child_res) {
    Ok(cc) => assert_eq!(cc, expected_center_child),
    Err(e) => panic!("cellToCenterChild failed: {:?}", e),
  }
}

#[test]
fn test_cli_cell_to_child_pos() {
  // "cellToChildPos -c 85283473fffffff -r 3" "25" (Child is res 5, parent_res for query is 3)
  let child = H3Index(0x85283473fffffff); // Res 5
  let parent_res = 3;
  assert_eq!(cell_to_child_pos(child, parent_res), Ok(25));
}

#[test]
fn test_cli_child_pos_to_cell() {
  // "childPosToCell -c 85283473fffffff -r 7 -p 42 -f newline" "872834730ffffff"
  // Parent is 85283473fffffff (res 5). Target child_res is 7. Position is 42.
  let parent = H3Index(0x85283473fffffff);
  let child_res = 7;
  let child_pos = 42;
  let expected_child = H3Index(0x872834730ffffff);
  match child_pos_to_cell(child_pos, parent, child_res) {
    Ok(c) => assert_eq!(c, expected_child),
    Err(e) => panic!("childPosToCell failed: {:?}", e),
  }
}

// TODO: Tests for cellToChildren, compactCells, uncompactCells using CLI data.
// These require parsing lists of H3Indexes from CLI output and potentially from fixture files.
