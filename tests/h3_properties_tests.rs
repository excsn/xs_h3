// tests/h3_properties_tests.rs

use xs_h3::{h3_index::{get_resolution, inspection::is_res_class_iii}, *};

#[test]
fn test_cli_get_resolution() {
  // "getResolution -c 85283473fffffff" "5"
  let cell = H3Index(0x85283473fffffff);
  assert_eq!(get_resolution(cell), 5);
}

#[test]
fn test_cli_get_base_cell_number() {
  // From CLI test: "getBaseCellNumber -c 85283473fffffff" "20"
  let cell = H3Index(0x85283473fffffff);

  // Call the *actual library function* being tested
  let base_cell_result = get_base_cell_number(cell);

  assert_eq!(
    base_cell_result, 20,
    "Expected base cell 20 for H3Index {:x}, got {}",
    cell.0, base_cell_result
  );
}

#[test]
fn test_cli_is_pentagon() {
  // "isPentagon -c 85283473fffffff" "false"
  let cell_hex = H3Index(0x85283473fffffff);
  assert_eq!(is_pentagon(cell_hex), false);

  let cell_pent = H3Index(0x8009fffffffffff); // BC 4 is a pentagon
  assert_eq!(is_pentagon(cell_pent), true);
}

#[test]
fn test_cli_is_res_class_iii() {
  // "isResClassIII -c 85283473fffffff" "true" (Res 5 is Class III)
  let cell_res5 = H3Index(0x85283473fffffff);
  assert_eq!(is_res_class_iii(cell_res5), true);

  let cell_res4 = H3Index(0x8428347ffffffff); // Res 4 is Class II
  assert_eq!(is_res_class_iii(cell_res4), false);
}

#[test]
fn test_cli_is_valid_cell() {
  // "isValidCell -c 85283473fffffff" "true"
  let cell_valid = H3Index(0x85283473fffffff);
  assert_eq!(is_valid_cell(cell_valid), true);

  // "isValidCell -c 85283473ffff" "false" (too short, implies invalid structure)
  // Our stringToH3 would fail for "85283473ffff".
  // Let's test a structurally invalid H3Index if constructed directly.
  let cell_invalid_struct = H3Index(0x05283473fffffff); // Mode 0 (invalid for cell)
  assert_eq!(is_valid_cell(cell_invalid_struct), false);
}

#[test]
fn test_cli_string_to_h3_and_h3_to_string() {
  // "stringToInt -c 85283473fffffff" "599686042433355775"
  // "intToString -c 599686042433355775" "85283473fffffff"
  let h3_str = "85283473fffffff";
  let h3_val_u64 = 599686042433355775u64;
  let h3_index = H3Index(h3_val_u64);

  match string_to_h3(h3_str) {
    Ok(parsed_h3) => assert_eq!(parsed_h3, h3_index, "string_to_h3 mismatch"),
    Err(e) => panic!("string_to_h3 failed: {:?}", e),
  }

  let mut buffer = [0u8; 17];
  match h3_to_string(h3_index, &mut buffer) {
    Ok(()) => {
      let result_str = std::str::from_utf8(&buffer[0..h3_str.len()]).unwrap();
      assert_eq!(result_str, h3_str, "h3_to_string mismatch");
    }
    Err(e) => panic!("h3_to_string failed: {:?}", e),
  }

  assert_eq!(h3_to_string_alloc(h3_index), h3_str, "h3_to_string_alloc mismatch");
}

// getNumCells, getPentagons, getRes0Cells, pentagonCount, getIcosahedronFaces
// These will be tested when their respective functions are fully ported.
