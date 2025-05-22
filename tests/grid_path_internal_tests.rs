use xs_h3::coords::ijk::{_ijk_normalize, cube_to_ijk, ijk_to_cube}; // Make cube_to_ijk pub(crate) if not already
use xs_h3::indexing::lat_lng_to_cell;
use xs_h3::local_ij::{cell_to_local_ijk, local_ijk_to_cell};
use xs_h3::traversal::grid_distance; // Assuming grid_distance is pub
use xs_h3::traversal::grid_path::{grid_path_cells, grid_path_cells_size}; // Assuming grid_distance is pub
use xs_h3::types::{CoordIJK, H3Index, LatLng, H3_NULL};

// --- Test Fixtures ---
// Values from the failing test case
const START_H3_FAIL_CASE: H3Index = H3Index(0x855943cbfffffff);
const END_H3_FAIL_CASE: H3Index = H3Index(0x85594303fffffff);

// Expected intermediate values for the failing case
const LOCAL_IJK_START_ORIGIN_EXPECTED: CoordIJK = CoordIJK { i: 0, j: 0, k: 0 };
const LOCAL_IJK_END_ORIGIN_EXPECTED: CoordIJK = CoordIJK { i: -3, j: -3, k: 0 };
const AXIAL_START_EXPECTED: CoordIJK = CoordIJK { i: 0, j: 0, k: 0 };
const AXIAL_END_EXPECTED: CoordIJK = CoordIJK { i: -3, j: -3, k: 6 };
const DISTANCE_EXPECTED: i64 = 3;

// Expected H3 IJK+ coords for each step in grid_path_cells for the failing case (Rust's current output)
const H3_IJK_PLUS_N0_RUST: CoordIJK = CoordIJK { i: 0, j: 0, k: 0 };
const H3_IJK_PLUS_N1_RUST: CoordIJK = CoordIJK { i: 2, j: 0, k: 1 };
const H3_IJK_PLUS_N2_RUST: CoordIJK = CoordIJK { i: 4, j: 0, k: 2 };
const H3_IJK_PLUS_N3_RUST: CoordIJK = CoordIJK { i: -3, j: -3, k: 0 }; // This is end_ijk_local_h3plus

// Expected H3 cells for each step (Rust's current output)
const PATH_N0_RUST: H3Index = H3Index(0x855943cbfffffff);
const PATH_N1_RUST: H3Index = H3Index(0x855943d3fffffff);
const PATH_N2_RUST: H3Index = H3Index(0x85594063fffffff);
const PATH_N3_RUST: H3Index = H3Index(0x85594303fffffff);

// C Path cells for comparison
const C_PATH_N1: H3Index = H3Index(0x855943cffffffff);
const C_PATH_N2: H3Index = H3Index(0x8559431bfffffff);

// Helper function from grid_path.rs (make it pub(crate) or copy for test)
// For testing, it's fine to copy its definition here.
fn lround_c99_style_test(val: f64) -> f64 {
  if val == 0.0 {
    0.0
  } else if val > 0.0 {
    (val + 0.5).floor()
  } else {
    (val - 0.5).ceil()
  }
}

fn ijk_round_to_axial_hex_center_test_version(i_f: f64, j_f: f64, k_f: f64, out_ijk: &mut CoordIJK) {
  let mut ri = lround_c99_style_test(i_f);
  let mut rj = lround_c99_style_test(j_f);
  let mut rk = lround_c99_style_test(k_f);
  let i_diff = (ri - i_f).abs();
  let j_diff = (rj - j_f).abs();
  let k_diff = (rk - k_f).abs();
  if i_diff > j_diff && i_diff > k_diff {
    ri = -rj - rk;
  } else if j_diff > k_diff {
    rj = -ri - rk;
  } else {
    rk = -ri - rj;
  }
  out_ijk.i = ri as i32;
  out_ijk.j = rj as i32;
  out_ijk.k = rk as i32;
}

#[test]
fn test_step1_cell_to_local_ijk() {
  let mut local_ijk_start_origin = CoordIJK::default();
  assert!(cell_to_local_ijk(START_H3_FAIL_CASE, START_H3_FAIL_CASE, &mut local_ijk_start_origin).is_ok());
  assert_eq!(local_ijk_start_origin, LOCAL_IJK_START_ORIGIN_EXPECTED);

  let mut local_ijk_end_origin = CoordIJK::default();
  assert!(cell_to_local_ijk(START_H3_FAIL_CASE, END_H3_FAIL_CASE, &mut local_ijk_end_origin).is_ok());
  assert_eq!(local_ijk_end_origin, LOCAL_IJK_END_ORIGIN_EXPECTED);
}

#[test]
fn test_step2_ijk_to_cube() {
  let mut axial_start = LOCAL_IJK_START_ORIGIN_EXPECTED; // This is H3 IJK+
  ijk_to_cube(&mut axial_start);
  assert_eq!(axial_start, AXIAL_START_EXPECTED);

  let mut axial_end = LOCAL_IJK_END_ORIGIN_EXPECTED; // This is H3 IJK+
  ijk_to_cube(&mut axial_end);
  assert_eq!(axial_end, AXIAL_END_EXPECTED);
}

#[test]
fn test_step3_line_interpolation_steps() {
  let distance = grid_distance(START_H3_FAIL_CASE, END_H3_FAIL_CASE).unwrap();
  assert_eq!(distance, DISTANCE_EXPECTED);
  let inv_distance_f = 1.0 / (distance as f64);
  let i_step = (AXIAL_END_EXPECTED.i - AXIAL_START_EXPECTED.i) as f64 * inv_distance_f;
  let j_step = (AXIAL_END_EXPECTED.j - AXIAL_START_EXPECTED.j) as f64 * inv_distance_f;
  let k_step = (AXIAL_END_EXPECTED.k - AXIAL_START_EXPECTED.k) as f64 * inv_distance_f;
  assert_eq!(i_step, -1.0);
  assert_eq!(j_step, -1.0);
  assert_eq!(k_step, 2.0);
}

#[test]
fn test_step4_ijk_round_to_axial_hex_center() {
  let mut rounded_axial = CoordIJK::default();
  // n=0
  ijk_round_to_axial_hex_center_test_version(0.0, 0.0, 0.0, &mut rounded_axial);
  assert_eq!(rounded_axial, AXIAL_START_EXPECTED); // {0,0,0}
                                                   // n=1 (inputs are -1, -1, 2 from i_f, j_f, k_f calc)
  ijk_round_to_axial_hex_center_test_version(-1.0, -1.0, 2.0, &mut rounded_axial);
  assert_eq!(rounded_axial, CoordIJK { i: -1, j: -1, k: 2 });
  // n=2
  ijk_round_to_axial_hex_center_test_version(-2.0, -2.0, 4.0, &mut rounded_axial);
  assert_eq!(rounded_axial, CoordIJK { i: -2, j: -2, k: 4 });
  // n=3
  ijk_round_to_axial_hex_center_test_version(-3.0, -3.0, 6.0, &mut rounded_axial);
  assert_eq!(rounded_axial, AXIAL_END_EXPECTED); // {-3,-3,6}
}

#[test]
fn test_step4_ijk_round_to_axial_hex_center_tie_breaking() {
  let mut ijk = CoordIJK::default();
  // Test 0.5 cases with lround_c99_style
  // Example: line point is (0.5, -0.5, 0.0)
  // ri=lround(0.5)=1, rj=lround(-0.5)=-1, rk=lround(0)=0. Sum = 0. Correct.
  ijk_round_to_axial_hex_center_test_version(0.5, -0.5, 0.0, &mut ijk);
  assert_eq!(ijk, CoordIJK { i: 1, j: -1, k: 0 });

  // Example: (0.5, 0.5, -1.0) -> sum = 0
  ijk_round_to_axial_hex_center_test_version(0.5, 0.5, -1.0, &mut ijk);
  assert_eq!(ijk, CoordIJK { i: 1, j: 0, k: -1 });

  // Example requiring correction logic
  // (0.6, 0.7, -1.3) -> sum 0. Point is further from integer grid.
  // ri=1, rj=1, rk=-1. i_diff=0.4, j_diff=0.3, k_diff=0.3
  // j_diff == k_diff, and i_diff > j_diff. So first if branch: ri = -rj-rk = -1 -(-1) = 0.
  // Result: {0, 1, -1}
  ijk_round_to_axial_hex_center_test_version(0.6, 0.7, -1.3, &mut ijk); // Input sum is 0
                                                                        // ri=lround(0.6)=1, rj=lround(0.7)=1, rk=lround(-1.3)=-1.
                                                                        // i_diff = |1-0.6|=0.4. j_diff=|1-0.7|=0.3. k_diff=|-1 - (-1.3)| = 0.3.
                                                                        // i_diff is largest. So ri = -rj - rk = -(1) - (-1) = 0.
                                                                        // Output: {0, 1, -1}. Sum is 0.
  assert_eq!(ijk, CoordIJK { i: 0, j: 1, k: -1 });
}

#[test]
fn test_step5_axial_to_h3_ijk_plus_conversion() {
  // From axial {0,0,0} (for n=0)
  let mut h3_ijk_plus = AXIAL_START_EXPECTED; // {0,0,0}
  h3_ijk_plus.i = AXIAL_START_EXPECTED.i.saturating_neg();
  h3_ijk_plus.j = AXIAL_START_EXPECTED.j;
  h3_ijk_plus.k = 0;
  _ijk_normalize(&mut h3_ijk_plus);
  assert_eq!(h3_ijk_plus, H3_IJK_PLUS_N0_RUST);

  // From axial {-1,-1,2} (for n=1)
  let axial_n1 = CoordIJK { i: -1, j: -1, k: 2 };
  let mut h3_ijk_plus_n1 = axial_n1;
  h3_ijk_plus_n1.i = axial_n1.i.saturating_neg();
  h3_ijk_plus_n1.j = axial_n1.j;
  h3_ijk_plus_n1.k = 0;
  _ijk_normalize(&mut h3_ijk_plus_n1);
  assert_eq!(h3_ijk_plus_n1, H3_IJK_PLUS_N1_RUST);

  // From axial {-2,-2,4} (for n=2)
  let axial_n2 = CoordIJK { i: -2, j: -2, k: 4 };
  let mut h3_ijk_plus_n2 = axial_n2;
  h3_ijk_plus_n2.i = axial_n2.i.saturating_neg();
  h3_ijk_plus_n2.j = axial_n2.j;
  h3_ijk_plus_n2.k = 0;
  _ijk_normalize(&mut h3_ijk_plus_n2);
  assert_eq!(h3_ijk_plus_n2, H3_IJK_PLUS_N2_RUST);

  // For n=3 (end point), grid_path_cells uses LOCAL_IJK_END_ORIGIN_EXPECTED directly
  assert_eq!(LOCAL_IJK_END_ORIGIN_EXPECTED, H3_IJK_PLUS_N3_RUST);
}

#[test]
fn test_step6_local_ijk_to_cell_rust_path() {
  let mut h3_out = H3_NULL;
  assert!(local_ijk_to_cell(START_H3_FAIL_CASE, &H3_IJK_PLUS_N0_RUST, &mut h3_out).is_ok());
  assert_eq!(h3_out, PATH_N0_RUST);

  assert!(local_ijk_to_cell(START_H3_FAIL_CASE, &H3_IJK_PLUS_N1_RUST, &mut h3_out).is_ok());
  assert_eq!(h3_out, PATH_N1_RUST);

  assert!(local_ijk_to_cell(START_H3_FAIL_CASE, &H3_IJK_PLUS_N2_RUST, &mut h3_out).is_ok());
  assert_eq!(h3_out, PATH_N2_RUST);

  assert!(local_ijk_to_cell(START_H3_FAIL_CASE, &H3_IJK_PLUS_N3_RUST, &mut h3_out).is_ok());
  assert_eq!(h3_out, PATH_N3_RUST);
}

#[test]
fn test_step7_compare_rust_rounded_axial_with_c_implied_axial() {
  // Rust's rounded axial for n=1:
  let rust_rounded_axial_n1 = CoordIJK { i: -1, j: -1, k: 2 };
  // Rust's rounded axial for n=2:
  let rust_rounded_axial_n2 = CoordIJK { i: -2, j: -2, k: 4 };

  // Get local IJK+ for C path cells
  let mut local_ijk_c1 = CoordIJK::default();
  assert!(cell_to_local_ijk(START_H3_FAIL_CASE, C_PATH_N1, &mut local_ijk_c1).is_ok());
  // Convert this local IJK+ to axial to see what C's _cubeRound must have produced
  let mut c_implied_axial_n1 = local_ijk_c1;
  ijk_to_cube(&mut c_implied_axial_n1);
  eprintln!(
    "C path[1] (0x{:x}) -> local IJK+ {:?} -> implied axial {:?}",
    C_PATH_N1.0, local_ijk_c1, c_implied_axial_n1
  );

  let mut local_ijk_c2 = CoordIJK::default();
  assert!(cell_to_local_ijk(START_H3_FAIL_CASE, C_PATH_N2, &mut local_ijk_c2).is_ok());
  let mut c_implied_axial_n2 = local_ijk_c2;
  ijk_to_cube(&mut c_implied_axial_n2);
  eprintln!(
    "C path[2] (0x{:x}) -> local IJK+ {:?} -> implied axial {:?}",
    C_PATH_N2.0, local_ijk_c2, c_implied_axial_n2
  );

  // THE CRITICAL ASSERTIONS:
  assert_eq!(
    rust_rounded_axial_n1, c_implied_axial_n1,
    "Mismatch in rounded axial coords for n=1. Rust produced {:?}, C path implies {:?}",
    rust_rounded_axial_n1, c_implied_axial_n1
  );
  assert_eq!(
    rust_rounded_axial_n2, c_implied_axial_n2,
    "Mismatch in rounded axial coords for n=2. Rust produced {:?}, C path implies {:?}",
    rust_rounded_axial_n2, c_implied_axial_n2
  );
}

#[test]
fn test_grid_path_internals_for_failing_case() {
  let start_h3 = H3Index(0x855943cbfffffff);
  let end_h3 = H3Index(0x85594303fffffff);

  // 1. cell_to_local_ijk Verification
  let mut local_ijk_start_origin = CoordIJK::default();
  assert!(cell_to_local_ijk(start_h3, start_h3, &mut local_ijk_start_origin).is_ok());
  assert_eq!(local_ijk_start_origin, CoordIJK { i: 0, j: 0, k: 0 });

  let mut local_ijk_end_origin = CoordIJK::default();
  assert!(cell_to_local_ijk(start_h3, end_h3, &mut local_ijk_end_origin).is_ok());
  assert_eq!(local_ijk_end_origin, CoordIJK { i: -3, j: -3, k: 0 });

  // 2. ijk_to_cube Verification
  let mut axial_start = local_ijk_start_origin;
  ijk_to_cube(&mut axial_start);
  assert_eq!(axial_start, CoordIJK { i: 0, j: 0, k: 0 });

  let mut axial_end = local_ijk_end_origin;
  ijk_to_cube(&mut axial_end);
  assert_eq!(axial_end, CoordIJK { i: -3, j: -3, k: 6 });

  // 3. Line Interpolation Step Values
  let distance = grid_distance(start_h3, end_h3).unwrap();
  assert_eq!(distance, 3);
  let inv_distance_f = 1.0 / (distance as f64);
  let i_step = (axial_end.i - axial_start.i) as f64 * inv_distance_f;
  let j_step = (axial_end.j - axial_start.j) as f64 * inv_distance_f;
  let k_step = (axial_end.k - axial_start.k) as f64 * inv_distance_f;
  assert_eq!(i_step, -1.0);
  assert_eq!(j_step, -1.0);
  assert_eq!(k_step, 2.0);

  // 4. ijk_round_to_axial_hex_center Verification
  let mut rounded_axial_n0 = CoordIJK::default();
  ijk_round_to_axial_hex_center_test_version(0.0, 0.0, 0.0, &mut rounded_axial_n0);
  assert_eq!(rounded_axial_n0, CoordIJK { i: 0, j: 0, k: 0 });

  let mut rounded_axial_n1 = CoordIJK::default();
  ijk_round_to_axial_hex_center_test_version(-1.0, -1.0, 2.0, &mut rounded_axial_n1);
  assert_eq!(rounded_axial_n1, CoordIJK { i: -1, j: -1, k: 2 });

  let mut rounded_axial_n2 = CoordIJK::default();
  ijk_round_to_axial_hex_center_test_version(-2.0, -2.0, 4.0, &mut rounded_axial_n2);
  assert_eq!(rounded_axial_n2, CoordIJK { i: -2, j: -2, k: 4 });

  let mut rounded_axial_n3 = CoordIJK::default();
  ijk_round_to_axial_hex_center_test_version(-3.0, -3.0, 6.0, &mut rounded_axial_n3);
  assert_eq!(rounded_axial_n3, CoordIJK { i: -3, j: -3, k: 6 });

  // 5. Axial to H3 IJK+ Conversion
  let mut h3_ijk_plus_n0 = rounded_axial_n0;
  h3_ijk_plus_n0.i = rounded_axial_n0.i.saturating_neg();
  h3_ijk_plus_n0.j = rounded_axial_n0.j;
  h3_ijk_plus_n0.k = 0;
  _ijk_normalize(&mut h3_ijk_plus_n0);
  assert_eq!(h3_ijk_plus_n0, CoordIJK { i: 0, j: 0, k: 0 });

  let mut h3_ijk_plus_n1 = rounded_axial_n1;
  h3_ijk_plus_n1.i = rounded_axial_n1.i.saturating_neg(); // -(-1) = 1
  h3_ijk_plus_n1.j = rounded_axial_n1.j; // -1
  h3_ijk_plus_n1.k = 0; // 0
                        // normalize {1, -1, 0} -> {2,0,1}
  _ijk_normalize(&mut h3_ijk_plus_n1);
  assert_eq!(h3_ijk_plus_n1, CoordIJK { i: 2, j: 0, k: 1 });

  let mut h3_ijk_plus_n2 = rounded_axial_n2;
  h3_ijk_plus_n2.i = rounded_axial_n2.i.saturating_neg(); // -(-2) = 2
  h3_ijk_plus_n2.j = rounded_axial_n2.j; // -2
  h3_ijk_plus_n2.k = 0; // 0
                        // normalize {2, -2, 0} -> {4,0,2}
  _ijk_normalize(&mut h3_ijk_plus_n2);
  assert_eq!(h3_ijk_plus_n2, CoordIJK { i: 4, j: 0, k: 2 });

  // For n=3 (end point), grid_path_cells uses local_ijk_end_origin directly:
  let h3_ijk_plus_n3_direct = local_ijk_end_origin;
  assert_eq!(h3_ijk_plus_n3_direct, CoordIJK { i: -3, j: -3, k: 0 });

  // 6. local_ijk_to_cell Verification (matches Rust path)
  let mut path_n0 = H3_NULL;
  assert!(local_ijk_to_cell(start_h3, &h3_ijk_plus_n0, &mut path_n0).is_ok());
  assert_eq!(path_n0, H3Index(0x855943cbfffffff)); // start

  let mut path_n1 = H3_NULL;
  assert!(local_ijk_to_cell(start_h3, &h3_ijk_plus_n1, &mut path_n1).is_ok());
  assert_eq!(path_n1, H3Index(0x855943d3fffffff)); // rust path[1]

  let mut path_n2 = H3_NULL;
  assert!(local_ijk_to_cell(start_h3, &h3_ijk_plus_n2, &mut path_n2).is_ok());
  assert_eq!(path_n2, H3Index(0x85594063fffffff)); // rust path[2]

  let mut path_n3 = H3_NULL;
  assert!(local_ijk_to_cell(start_h3, &h3_ijk_plus_n3_direct, &mut path_n3).is_ok());
  assert_eq!(path_n3, H3Index(0x85594303fffffff)); // end

  // 7. Compare with C Path Intermediate Local IJKs
  let c_path1_h3 = H3Index(0x855943cffffffff);
  let c_path2_h3 = H3Index(0x8559431bfffffff);

  let mut local_ijk_c1 = CoordIJK::default();
  assert!(cell_to_local_ijk(start_h3, c_path1_h3, &mut local_ijk_c1).is_ok());
  eprintln!(
    "Local IJK+ for C_path[1] (0x{:x}) relative to start: {:?}",
    c_path1_h3.0, local_ijk_c1
  );
  // Expected local IJK+ for C_path1 should be a distance 1 vector, e.g., {1,0,0} or similar after norm.
  // From the C output of h3Distance for start and C_path[1], distance is 1.
  // If C_path[1] is {i:X, j:Y, k:Z} (axial from its _cubeRound), then after conversion to local IJK+
  // for localIjkToCell, this local IJK+ should be `local_ijk_c1`.

  let mut local_ijk_c2 = CoordIJK::default();
  assert!(cell_to_local_ijk(start_h3, c_path2_h3, &mut local_ijk_c2).is_ok());
  eprintln!(
    "Local IJK+ for C_path[2] (0x{:x}) relative to start: {:?}",
    c_path2_h3.0, local_ijk_c2
  );

  // Convert these to axial to see what C's _cubeRound *must have* produced for n=1 and n=2
  let mut axial_from_local_c1 = local_ijk_c1;
  ijk_to_cube(&mut axial_from_local_c1);
  eprintln!("Axial for C_path[1]'s local IJK+: {:?}", axial_from_local_c1);

  let mut axial_from_local_c2 = local_ijk_c2;
  ijk_to_cube(&mut axial_from_local_c2);
  eprintln!("Axial for C_path[2]'s local IJK+: {:?}", axial_from_local_c2);

  // Now, compare `axial_from_local_c1` with Rust's `rounded_axial_n1` (which was {-1,-1,2})
  // And `axial_from_local_c2` with Rust's `rounded_axial_n2` (which was {-2,-2,4})
  // If they are different, then _cubeRound (ijk_round_to_axial_hex_center_test) is the point of divergence.

  assert_eq!(
    axial_from_local_c1, rounded_axial_n1,
    "Mismatch in rounded axial coords for n=1. C path implies axial {:?}, Rust produced {:?}",
    axial_from_local_c1, rounded_axial_n1
  );
  assert_eq!(
    axial_from_local_c2, rounded_axial_n2,
    "Mismatch in rounded axial coords for n=2. C path implies axial {:?}, Rust produced {:?}",
    axial_from_local_c2, rounded_axial_n2
  );
}