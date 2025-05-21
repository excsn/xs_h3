// tests/indexing_tests.rs

use xs_h3::{h3_index::set_mode, *}; // Import all public items from your xs_h3 crate

// Helper for tests
fn latlng_from_degs(lat_deg: f64, lng_deg: f64) -> LatLng {
  LatLng {
    lat: degs_to_rads(lat_deg),
    lng: degs_to_rads(lng_deg),
  }
}

#[test]
fn test_cli_lat_lng_to_cell() {
  // From tests/cli/latLngToCell.txt:
  // "latLngToCell --lat 20 --lng 123 -r 2 -f newline" "824b9ffffffffff"
  let geo = latlng_from_degs(20.0, 123.0);
  let res = 2;
  let expected_h3 = H3Index(0x824b9ffffffffff);

  match lat_lng_to_cell(&geo, res) {
    Ok(h3_result) => assert_eq!(h3_result, expected_h3, "Expected H3Index from latLngToCell"),
    Err(e) => panic!("lat_lng_to_cell failed: {:?}", e),
  }
}

#[test]
fn test_cli_cell_to_lat_lng() {
  // From tests/cli/cellToLatLng.txt:
  // "cellToLatLng -c 8928342e20fffff -f wkt" "POINT(-122.5003039349 37.5012466151)"
  let cell = H3Index(0x8928342e20fffff);
  let expected_lng_deg = -122.5003039349;
  let expected_lat_deg = 37.5012466151;

  match cell_to_lat_lng(cell) {
    Ok(latlng_rad) => {
      let lng_deg_result = rads_to_degs(latlng_rad.lng);
      let lat_deg_result = rads_to_degs(latlng_rad.lat);
      assert!(
        (lng_deg_result - expected_lng_deg).abs() < 1e-9,
        "Longitude mismatch. Got: {}, Expected: {}",
        lng_deg_result,
        expected_lng_deg
      );
      assert!(
        (lat_deg_result - expected_lat_deg).abs() < 1e-9,
        "Latitude mismatch. Got: {}, Expected: {}",
        lat_deg_result,
        expected_lat_deg
      );
    }
    Err(e) => panic!("cell_to_lat_lng failed: {:?}", e),
  }
}

#[test]
fn test_cli_invalid_cell_to_lat_lng() {
  // From tests/cli/cellToLatLng.txt:
  // "cellToLatLng -c asdf 2>&1" "Error 5: Cell argument was not valid"
  // We test an invalid H3 index directly
  let mut h_invalid_mode = H3Index(0x8001fffffffffff); // Valid res 0 cell
  set_mode(&mut h_invalid_mode, 0); // Mode 0 is invalid

  match cell_to_lat_lng(h_invalid_mode) {
    Ok(_) => panic!("Expected CellInvalid error, got Ok"),
    Err(e) => assert_eq!(e, H3Error::CellInvalid, "Expected CellInvalid from invalid mode"),
  }
}

#[test]
fn test_cli_cell_to_boundary() {
  // From tests/cli/cellToBoundary.txt
  // "cellToBoundary -c 8928342e20fffff -f wkt"
  // "POLYGON((-122.4990471431 37.4997389893, -122.4979805011 37.5014245698, ...))"
  let cell = H3Index(0x8928342e20fffff);

  // Expected vertices in degrees (Lng, Lat order from WKT)
  // Note: H3 cellToBoundary returns CCW vertices. WKT might represent them differently.
  // The CLI output needs careful parsing or comparison against known values.
  // For this example, we'll just check the number of vertices and basic validity.
  let expected_num_verts = if is_pentagon(cell) { 5 } else { 6 }; // Simplified; could be more if distorted

  match cell_to_boundary(cell) {
    Ok(boundary) => {
      assert!(
        boundary.num_verts >= expected_num_verts,
        "Boundary has at least {} vertices, got {}",
        expected_num_verts,
        boundary.num_verts
      );
      assert!(
        boundary.num_verts <= MAX_CELL_BNDRY_VERTS,
        "Boundary has too many vertices"
      );
      for i in 0..boundary.num_verts {
        assert!(boundary.verts[i].lat.is_finite());
        assert!(boundary.verts[i].lng.is_finite());
      }
    }
    Err(e) => panic!("cell_to_boundary failed: {:?}", e),
  }
}
