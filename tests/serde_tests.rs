// tests/serde_tests.rs

// Only compile and run these tests if the "serde" feature is enabled.
#![cfg(feature = "serde")]

use serde_json;
use xs_h3::{types::Direction, *}; // Import all public items from your xs_h3 crate // For easy JSON serialization/deserialization in tests

#[test]
fn test_h3index_serde() {
  let h = H3Index(0x8928308280fffff);
  let serialized = serde_json::to_string(&h).unwrap();
  // H3Index is repr(transparent) over u64, so it serializes as its inner u64 value directly.
  assert_eq!(serialized, "617700169958293503");
  let deserialized: H3Index = serde_json::from_str(&serialized).unwrap();
  assert_eq!(h, deserialized);

  let h_null = H3_NULL;
  let ser_null = serde_json::to_string(&h_null).unwrap();
  assert_eq!(ser_null, "0");
  let de_null: H3Index = serde_json::from_str(&ser_null).unwrap();
  assert_eq!(h_null, de_null);
}

#[test]
fn test_latlng_serde() {
  let ll = LatLng { lat: 0.5, lng: -1.2 };
  let serialized = serde_json::to_string(&ll).unwrap();
  assert_eq!(serialized, r#"{"lat":0.5,"lng":-1.2}"#); // Note the raw string for JSON
  let deserialized: LatLng = serde_json::from_str(&serialized).unwrap();
  assert_eq!(ll, deserialized);
}

#[test]
fn test_h3error_serde() {
  let err = H3Error::CellInvalid; // This has repr(u32) and value 5
  let serialized = serde_json::to_string(&err).unwrap();
  assert_eq!(serialized, "5"); // serde_repr serializes to the number
  let deserialized: H3Error = serde_json::from_str(&serialized).unwrap();
  assert_eq!(err, deserialized);

  let err_success = H3Error::Success;
  let ser_success = serde_json::to_string(&err_success).unwrap();
  assert_eq!(ser_success, "0");
  let de_success: H3Error = serde_json::from_str(&ser_success).unwrap();
  assert_eq!(err_success, de_success);
}

#[test]
fn test_direction_serde() {
  let dir = Direction::KAxes; // This has repr(u8) and value 1
  let serialized = serde_json::to_string(&dir).unwrap();
  assert_eq!(serialized, "1"); // serde_repr serializes to the number
  let deserialized: Direction = serde_json::from_str(&serialized).unwrap();
  assert_eq!(dir, deserialized);
}

#[test]
fn test_geopolygon_serde() {
  let outer_loop = GeoLoop {
    num_verts: 2,
    verts: vec![LatLng { lat: 1.0, lng: 1.0 }, LatLng { lat: 2.0, lng: 2.0 }],
  };
  let hole_loop = GeoLoop {
    num_verts: 1,
    verts: vec![LatLng { lat: 1.5, lng: 1.5 }],
  };
  let polygon = GeoPolygon {
    geoloop: outer_loop.clone(), // Clone if outer_loop is used elsewhere
    num_holes: 1,
    holes: vec![hole_loop.clone()], // Clone if hole_loop is used elsewhere
  };

  let serialized = serde_json::to_string(&polygon).unwrap();
  // Example: {"geoloop":{"num_verts":2,"verts":[{"lat":1.0,"lng":1.0},{"lat":2.0,"lng":2.0}]},"num_holes":1,"holes":[{"num_verts":1,"verts":[{"lat":1.5,"lng":1.5}]}]}
  // println!("Serialized GeoPolygon: {}", serialized);

  let deserialized: GeoPolygon = serde_json::from_str(&serialized).unwrap();
  assert_eq!(polygon, deserialized);
}

#[test]
fn test_cell_boundary_serde() {
  let mut cb = CellBoundary::default();
  cb.num_verts = 2;
  cb.verts[0] = LatLng { lat: 1.0, lng: 1.0 };
  cb.verts[1] = LatLng { lat: 2.0, lng: 2.0 };
  // cb.verts[2..] are LatLng::default()

  let serialized = serde_json::to_string(&cb).unwrap();
  // println!("Serialized CellBoundary: {}", serialized);

  let deserialized: CellBoundary = serde_json::from_str(&serialized).unwrap();
  assert_eq!(cb, deserialized);
}

#[test]
fn test_vec_h3index_serde() {
  let indices = vec![
    H3Index(0x8928308280fffff),
    H3Index(0x8928308281fffff),
    H3_NULL, // H3Index(0)
  ];
  let serialized = serde_json::to_string(&indices).unwrap();
  assert_eq!(serialized, "[617700169958293503,617700169959342079,0]");
  let deserialized: Vec<H3Index> = serde_json::from_str(&serialized).unwrap();
  assert_eq!(indices, deserialized);
}

#[test]
fn test_polygon_rust_serde() {
  // Ensure PolygonRust is accessible, e.g. by `use xs_h3::regions::PolygonRust;`
  // if it's re-exported from lib.rs, or direct path.
  // If it's pub(crate) in regions::to_polygon, this test needs to be in that module.
  // Assuming it's made public for this test as per notes.
  use xs_h3::PolygonRust; // Assuming re-exported or made public

  let poly_rust = PolygonRust {
    outer: vec![LatLng { lat: 1.0, lng: 1.0 }, LatLng { lat: 2.0, lng: 2.0 }],
    holes: vec![vec![LatLng { lat: 1.5, lng: 1.5 }]],
  };
  let serialized = serde_json::to_string(&poly_rust).unwrap();
  let deserialized: PolygonRust = serde_json::from_str(&serialized).unwrap();

  assert_eq!(poly_rust.outer, deserialized.outer);
  assert_eq!(poly_rust.holes, deserialized.holes);
}

#[test]
fn test_multi_polygon_rust_serde() {
  use xs_h3::{MultiPolygonRust, PolygonRust}; // Assuming re-exported

  let poly1 = PolygonRust {
    outer: vec![LatLng { lat: 1.0, lng: 1.0 }, LatLng { lat: 2.0, lng: 2.0 }],
    holes: vec![],
  };
  let poly2 = PolygonRust {
    outer: vec![LatLng { lat: 3.0, lng: 3.0 }, LatLng { lat: 4.0, lng: 4.0 }],
    holes: vec![vec![LatLng { lat: 3.5, lng: 3.5 }]],
  };
  let multi_poly: MultiPolygonRust = vec![poly1.clone(), poly2.clone()];

  let serialized = serde_json::to_string(&multi_poly).unwrap();
  let deserialized: MultiPolygonRust = serde_json::from_str(&serialized).unwrap();

  assert_eq!(multi_poly.len(), deserialized.len());
  assert_eq!(multi_poly[0].outer, deserialized[0].outer);
  assert_eq!(multi_poly[1].outer, deserialized[1].outer);
  assert_eq!(multi_poly[1].holes, deserialized[1].holes);
}
