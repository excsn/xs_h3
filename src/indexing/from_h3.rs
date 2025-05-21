// src/indexing/from_h3.rs

use crate::coords::face_ijk::{_face_ijk_pent_to_cell_boundary, _face_ijk_to_cell_boundary, _face_ijk_to_geo};
use crate::h3_index::inspection::{is_pentagon, is_valid_cell};
use crate::h3_index::{_h3_to_face_ijk, get_resolution};
use crate::types::{CellBoundary, FaceIJK, H3Error, H3Index, LatLng}; // From inspection submodule

/// Finds the center `LatLng` point of the given H3 cell.
///
/// # Arguments
/// * `cell` - The H3 cell.
///
/// # Returns
/// `Ok(LatLng)` on success, or an `H3Error` if the cell is invalid.
pub fn cell_to_lat_lng(cell: H3Index) -> Result<LatLng, H3Error> {
  if !is_valid_cell(cell) {
    // Using the public is_valid_cell
    return Err(H3Error::CellInvalid);
  }

  let mut fijk = FaceIJK::default();
  // _h3_to_face_ijk can fail if the cell, despite passing basic isValidCell,
  // has an invalid internal structure (e.g., bad base cell for digits)
  // that only _h3ToFaceIjk's deeper logic catches.
  _h3_to_face_ijk(cell, &mut fijk)?; // Propagate error if any

  let mut geo = LatLng::default();
  _face_ijk_to_geo(&fijk, get_resolution(cell), &mut geo);
  Ok(geo)
}

/// Finds the boundary of the given H3 cell.
///
/// # Arguments
/// * `cell` - The H3 cell.
///
/// # Returns
/// `Ok(CellBoundary)` on success, or an `H3Error` if the cell is invalid.
pub fn cell_to_boundary(cell: H3Index) -> Result<CellBoundary, H3Error> {
  if !is_valid_cell(cell) {
    return Err(H3Error::CellInvalid);
  }

  let mut fijk = FaceIJK::default();
  _h3_to_face_ijk(cell, &mut fijk)?; // Propagate error

  let mut boundary = CellBoundary::default();
  let res = get_resolution(cell);

  if is_pentagon(cell) {
    // Using the public is_pentagon
    _face_ijk_pent_to_cell_boundary(&fijk, res, 0, crate::constants::NUM_PENT_VERTS as i32, &mut boundary);
  } else {
    _face_ijk_to_cell_boundary(&fijk, res, 0, crate::constants::NUM_HEX_VERTS as i32, &mut boundary);
  }
  Ok(boundary)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::constants::MAX_CELL_BNDRY_VERTS;
  use crate::h3_index::inspection::is_valid_cell;
  use crate::indexing::to_h3::lat_lng_to_cell; // For setting up test cells
  use crate::latlng::_set_geo_degs;
  use crate::types::{H3_NULL};

  #[test]
  fn test_cell_to_lat_lng_invalid_input() {
    assert_eq!(cell_to_lat_lng(H3_NULL), Err(H3Error::CellInvalid));
    let mut invalid_h = H3Index(0x85283473fffffff); // A valid cell
    crate::h3_index::set_mode(&mut invalid_h, crate::constants::H3_DIRECTEDEDGE_MODE as u8);
    assert_eq!(cell_to_lat_lng(invalid_h), Err(H3Error::CellInvalid));
  }

  #[test]
  fn test_cell_to_boundary_invalid_input() {
    assert_eq!(cell_to_boundary(H3_NULL).map_err(|e| e), Err(H3Error::CellInvalid));
    let mut invalid_h = H3Index(0x85283473fffffff); // A valid cell
    crate::h3_index::set_mode(&mut invalid_h, crate::constants::H3_DIRECTEDEDGE_MODE as u8);
    assert_eq!(cell_to_boundary(invalid_h).map_err(|e| e), Err(H3Error::CellInvalid));
  }

  #[test]
  fn test_cell_to_lat_lng_and_boundary_roundtrip() {
    // Test a few resolutions for a known point
    let mut geo_orig = LatLng::default();
    _set_geo_degs(&mut geo_orig, 37.779, -122.419); // Near SF

    for res in 0..=10 {
      // Test up to res 10
      let cell = lat_lng_to_cell(&geo_orig, res).unwrap();
      assert!(is_valid_cell(cell));

      // Test cellToLatLng
      let center_geo = cell_to_lat_lng(cell).unwrap();
      // The center of a cell might not be the original point, but it should be close
      // and re-indexing it should yield the same cell.
      let cell_from_center = lat_lng_to_cell(&center_geo, res).unwrap();
      assert_eq!(
        cell_from_center, cell,
        "Center of cell should re-index to itself at res {}",
        res
      );

      // Test cellToBoundary
      let boundary = cell_to_boundary(cell).unwrap();
      if crate::h3_index::inspection::is_pentagon(cell) {
        assert!(
          boundary.num_verts >= 5 && boundary.num_verts <= MAX_CELL_BNDRY_VERTS,
          "Pentagon boundary has between 5 and {} verts, got {} at res {}",
          MAX_CELL_BNDRY_VERTS,
          boundary.num_verts,
          res
        );
      } else {
        assert!(
          boundary.num_verts >= 6 && boundary.num_verts <= MAX_CELL_BNDRY_VERTS,
          "Hexagon boundary has between 6 and {} verts, got {} at res {}",
          MAX_CELL_BNDRY_VERTS,
          boundary.num_verts,
          res
        );
      }

      // Check that the original point is (likely) within the boundary produced
      // This is a weaker check as point-in-poly is complex.
      // A simple check: is the original point closer to this cell's center than to any other cell's center?
      // This is implicitly tested by lat_lng_to_cell.

      // Check that all boundary vertices are valid lat/lngs
      for i in 0..boundary.num_verts {
        assert!(
          boundary.verts[i].lat.is_finite() && boundary.verts[i].lng.is_finite(),
          "Boundary vertex {} is not finite for cell {:x} res {}",
          i,
          cell.0,
          res
        );
        assert!(
          boundary.verts[i].lat.abs() <= crate::constants::M_PI_2 + crate::constants::EPSILON_RAD,
          "Boundary vertex {} lat out of range for cell {:x} res {}",
          i,
          cell.0,
          res
        );
      }
    }
  }
}
