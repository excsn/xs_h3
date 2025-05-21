// src/traversal/distance.rs

use crate::coords::ijk::ijk_distance; // Assuming ijk_distance is pub(crate) in coords::ijk
use crate::h3_index::inspection::is_valid_cell;
use crate::h3_index::{get_mode, get_resolution};
use crate::local_ij::cell_to_local_ijk; // Assuming cellToLocalIjk is here and pub(crate)
use crate::types::{H3Error, H3Index};

/// Produces the grid distance between the two H3 cells.
///
/// This function may fail to find the distance between two cells, for
/// example if they are very far apart or on opposite sides of a pentagon.
///
/// # Arguments
/// * `origin` - The origin H3 cell.
/// * `destination` - The destination H3 cell.
///
/// # Returns
/// `Ok(distance)` on success, or an `H3Error` if the distance cannot be computed.
pub fn grid_distance(origin: H3Index, destination: H3Index) -> Result<i64, H3Error> {
  if get_mode(origin) != crate::constants::H3_CELL_MODE || get_mode(destination) != crate::constants::H3_CELL_MODE {
    return Err(H3Error::CellInvalid); // Or a more specific mode error
  }

  if get_resolution(origin) != get_resolution(destination) {
    return Err(H3Error::ResMismatch);
  }

  // Validate cells more thoroughly if a public API
  if !is_valid_cell(origin) || !is_valid_cell(destination) {
    return Err(H3Error::CellInvalid);
  }

  // This is the core logic from C's gridDistance:
  // Convert to local IJK coordinates centered on the origin.
  let mut origin_ijk = crate::types::CoordIJK::default();
  // The IJK for origin relative to itself is always (0,0,0)
  // cell_to_local_ijk(origin, origin, &mut origin_ijk)?; // This should yield (0,0,0) or error
  // More directly, the origin IJK in its own local coordinate system is (0,0,0)
  // if the subsequent call to cell_to_local_ijk for destination uses `origin` as its anchor.
  // The C code: `_h3ToLocalIJK(origin, origin, &originIjk)` and `_h3ToLocalIJK(origin, h3, &h3Ijk)`.
  // The first call is effectively setting originIjk to {0,0,0} in its own frame.
  // We can skip it if ijk_distance is robust to that or if we assume that.
  // For safety, let's do what C does:
  cell_to_local_ijk(origin, origin, &mut origin_ijk)?;
  // This should result in origin_ijk = {0,0,0} if origin is valid.

  let mut destination_ijk = crate::types::CoordIJK::default();
  cell_to_local_ijk(origin, destination, &mut destination_ijk)?;

  // The distance is the IJK distance between these two points.
  Ok(ijk_distance(&origin_ijk, &destination_ijk) as i64)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::indexing::lat_lng_to_cell;
  use crate::latlng::{_set_geo_degs};
  use crate::traversal::grid_disk::grid_disk;
  use crate::types::{H3_NULL, LatLng}; // For getting neighbors

  #[test]
  fn test_grid_distance_identity() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let h = lat_lng_to_cell(&geo, 5).unwrap();
    assert_eq!(grid_distance(h, h), Ok(0));
  }

  #[test]
  fn test_grid_distance_direct_neighbors() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let origin = lat_lng_to_cell(&geo, 5).unwrap();

    let mut ring = [H3_NULL; 7];
    grid_disk(origin, 1, &mut ring).unwrap(); // Assumes grid_disk is available

    for neighbor_h3 in ring {
      if neighbor_h3 != H3_NULL && neighbor_h3 != origin {
        assert_eq!(
          grid_distance(origin, neighbor_h3),
          Ok(1),
          "Direct neighbor distance is 1. Neighbor: {:x}",
          neighbor_h3.0
        );
      }
    }
  }

  #[test]
  fn test_grid_distance_res_mismatch() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let h_res5 = lat_lng_to_cell(&geo, 5).unwrap();
    let h_res6 = lat_lng_to_cell(&geo, 6).unwrap();
    assert_eq!(grid_distance(h_res5, h_res6), Err(H3Error::ResMismatch));
  }

  #[test]
  fn test_grid_distance_invalid_input() {
    let valid_h = lat_lng_to_cell(&LatLng { lat: 0.0, lng: 0.0 }, 5).unwrap();
    assert_eq!(grid_distance(H3_NULL, valid_h), Err(H3Error::CellInvalid));
    assert_eq!(grid_distance(valid_h, H3_NULL), Err(H3Error::CellInvalid));

    let mut invalid_mode_h = valid_h;
    crate::h3_index::set_mode(&mut invalid_mode_h, crate::constants::H3_DIRECTEDEDGE_MODE as u8);
    assert_eq!(grid_distance(invalid_mode_h, valid_h), Err(H3Error::CellInvalid));
  }

  // More tests are needed for pentagon involvement and long distances,
  // especially those that might fail in cellToLocalIjk.
  // The C tests for gridDistance (e.g., testGridDistance.c) are good sources.
}
