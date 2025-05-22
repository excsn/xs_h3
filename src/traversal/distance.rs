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
    // Initial debug print for the function call
    eprintln!(
        "DEBUG: grid_distance called with origin=0x{:x}, destination=0x{:x}",
        origin.0, destination.0
    );

    // Initial checks from C version of gridDistance (before _h3Distance)
    if get_mode(origin) != crate::constants::H3_CELL_MODE
        || get_mode(destination) != crate::constants::H3_CELL_MODE
    {
        eprintln!("  DEBUG: grid_distance -> returning CellInvalid (mode check)");
        return Err(H3Error::CellInvalid);
    }

    if get_resolution(origin) != get_resolution(destination) {
        eprintln!("  DEBUG: grid_distance -> returning ResMismatch");
        return Err(H3Error::ResMismatch);
    }

    // The C public API H3_EXPORT(gridDistance) also calls isValidCell.
    if !is_valid_cell(origin) || !is_valid_cell(destination) {
        eprintln!("  DEBUG: grid_distance -> returning CellInvalid (isValidCell check)");
        return Err(H3Error::CellInvalid);
    }
    
    // If origin and destination are the same, distance is 0.
    if origin == destination {
        eprintln!("  DEBUG: grid_distance -> origin == destination, returning Ok(0)");
        return Ok(0);
    }

    // Get local IJK coordinates.
    // origin_ijk is implicitly {0,0,0} in the origin's local coordinate system.
    let origin_ijk_in_own_frame = crate::types::CoordIJK { i: 0, j: 0, k: 0 };

    let mut destination_ijk_in_origin_frame = crate::types::CoordIJK::default();
    // `cell_to_local_ijk` will have its own eprintln! for inputs and outputs
    match cell_to_local_ijk(origin, destination, &mut destination_ijk_in_origin_frame) {
        Ok(()) => {
            eprintln!(
                "  DEBUG: grid_distance: cell_to_local_ijk for destination (0x{:x} relative to 0x{:x}) yielded: {:?}",
                destination.0, origin.0, destination_ijk_in_origin_frame
            );
        }
        Err(e) => {
            eprintln!(
                "  DEBUG: grid_distance: cell_to_local_ijk for destination (0x{:x} relative to 0x{:x}) errored: {:?}",
                destination.0, origin.0, e
            );
            return Err(e); // Propagate error from cell_to_local_ijk
        }
    }

    // The distance is the IJK grid distance between {0,0,0} and the destination's local IJK.
    let dist_val = ijk_distance(&origin_ijk_in_own_frame, &destination_ijk_in_origin_frame) as i64;
    
    eprintln!(
        "  DEBUG: grid_distance: ijk_distance({:?}, {:?}) calculated distance = {}",
        origin_ijk_in_own_frame, destination_ijk_in_origin_frame, dist_val
    );
    eprintln!(
        "DEBUG: grid_distance for origin=0x{:x}, destination=0x{:x} FINISHED, returning Ok({})",
        origin.0, destination.0, dist_val
    );
    Ok(dist_val)
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
