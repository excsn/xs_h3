// src/traversal/grid_path.rs

use crate::coords::ijk::{cube_to_ijk, ijk_to_cube}; // Assuming these are pub(crate)
use crate::h3_index::get_mode;
use crate::h3_index::inspection::is_valid_cell;
use crate::local_ij::cell_to_local_ijk; // Assuming cellToLocalIjk is here
use crate::local_ij::local_ijk_to_cell; // Assuming localIjkToCell is here
use crate::traversal::distance::grid_distance; // From the distance module
use crate::types::{H3Error, H3Index, CoordIJK};

/// Number of H3 cells in a line from the `start` H3 cell to the `end` H3 cell.
///
/// # Arguments
/// * `start` - The start H3 cell.
/// * `end` - The end H3 cell.
///
/// # Returns
/// `Ok(size)` of the line, or an `H3Error` if the line cannot be computed.
pub fn grid_path_cells_size(start: H3Index, end: H3Index) -> Result<i64, H3Error> {
  match grid_distance(start, end) {
    Ok(dist) => Ok(dist + 1),
    Err(e) => Err(e), // Propagate error from grid_distance
  }
}

/// Helper: Rounds cube coordinates to the nearest integer cube coordinate.
/// Algorithm from https://www.redblobgames.com/grids/hexagons/#rounding
fn cube_round(i_f: f64, j_f: f64, k_f: f64, out_ijk: &mut CoordIJK) {
  let mut ri = i_f.round();
  let mut rj = j_f.round();
  let mut rk = k_f.round();

  let i_diff = (ri - i_f).abs();
  let j_diff = (rj - j_f).abs();
  let k_diff = (rk - k_f).abs();

  // Round, maintaining valid cube coords (sum of components is 0)
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

/// Given two H3 cells, returns the line of H3 cells between them (inclusive).
///
/// This function may fail to find the line between two cells, for
/// example if they are very far apart or on opposite sides of a pentagon.
///
/// Notes:
///  - The specific output of this function should not be considered stable
///    across library versions. The only guarantees are line length and that
///    each cell in the line is a neighbor of the preceding cell.
///  - Lines are drawn in grid space, and may not correspond exactly to
///    Cartesian lines or great arcs.
///
/// # Arguments
/// * `start` - The start H3 cell.
/// * `end` - The end H3 cell.
/// * `out_path` - Output array for the H3 cells in the path. Must be
///                pre-allocated to the size determined by `grid_path_cells_size`.
///
/// # Returns
/// `Ok(())` on success, or an `H3Error` if the path cannot be computed or
/// `out_path` is too small.
pub fn grid_path_cells(start: H3Index, end: H3Index, out_path: &mut [H3Index]) -> Result<(), H3Error> {
  let distance = match grid_distance(start, end) {
    Ok(dist) => dist,
    Err(e) => return Err(e),
  };

  let expected_size = distance + 1;
  if out_path.len() < expected_size as usize {
    return Err(H3Error::MemoryBounds);
  }

  // Get IJK coords for the start and end.
  // cell_to_local_ijk uses `start` as the anchor for the coordinate system.
  let mut start_ijk_local = crate::types::CoordIJK::default(); // This will be {0,0,0} in start's frame
  cell_to_local_ijk(start, start, &mut start_ijk_local)?;

  let mut end_ijk_local = crate::types::CoordIJK::default();
  cell_to_local_ijk(start, end, &mut end_ijk_local)?;

  // Convert local IJK to cube coordinates for linear interpolation.
  let mut start_cube = start_ijk_local;
  ijk_to_cube(&mut start_cube);
  let mut end_cube = end_ijk_local;
  ijk_to_cube(&mut end_cube);

  let inv_distance_f = if distance > 0 { 1.0 / (distance as f64) } else { 0.0 };
  let i_step = (end_cube.i - start_cube.i) as f64 * inv_distance_f;
  let j_step = (end_cube.j - start_cube.j) as f64 * inv_distance_f;
  let k_step = (end_cube.k - start_cube.k) as f64 * inv_distance_f;

  let mut current_cube_f = CoordIJK {
    // Use CoordIJK to store f64 temp for cube coords
    i: start_cube.i as i32, // These are actually f64 conceptually
    j: start_cube.j as i32,
    k: start_cube.k as i32,
  };

  let mut current_ijk_rounded = CoordIJK::default();

  for n in 0..=distance {
    cube_round(
      start_cube.i as f64 + i_step * (n as f64),
      start_cube.j as f64 + j_step * (n as f64),
      start_cube.k as f64 + k_step * (n as f64),
      &mut current_ijk_rounded, // This is now integer cube coords
    );

    // Convert rounded cube coords back to H3 IJK+ (normalized)
    cube_to_ijk(&mut current_ijk_rounded);

    // Convert the IJK+ (which is local to `start`) back to an H3 index.
    local_ijk_to_cell(start, &current_ijk_rounded, &mut out_path[n as usize])?;
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::indexing::lat_lng_to_cell;
  use crate::latlng::_set_geo_degs;
  use crate::traversal::neighbors::are_neighbor_cells;
  use crate::types::{LatLng, H3_NULL}; // For validation

  #[test]
  fn test_grid_path_cells_size_and_path_identity() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let h = lat_lng_to_cell(&geo, 5).unwrap();

    assert_eq!(grid_path_cells_size(h, h), Ok(1));
    let mut path = [H3_NULL; 1];
    assert!(grid_path_cells(h, h, &mut path).is_ok());
    assert_eq!(path[0], h);
  }

  #[test]
  fn test_grid_path_cells_direct_neighbor() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let origin = lat_lng_to_cell(&geo, 5).unwrap();

    let mut k1_ring = [H3_NULL; 7];
    crate::traversal::grid_disk::grid_disk(origin, 1, &mut k1_ring).unwrap();
    let neighbor = k1_ring.iter().find(|&&cell| cell != H3_NULL && cell != origin).unwrap();

    let path_size = grid_path_cells_size(origin, *neighbor).unwrap();
    assert_eq!(path_size, 2, "Path size to direct neighbor is 2");

    let mut path = vec![H3_NULL; path_size as usize];
    assert!(grid_path_cells(origin, *neighbor, &mut path).is_ok());
    assert_eq!(path[0], origin);
    assert_eq!(path[1], *neighbor);
  }

  #[test]
  fn test_grid_path_cells_res_mismatch() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.779, -122.419);
    let h_res5 = lat_lng_to_cell(&geo, 5).unwrap();
    let h_res6 = lat_lng_to_cell(&geo, 6).unwrap();
    assert_eq!(grid_path_cells_size(h_res5, h_res6), Err(H3Error::ResMismatch));
    let mut path = [H3_NULL; 1]; // Dummy buffer
    assert_eq!(grid_path_cells(h_res5, h_res6, &mut path), Err(H3Error::ResMismatch));
  }

  #[test]
  fn test_grid_path_cells_properties_longer_path() {
    let start = lat_lng_to_cell(
      &LatLng {
        lat: (20.0_f64).to_radians(),
        lng: (10.0_f64).to_radians(),
      },
      5,
    )
    .unwrap();
    let end = lat_lng_to_cell(
      &LatLng {
        lat: (20.0_f64).to_radians(),
        lng: (10.5_f64).to_radians(),
      },
      5,
    )
    .unwrap();

    let path_size_res = grid_path_cells_size(start, end);
    assert!(path_size_res.is_ok(), "Path size calculation failed unexpectedly");
    let path_size = path_size_res.unwrap();
    assert!(path_size > 2, "Test path should be longer than direct neighbor");

    let mut path = vec![H3_NULL; path_size as usize];
    assert!(grid_path_cells(start, end, &mut path).is_ok());

    assert_eq!(path[0], start, "Path starts with start index");
    assert_eq!(path[(path_size - 1) as usize], end, "Path ends with end index");

    for i in 1..(path_size as usize) {
      assert!(is_valid_cell(path[i]), "Index in path is valid: {:x}", path[i].0);
      assert!(
        are_neighbor_cells(path[i], path[i - 1]).unwrap_or(false),
        "Index {:x} is neighbor of previous {:x}",
        path[i].0,
        path[i - 1].0
      );
      if i > 1 {
        assert!(
          !are_neighbor_cells(path[i], path[i - 2]).unwrap_or(true), // unwrap_or(true) to fail if are_neighbor_cells errors
          "Index {:x} is not neighbor of {:x} (should be on a line)",
          path[i].0,
          path[i - 2].0
        );
      }
    }
  }

  // Add tests for pentagon involvement and very long distances that might fail
  // cellToLocalIjk or localIjkToCell.
}
