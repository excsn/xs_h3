// src/indexing/to_h3.rs

use crate::constants::MAX_H3_RES;
use crate::coords::face_ijk::_geo_to_face_ijk; // pub(crate) in face_ijk.rs
use crate::h3_index::_face_ijk_to_h3; // pub(crate) in h3_index/mod.rs
use crate::types::{H3Error, H3Index, LatLng, H3_NULL};
use std::f64; // For is_finite

/// Finds the H3 index of the cell containing the given `LatLng` point at the specified resolution.
///
/// # Arguments
/// * `geo` - The latitude/longitude point (in radians).
/// * `res` - The H3 resolution (0-15).
///
/// # Returns
/// `Ok(H3Index)` on success, or an `H3Error` if the input is invalid.
pub fn lat_lng_to_cell(geo: &LatLng, res: i32) -> Result<H3Index, H3Error> {
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }
  // Check for invalid lat/lng
  // lat must be in [-PI/2, PI/2] approx, lng in [-PI, PI] approx.
  // isfinite also catches NaN and Infinity.
  if !geo.lat.is_finite()
    || !geo.lng.is_finite()
    || geo.lat.abs() > (crate::constants::M_PI_2 + crate::constants::EPSILON_RAD)
  {
    // Allow slight epsilon for pole checks
    return Err(H3Error::LatLngDomain);
  }
  // Longitude can be outside [-PI, PI] and will be normalized by geoToFaceIjk->geoToHex2d->geoAzimuthRads
  // but extreme non-finite values are bad.

  let mut fijk = crate::types::FaceIJK::default();
  _geo_to_face_ijk(geo, res, &mut fijk);

  let h = _face_ijk_to_h3(&fijk, res);
  if h == H3_NULL {
    // This can happen if _faceIjkToH3 determines the Fijk is out of range
    // for a valid H3 index (e.g., an IJK coordinate that's too large for the resolution
    // on a particular face, which might occur from extreme geo coordinates if not fully
    // caught by the is_finite checks or if _geo_to_face_ijk produces such values).
    Err(H3Error::Failed) // Or a more specific error if one emerges from patterns
  } else {
    Ok(h)
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::h3_index; // To access get_resolution, etc. for validation
  use crate::latlng::_set_geo_degs; // For easy test case setup
  use crate::types::H3_NULL;

  #[test]
  fn test_lat_lng_to_cell_res_domain() {
    let mut geo = LatLng::default();
    _set_geo_degs(&mut geo, 37.77, -122.4);
    assert_eq!(lat_lng_to_cell(&geo, -1), Err(H3Error::ResDomain));
    assert_eq!(lat_lng_to_cell(&geo, 16), Err(H3Error::ResDomain));
  }

  #[test]
  fn test_lat_lng_to_cell_coord_domain() {
    let mut geo_bad_lat = LatLng::default();
    _set_geo_degs(&mut geo_bad_lat, 100.0, -122.4); // Lat out of range
    assert_eq!(lat_lng_to_cell(&geo_bad_lat, 5), Err(H3Error::LatLngDomain));

    let mut geo_nan_lng = LatLng {
      lat: 0.0,
      lng: f64::NAN,
    };
    assert_eq!(lat_lng_to_cell(&geo_nan_lng, 5), Err(H3Error::LatLngDomain));

    let mut geo_inf_lat = LatLng {
      lat: f64::INFINITY,
      lng: 0.0,
    };
    assert_eq!(lat_lng_to_cell(&geo_inf_lat, 5), Err(H3Error::LatLngDomain));
  }

  #[test]
  fn test_lat_lng_to_cell_known_values() {
    // Use known values from C H3 tests if available, or generate some.
    // Example: San Francisco City Hall
    let mut sf_city_hall = LatLng::default();
    _set_geo_degs(&mut sf_city_hall, 37.779265, -122.419277);

    let h_res5 = lat_lng_to_cell(&sf_city_hall, 5).unwrap();
    assert_eq!(h_res5.0, 0x85283083fffffff, "SF City Hall res 5");
    assert_eq!(h3_index::get_resolution(h_res5), 5);

    let h_res10 = lat_lng_to_cell(&sf_city_hall, 10).unwrap();
    assert_eq!(h_res10.0, 0x8a2830828767fff, "SF City Hall res 10");
    assert_eq!(h3_index::get_resolution(h_res10), 10);

    // Test poles
    let mut north_pole = LatLng::default();
    _set_geo_degs(&mut north_pole, 90.0, 0.0);
    let h_pole_res3 = lat_lng_to_cell(&north_pole, 3).unwrap();
    // Expected from C: 830326fffffffff (Base cell 0, then center children)
    assert_eq!(h_pole_res3.0, 0x830326fffffffff, "North Pole res 3");

    let mut south_pole = LatLng::default();
    _set_geo_degs(&mut south_pole, -90.0, 0.0);
    let h_spole_res4 = lat_lng_to_cell(&south_pole, 4).unwrap();
    // Expected from C: 84f2939ffffffff (Base cell 117, then center children)
    assert_eq!(h_spole_res4.0, 0x84f2939ffffffff, "South Pole res 4");
  }
}
