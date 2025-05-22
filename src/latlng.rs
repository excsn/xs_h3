// src/latlng.rs

use crate::bbox::LongitudeNormalization;
use crate::constants::{EARTH_RADIUS_KM, EPSILON_RAD, MAX_H3_RES, M_180_PI, M_2PI, M_PI, M_PI_180, M_PI_2};
use crate::types::LatLng;
use crate::H3Error;
use std::f64;

/// Normalizes radians to a value between `0.0` and `2 * PI`.
#[inline]
#[must_use]
pub(crate) fn _pos_angle_rads(mut rads: f64) -> f64 {
  // First, shift up until non-negative
  while rads < 0.0 {
    rads += M_2PI;
  }
  // Then shift down until below 2π
  while rads >= M_2PI {
    rads -= M_2PI;
  }
  // Collapse negative zero to positive zero
  if rads == -0.0 {
    0.0
  } else {
    rads
  }
}

/// Determines if the components of two spherical coordinates are within some
/// threshold distance of each other.
#[inline]
#[must_use]
pub(crate) fn geo_almost_equal_threshold(p1: &LatLng, p2: &LatLng, threshold: f64) -> bool {
  (p1.lat - p2.lat).abs() < threshold && (p1.lng - p2.lng).abs() < threshold
}

/// Determines if the components of two spherical coordinates are within the
/// standard H3 epsilon distance of each other.
#[inline]
#[must_use]
pub(crate) fn geo_almost_equal(p1: &LatLng, p2: &LatLng) -> bool {
  geo_almost_equal_threshold(p1, p2, EPSILON_RAD)
}

/// Set the components of spherical coordinates in decimal degrees.
#[inline]
pub(crate) fn _set_geo_degs(p: &mut LatLng, lat_degs: f64, lng_degs: f64) {
  p.lat = lat_degs.to_radians();
  p.lng = lng_degs.to_radians();
}

/// Set the components of spherical coordinates in radians.
#[inline]
pub(crate) fn _set_geo_rads(p: &mut LatLng, lat_rads: f64, lng_rads: f64) {
  p.lat = lat_rads;
  p.lng = lng_rads;
}

/// Constrains latitude to the range `[-PI/2, PI/2]` by folding.
#[inline]
#[must_use]
pub(crate) fn constrain_lat(lat: f64) -> f64 {
  let mut x = (lat + M_PI).rem_euclid(M_2PI) - M_PI; // Normalize to [-PI, PI)
  if x > M_PI_2 {
    x = M_PI - x; // Fold from > PI/2 down
  } else if x < -M_PI_2 {
    x = -M_PI - x; // Fold from < -PI/2 up
  }
  x
}

/// Constrains longitude to the range `[-PI, PI)`.
#[inline]
#[must_use]
pub(crate) fn constrain_lng(mut lng: f64) -> f64 {
  while lng >= M_PI { // Changed > to >= to match common [-PI, PI) interval
    lng -= M_2PI;
  }
  while lng < -M_PI {
    lng += M_2PI;
  }
  lng
}

/// Determines the azimuth from p1 to p2 in radians.
/// Output is normalized to [0, 2PI).
#[inline]
#[must_use]
pub(crate) fn _geo_azimuth_rads(p1: &LatLng, p2: &LatLng) -> f64 {
    let y = p2.lat.cos() * (p2.lng - p1.lng).sin();
    let x = p1.lat.cos() * p2.lat.sin() - 
            p1.lat.sin() * p2.lat.cos() * (p2.lng - p1.lng).cos();
    _pos_angle_rads(y.atan2(x))
}

/// Computes the point on the sphere a specified azimuth and distance from another point.
#[inline]
pub(crate) fn _geo_az_distance_rads(p1: &LatLng, az: f64, distance: f64, p2: &mut LatLng) {
    const LOCAL_EPSILON: f64 = EPSILON_RAD;

    if distance < LOCAL_EPSILON {
        *p2 = *p1;
        return;
    }

    let az_norm = _pos_angle_rads(az);

    // Check for due north/south azimuth
    if az_norm < LOCAL_EPSILON || (az_norm - M_PI).abs() < LOCAL_EPSILON {
        let raw_lat = if az_norm < LOCAL_EPSILON { // due north
            p1.lat + distance
        } else { // due south
            p1.lat - distance
        };
        
        // p2.lat is not constrained by constrain_lat here to match C _geoAzDistanceRads behavior,
        // which can produce latitudes outside [-PI/2, PI/2] that are handled by later logic or tests.
        p2.lat = raw_lat;

        if (p2.lat - M_PI_2).abs() < LOCAL_EPSILON { // north pole
            p2.lat = M_PI_2;
            p2.lng = 0.0;
        } else if (p2.lat + M_PI_2).abs() < LOCAL_EPSILON { // south pole
            p2.lat = -M_PI_2;
            p2.lng = 0.0;
        } else {
            // C H3's _geoAzDistanceRads sets p2.lng = constrainLng(p1->lng) in this branch.
            // This doesn't account for longitude flipping if a pole was crossed.
            // However, C tests like _geoAzDistanceRads_dueNorthSouth ("due north to south pole")
            // expect an out-of-range latitude (e.g., 270 deg) and the *original* longitude.
            // `geoAlmostEqual` must handle these out-of-range lats for tests to pass, or tests use constrainLat.
            // For a direct port of this function as used by H3 internals:
            p2.lng = constrain_lng(p1.lng);
        }
    } else { // General Azimuth
        let mut sinlat = p1.lat.sin() * distance.cos() +
                       p1.lat.cos() * distance.sin() * az_norm.cos();
        
        if sinlat > 1.0 { sinlat = 1.0; }
        if sinlat < -1.0 { sinlat = -1.0; }
        p2.lat = sinlat.asin();

        if (p2.lat - M_PI_2).abs() < LOCAL_EPSILON { // Landed exactly on North pole
            p2.lat = M_PI_2;
            p2.lng = 0.0;
        } else if (p2.lat + M_PI_2).abs() < LOCAL_EPSILON { // Landed exactly on South pole
            p2.lat = -M_PI_2;
            p2.lng = 0.0;
        } else {
            let cos_p1_lat = p1.lat.cos();
            let invcosp2lat = 1.0 / p2.lat.cos(); // Safe, p2 is not a pole here

            let mut sinlng = az_norm.sin() * distance.sin() * invcosp2lat;
            
            let coslng_numerator = distance.cos() - p1.lat.sin() * p2.lat.sin();
            let coslng_denominator = cos_p1_lat * p2.lat.cos(); 
            
            let mut coslng;
            if coslng_denominator.abs() < LOCAL_EPSILON {
                // This implies p1 is a pole (since p2 is not a pole).
                // Standard formulas simplify significantly if p1 is a pole.
                // If p1 is a pole, delta_lng is related to az_norm directly.
                // C code's direct division would lead to Inf/NaN if not for clamping.
                // A robust approach might special-case p1 pole here:
                // if p1.lat.cos().abs() < LOCAL_EPSILON { p2.lng = constrain_lng(az_norm); return; }
                // However, to match C's _geoAzDistanceRads path:
                coslng = if coslng_numerator.abs() < LOCAL_EPSILON { 1.0 } // 0/0 type case
                         else if coslng_numerator > 0.0 { f64::INFINITY } 
                         else { f64::NEG_INFINITY };
            } else {
                coslng = coslng_numerator / coslng_denominator;
            }

            if sinlng > 1.0 { sinlng = 1.0; }
            if sinlng < -1.0 { sinlng = -1.0; }
            if coslng > 1.0 { coslng = 1.0; }
            if coslng < -1.0 { coslng = -1.0; }
            
            p2.lng = constrain_lng(p1.lng + sinlng.atan2(coslng));
        }
    }
}


// --- Public API Geo Functions (mirrors C H3 API) ---

/// The great circle distance in radians between two spherical coordinates.
pub fn great_circle_distance_rads(a: &LatLng, b: &LatLng) -> f64 {
  let sin_lat_half = ((b.lat - a.lat) * 0.5).sin();
  let sin_lng_half = ((b.lng - a.lng) * 0.5).sin();
  let component_a = sin_lat_half * sin_lat_half + a.lat.cos() * b.lat.cos() * sin_lng_half * sin_lng_half;
  let component_a_clamped = component_a.max(0.0).min(1.0);
  2.0 * component_a_clamped.sqrt().atan2((1.0 - component_a_clamped).sqrt())
}

/// The great circle distance in kilometers between two spherical coordinates.
pub fn great_circle_distance_km(a: &LatLng, b: &LatLng) -> f64 {
  great_circle_distance_rads(a, b) * EARTH_RADIUS_KM
}

/// The great circle distance in meters between two spherical coordinates.
pub fn great_circle_distance_m(a: &LatLng, b: &LatLng) -> f64 {
  great_circle_distance_km(a, b) * 1000.0
}

/// Converts degrees to radians.
pub fn degs_to_rads(degrees: f64) -> f64 {
  degrees * M_PI_180
}

/// Converts radians to degrees.
pub fn rads_to_degs(radians: f64) -> f64 {
  radians * M_180_PI
}

/// Average hexagon area in square kilometers (excludes pentagons) at the given resolution.
pub fn get_hexagon_area_avg_km2(res: i32) -> Result<f64, H3Error> {
  const AREAS_KM2: [f64; (MAX_H3_RES + 1) as usize] = [
    4.357_449_416_078_383e+06, 6.097_884_417_941_332e+05, 8.680_178_039_899_720e+04,
    1.239_343_465_508_816e+04, 1.770_347_654_491_307e+03, 2.529_038_581_819_449e+02,
    3.612_906_216_441_245e+01, 5.161_293_359_717_191e+00, 7.373_275_975_944_177e-01,
    1.053_325_134_272_067e-01, 1.504_750_190_766_435e-02, 2.149_643_129_451_879e-03,
    3.070_918_756_316_060e-04, 4.387_026_794_728_296e-05, 6.267_181_135_324_313e-06,
    8.953_115_907_605_790e-07,
  ];
  if res < 0 || res > MAX_H3_RES { return Err(H3Error::ResDomain); }
  Ok(AREAS_KM2[res as usize])
}

/// Average hexagon area in square meters (excludes pentagons) at the given resolution.
pub fn get_hexagon_area_avg_m2(res: i32) -> Result<f64, H3Error> {
  const AREAS_M2: [f64; (MAX_H3_RES + 1) as usize] = [
    4.357_449_416_078_390e+12, 6.097_884_417_941_339e+11, 8.680_178_039_899_731e+10,
    1.239_343_465_508_818e+10, 1.770_347_654_491_309e+09, 2.529_038_581_819_452e+08,
    3.612_906_216_441_250e+07, 5.161_293_359_717_198e+06, 7.373_275_975_944_188e+05,
    1.053_325_134_272_069e+05, 1.504_750_190_766_437e+04, 2.149_643_129_451_882e+03,
    3.070_918_756_316_063e+02, 4.387_026_794_728_301e+01, 6.267_181_135_324_322,
    8.953_115_907_605_802e-01,
  ];
  if res < 0 || res > MAX_H3_RES { return Err(H3Error::ResDomain); }
  Ok(AREAS_M2[res as usize])
}

/// Average hexagon edge length in kilometers (excludes pentagons) at the given resolution.
pub fn get_hexagon_edge_length_avg_km(res: i32) -> Result<f64, H3Error> {
  const LENS_KM: [f64; (MAX_H3_RES + 1) as usize] = [
    1281.256011, 483.0568391, 182.5129565, 68.97922179, 26.07175968, 9.854090990,
    3.724532667, 1.406475763, 0.531414010, 0.200786148, 0.075863783, 0.028663897,
    0.010830188, 0.004092010, 0.001546100, 0.000584169,
  ];
  if res < 0 || res > MAX_H3_RES { return Err(H3Error::ResDomain); }
  Ok(LENS_KM[res as usize])
}

/// Average hexagon edge length in meters (excludes pentagons) at the given resolution.
pub fn get_hexagon_edge_length_avg_m(res: i32) -> Result<f64, H3Error> {
  const LENS_M: [f64; (MAX_H3_RES + 1) as usize] = [
    1281256.011, 483056.8391, 182512.9565, 68979.22179, 26071.75968, 9854.090990,
    3724.532667, 1406.475763, 531.4140101, 200.7861476, 75.86378287, 28.66389748,
    10.83018784, 4.092010473, 1.546099657, 0.584168630,
  ];
  if res < 0 || res > MAX_H3_RES { return Err(H3Error::ResDomain); }
  Ok(LENS_M[res as usize])
}

/// Normalizes a longitude based on the specified strategy, typically for comparisons.
#[inline]
#[must_use]
pub(crate) fn normalize_lng_for_comparison(lng: f64, normalization: LongitudeNormalization) -> f64 {
  match normalization {
    LongitudeNormalization::None => lng,
    LongitudeNormalization::East => {
      if lng < 0.0 { lng + M_2PI } else { lng }
    }
    LongitudeNormalization::West => {
      if lng > 0.0 { lng - M_2PI } else { lng }
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::constants::{EPSILON_DEG, M_PI_180}; 
  // use crate::types::LatLng; // Already imported by use super::*;

  #[test]
  fn test_pos_angle_rads() {
    assert!((_pos_angle_rads(0.0) - 0.0).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI) - M_PI).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI * 1.5) - M_PI * 1.5).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI * 2.0) - 0.0).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI * 2.5) - M_PI * 0.5).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI_2) - (M_PI * 1.5)).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI) - M_PI).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI * 1.5) - (M_PI * 0.5)).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI * 2.0) - 0.0).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI * 4.0) - 0.0).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI * 4.0) - 0.0).abs() < f64::EPSILON);
  }

  #[test]
  fn test_geo_almost_equal_threshold() {
    let a = LatLng { lat: 15.0 * M_PI_180, lng: 10.0 * M_PI_180 };
    let mut b = LatLng { lat: 15.0 * M_PI_180, lng: 10.0 * M_PI_180 };
    assert!(geo_almost_equal_threshold(&a, &b, EPSILON_RAD / 2.0));
    b.lat = (15.0 + EPSILON_DEG * 2.0) * M_PI_180;
    b.lng = (10.0 + EPSILON_DEG * 2.0) * M_PI_180;
    assert!(geo_almost_equal_threshold(&a, &b, EPSILON_RAD * 3.0));
    b.lat = (15.0 + EPSILON_DEG * 2.0) * M_PI_180;
    b.lng = 10.0 * M_PI_180;
    assert!(!geo_almost_equal_threshold(&a, &b, EPSILON_RAD));
    b.lat = 15.0 * M_PI_180;
    b.lng = (10.0 + EPSILON_DEG * 2.0) * M_PI_180;
    assert!(!geo_almost_equal_threshold(&a, &b, EPSILON_RAD));
  }

  #[test]
  fn test_constrain_lat() {
    assert!((constrain_lat(0.0) - 0.0).abs() < f64::EPSILON);
    assert!((constrain_lat(1.0) - 1.0).abs() < f64::EPSILON);
    assert!((constrain_lat(M_PI_2) - M_PI_2).abs() < f64::EPSILON);
    assert!((constrain_lat(M_PI) - 0.0).abs() < f64::EPSILON); // Folds to 0
    assert!((constrain_lat(M_PI + 1.0) - (M_PI - (M_PI + 1.0))).abs() < f64::EPSILON); // Folds: PI - (lat - PI) = 2PI - lat
                                                                                    // (M_PI + 1.0) becomes (M_PI + 1.0 + M_PI) % M_2PI - M_PI = (1.0) % M_2PI - M_PI = 1.0 - M_PI. (Incorrect expectation)
                                                                                    // Correct: M_PI + 1 => (1.0 + 2M_PI) % M_2PI - M_PI = 1.0 - M_PI
                                                                                    // if 1.0-M_PI < -M_PI_2, then -M_PI - (1.0-M_PI) = -1.0
                                                                                    // if 1.0-M_PI > M_PI_2, then M_PI - (1.0-M_PI) = 2M_PI - 1.0
                                                                                    // Lat=PI+1 (4.14). Norm to [-PI,PI) -> 4.14-2PI = -2.14. Fold: -PI - (-2.14) = -PI+2.14 = -1.0
    assert!((constrain_lat(M_PI + 1.0) - -1.0).abs() < f64::EPSILON);


    assert!((constrain_lat(M_2PI + 1.0) - 1.0).abs() < f64::EPSILON);
    assert!((constrain_lat(-M_PI_2) - -M_PI_2).abs() < f64::EPSILON);
    assert!((constrain_lat(-M_PI) - 0.0).abs() < f64::EPSILON);
    assert!((constrain_lat(-M_PI - 1.0) - 1.0).abs() < f64::EPSILON); // -PI-1 (-4.14) -> norm: -4.14+2PI = 2.14. Fold: M_PI - 2.14 = 1.0
  }

  #[test]
  fn test_constrain_lng() {
    assert!((constrain_lng(0.0) - 0.0).abs() < f64::EPSILON);
    assert!((constrain_lng(1.0) - 1.0).abs() < f64::EPSILON);
    // constrain_lng aims for [-PI, PI). M_PI becomes -M_PI.
    assert!((constrain_lng(M_PI) - (-M_PI)).abs() < f64::EPSILON); // wraps to -PI if input is PI
    assert!((constrain_lng(M_2PI) - 0.0).abs() < f64::EPSILON);
    assert!((constrain_lng(M_PI * 3.0) - (-M_PI)).abs() < f64::EPSILON);
    assert!((constrain_lng(M_PI * 4.0) - 0.0).abs() < f64::EPSILON);
    assert!((constrain_lng(-M_PI) - (-M_PI)).abs() < f64::EPSILON);
    assert!((constrain_lng(-M_2PI) - 0.0).abs() < f64::EPSILON);
  }
  
  #[test]
  fn test_geo_azimuth_rads() {
      let p1 = LatLng { lat: 0.0, lng: 0.0 };
      let p2_north = LatLng { lat: M_PI_2, lng: 0.0 };
      let p2_east = LatLng { lat: 0.0, lng: M_PI_2 };
      let p2_south = LatLng { lat: -M_PI_2, lng: 0.0 };
      let p2_west = LatLng { lat: 0.0, lng: -M_PI_2 };

      assert!((_geo_azimuth_rads(&p1, &p2_north) - 0.0).abs() < EPSILON_RAD, "Azimuth to North"); // Due North
      assert!((_geo_azimuth_rads(&p1, &p2_east) - M_PI_2).abs() < EPSILON_RAD, "Azimuth to East"); // Due East
      assert!((_geo_azimuth_rads(&p1, &p2_south) - M_PI).abs() < EPSILON_RAD, "Azimuth to South"); // Due South
      assert!((_geo_azimuth_rads(&p1, &p2_west) - (3.0 * M_PI_2)).abs() < EPSILON_RAD, "Azimuth to West"); // Due West (positive angle)
  }

  #[test]
  fn test_geo_az_distance_rads_noop() {
    let start = LatLng { lat: (15.0_f64).to_radians(), lng: (10.0_f64).to_radians() };
    let mut out = LatLng::default();
    _geo_az_distance_rads(&start, 0.0, 0.0, &mut out);
    assert!(geo_almost_equal(&start, &out));
  }

  // Using C H3 test values for _geoAzDistanceRads_invertible
  #[test]
  fn test_geo_az_distance_rads_invertible() {
      let mut start = LatLng::default();
      _set_geo_degs(&mut start, 15.0, 10.0);
      let mut out = LatLng::default();

      let azimuth = degs_to_rads(20.0);
      let degrees180 = degs_to_rads(180.0);
      let distance = degs_to_rads(15.0);

      _geo_az_distance_rads(&start, azimuth, distance, &mut out);
      assert!( (great_circle_distance_rads(&start, &out) - distance).abs() < EPSILON_RAD,
               "moved distance is as expected. Dist: {}, Expected: {}", great_circle_distance_rads(&start, &out), distance);

      let start2 = out; // Use the point we moved to
      _geo_az_distance_rads(&start2, azimuth + degrees180, distance, &mut out);
      // Epsilon for this check in C is < 0.01 radians, which is quite large.
      // This suggests the inverse operation isn't perfectly precise due to formula choices or float issues.
      assert!(great_circle_distance_rads(&start, &out) < 0.01,
               "moved back to origin (within 0.01 rad). Actual dist: {}", great_circle_distance_rads(&start, &out));
  }
  
  // C test `_geoAzDistanceRads_dueNorthSouth` "due north to south pole" expected:
  // start: (45,1), az:0, dist: (45+180)rads => out_lat = 45+45+180 = 270 deg. out_lng = 1 deg.
  // This confirms C's _geoAzDistanceRads can produce lat outside [-PI/2, PI/2].
  #[test]
  fn test_geo_az_distance_rads_c_north_to_south_pole_case() {
    let mut start = LatLng::default();
    _set_geo_degs(&mut start, 45.0, 1.0);
    let mut out = LatLng::default();
    let distance = degs_to_rads(45.0 + 180.0); // To cross North Pole, pass South Pole, and end up on other side
    
    _geo_az_distance_rads(&start, 0.0, distance, &mut out);

    // Expected from C test (lat=270 deg, lng=1 deg)
    // Lat 270 deg = 3*PI/2 rads.
    let expected_lat_raw = degs_to_rads(270.0); // This is 4.712...
    let expected_lng_raw = degs_to_rads(1.0);

    assert!((out.lat - expected_lat_raw).abs() < EPSILON_RAD, "Latitude for due north past south pole. Got {}, Exp {}", out.lat, expected_lat_raw);
    assert!((out.lng - expected_lng_raw).abs() < EPSILON_RAD, "Longitude for due north past south pole. Got {}, Exp {}", out.lng, expected_lng_raw);
  }


}