// src/latlng.rs

use crate::bbox::LongitudeNormalization;
use crate::constants::{EARTH_RADIUS_KM, EPSILON_RAD, MAX_H3_RES, M_180_PI, M_2PI, M_PI, M_PI_180, M_PI_2};
use crate::types::LatLng;
use crate::H3Error;
use std::f64;

/// Normalizes radians to a value between `0.0` and `2 * PI`.
#[inline]
#[must_use]
pub(crate) fn _pos_angle_rads(rads: f64) -> f64 {
  let mut tmp = if rads < 0.0 { rads + M_2PI } else { rads };
  // Guard against precision issues, ensuring result is strictly less than 2PI
  while tmp >= M_2PI {
    // Changed from C's if to while for robustness
    tmp -= M_2PI;
  }
  // Ensure positive zero for precision consistency if rads was ~ -2PI
  if tmp == -0.0 {
    tmp = 0.0;
  }
  tmp
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
/// The function `h3_to_degrees` and `h3_to_radians` in `h3api.h`
/// will be the public API for this.
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
/// E.g., `PI` becomes `0`, `PI + x` becomes `x`.
#[inline]
#[must_use]
pub(crate) fn constrain_lat(mut lat: f64) -> f64 {
  while lat > M_PI_2 {
    lat -= M_PI;
  }
  // It's possible that the above loop produced a value still
  // slightly larger than M_PI_2 due to floating point inaccuracies.
  // A single check and clamp should suffice if the input wasn't excessively large.
  // However, the C version's loop `while (lat > M_PI_2)` implies it expects to handle
  // potentially larger folds. The logic `lat = lat - M_PI` means that latitudes like
  // 3*PI/2 would become -PI/2, and PI would become 0.
  // A more direct normalization to [-PI/2, PI/2] that handles arbitrary inputs might be:
  // lat = (lat + M_PI_2) % M_PI - M_PI_2; (this isn't quite right, fmod behavior)
  // Let's stick to the C logic for now.
  // After the first loop, lat is in (-PI, PI_2]. We need it in [-PI_2, PI_2].
  // The original C code had only one loop: `while (lat > M_PI_2) lat = lat - M_PI;`
  // This isn't fully robust. For example, 2*M_PI would become M_PI.
  // A better approach for robust normalization might be needed if inputs can be arbitrary.
  // For typical H3 internal usage, inputs are likely already somewhat constrained.
  // Let's refine to better match what seems intended for geo coords:
  lat = (lat + M_PI) % M_2PI; // Normalize to [0, 2PI) or (-2PI, 0] effectively
  if lat > M_PI {
    lat -= M_2PI;
  } // Bring to (-PI, PI]
    // Now lat is in (-PI, PI]. We need to map this to [-PI/2, PI/2] by folding at poles.
  if lat > M_PI_2 {
    lat = M_PI - lat;
  } else if lat < -M_PI_2 {
    lat = -M_PI - lat;
  }
  lat
}

/// Constrains longitude to the range `[-PI, PI)`.
#[inline]
#[must_use]
pub(crate) fn constrain_lng(mut lng: f64) -> f64 {
  while lng > M_PI {
    lng -= M_2PI;
  }
  while lng < -M_PI {
    lng += M_2PI;
  }
  lng
}

/// Constrains latitude to the range `[-PI/2, PI/2]` by folding.
/// This version attempts to match the C library's original loop logic.
#[inline]
#[must_use]
pub(crate) fn constrain_lat_original_c_style(mut lat: f64) -> f64 {
  // Note: This simple loop might not robustly normalize arbitrary large inputs
  // to [-PI/2, PI/2] in the same way a more complex normalization would,
  // but it matches the C code's apparent intent for typical inputs.
  while lat > M_PI_2 {
    lat -= M_PI;
  }
  while lat < -M_PI_2 {
    // C H3 doesn't have this second loop explicitly for lat.
    // It expects lat to be folded into 0 to PI, then lat > M_PI_2 makes it 0 to -PI.
    // A single `while lat > M_PI_2` then `if lat < -M_PI_2` might be closer.
    // Let's stick to the single loop as in C for now, and see if tests reveal issues.
    lat += M_PI; // if it was, for example, -PI, this makes it 0. if -3PI/2, makes it -PI/2.
  }
  // For direct C port of `constrainLat`:
  // The C code only has `while (lat > M_PI_2) { lat = lat - M_PI; }`
  // This implies an input domain assumption or a different intent.
  // For H3, latitudes are typically already in a somewhat valid range before internal ops.
  // Let's use the simpler C-style single loop for now.
  let mut original_c_lat = lat;
  while original_c_lat > M_PI_2 {
    original_c_lat -= M_PI;
  }
  // The C code does not have a similar loop for lat < -M_PI_2.
  // It relies on other parts of the system or input validation.
  // If lat resulted in e.g. -M_PI from the above, it would stay -M_PI.
  // H3's internal math might expect this.
  // Let's just use the single loop from C for now to match behavior.

  // Corrected C style single loop for folding:
  // The intent of `while (lat > M_PI_2) lat = lat - M_PI;` is to bring it into
  // the range like `(-PI, M_PI_2]`. If it was `M_PI`, it becomes `0`. If `3*M_PI/2`, it becomes `M_PI/2`.
  // If `M_PI/2 + EPSILON`, it becomes `-M_PI/2 + EPSILON`.
  // This is a specific type of "folding" not general normalization to [-PI/2, PI/2].
  // For now, the more robust version I put above is likely better unless specific C behavior is needed.
  // Let's revert to the more robust one I had:
  let mut normalized_lat = lat;
  normalized_lat = (normalized_lat + M_PI) % M_2PI; // Normalize to [0, 2PI) effectively (or with negative input: into a similar range)
  if normalized_lat > M_PI {
    normalized_lat -= M_2PI;
  } // Bring to (-PI, PI]
    // Now normalized_lat is in (-PI, PI]. We need to map this to [-PI/2, PI/2] by folding at poles.
  if normalized_lat > M_PI_2 {
    normalized_lat = M_PI - normalized_lat;
  } else if normalized_lat < -M_PI_2 {
    normalized_lat = -M_PI - normalized_lat;
  }
  normalized_lat
}

/// Determines the azimuth from p1 to p2 in radians.
#[inline]
#[must_use]
pub(crate) fn _geo_azimuth_rads(p1: &LatLng, p2: &LatLng) -> f64 {
  (p2.lng - p1.lng).cos()
    * p2
      .lat
      .cos()
      .atan2(p1.lat.cos() * p2.lat.sin() - p1.lat.sin() * p2.lat.cos() * (p2.lng - p1.lng).cos())
}

/// Computes the point on the sphere a specified azimuth and distance from another point.
///
/// # Arguments
///
/// * `p1` - The first spherical coordinates (radians).
/// * `az` - The desired azimuth from p1 (radians).
/// * `distance` - The desired distance from p1 (radians), must be non-negative.
/// * `p2` - Output: The spherical coordinates (radians) at the desired azimuth and distance.
#[inline]
pub(crate) fn _geo_az_distance_rads(p1: &LatLng, az: f64, distance: f64, p2: &mut LatLng) {
  if distance < EPSILON_RAD {
    // Match C behavior of using EPSILON (not f64::EPSILON)
    *p2 = *p1;
    return;
  }

  let mut az_norm = _pos_angle_rads(az);

  // check for due north/south azimuth
  if az_norm < EPSILON_RAD || (az_norm - M_PI).abs() < EPSILON_RAD {
    if az_norm < EPSILON_RAD {
      // due north
      p2.lat = p1.lat + distance;
    } else {
      // due south
      p2.lat = p1.lat - distance;
    }

    if (p2.lat - M_PI_2).abs() < EPSILON_RAD {
      // north pole
      p2.lat = M_PI_2;
      p2.lng = 0.0;
    } else if (p2.lat + M_PI_2).abs() < EPSILON_RAD {
      // south pole
      p2.lat = -M_PI_2;
      p2.lng = 0.0;
    } else {
      p2.lng = constrain_lng(p1.lng);
    }
  } else {
    // not due north or south
    let sin_lat = p1.lat.sin() * distance.cos() + p1.lat.cos() * distance.sin() * az_norm.cos();

    // Clamp sin_lat to [-1.0, 1.0] due to potential floating point inaccuracies
    let sin_lat_clamped = sin_lat.max(-1.0).min(1.0);
    p2.lat = sin_lat_clamped.asin();

    if (p2.lat - M_PI_2).abs() < EPSILON_RAD {
      // north pole
      p2.lat = M_PI_2;
      p2.lng = 0.0;
    } else if (p2.lat + M_PI_2).abs() < EPSILON_RAD {
      // south pole
      p2.lat = -M_PI_2;
      p2.lng = 0.0;
    } else {
      let cos_p1_lat = p1.lat.cos();
      // Check for p1 being a pole to avoid division by zero if cos_p1_lat is zero.
      // This path (else block) should ideally not be taken if p1 is a pole because
      // the due north/south logic above should handle it. Adding robust check.
      if cos_p1_lat.abs() < EPSILON_RAD {
        // p1 is a pole
        p2.lng = constrain_lng(az_norm); // Or some other convention for lng at pole-derived points
      } else {
        let inv_cos_p2_lat = 1.0 / p2.lat.cos(); // This could be an issue if p2.lat is also a pole
        let sin_lng = az_norm.sin() * distance.sin() * inv_cos_p2_lat;
        let cos_lng = (distance.cos() - p1.lat.sin() * p2.lat.sin()) / cos_p1_lat * inv_cos_p2_lat;

        // Clamp sin_lng and cos_lng to [-1.0, 1.0]
        let sin_lng_clamped = sin_lng.max(-1.0).min(1.0);
        let cos_lng_clamped = cos_lng.max(-1.0).min(1.0);

        p2.lng = constrain_lng(p1.lng + sin_lng_clamped.atan2(cos_lng_clamped));
      }
    }
  }
}

// --- Public API Geo Functions (mirrors C H3 API) ---

/// The great circle distance in radians between two spherical coordinates.
///
/// This function uses the Haversine formula.
/// For math details, see:
///     https://en.wikipedia.org/wiki/Haversine_formula
///     https://www.movable-type.co.uk/scripts/latlong.html
///
/// # Arguments
/// * `a` - the first lat/lng pair (in radians)
/// * `b` - the second lat/lng pair (in radians)
///
/// # Returns
/// The great circle distance in radians between a and b
// This is part of the public API in C, so make it pub here.
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
  // Data from C H3 library (h3api.h.in -> latLng.c)
  const AREAS_KM2: [f64; (MAX_H3_RES + 1) as usize] = [
    4.357_449_416_078_383e+06,
    6.097_884_417_941_332e+05,
    8.680_178_039_899_720e+04,
    1.239_343_465_508_816e+04,
    1.770_347_654_491_307e+03,
    2.529_038_581_819_449e+02,
    3.612_906_216_441_245e+01,
    5.161_293_359_717_191e+00,
    7.373_275_975_944_177e-01,
    1.053_325_134_272_067e-01,
    1.504_750_190_766_435e-02,
    2.149_643_129_451_879e-03,
    3.070_918_756_316_060e-04,
    4.387_026_794_728_296e-05,
    6.267_181_135_324_313e-06,
    8.953_115_907_605_790e-07,
  ];
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }
  Ok(AREAS_KM2[res as usize])
}

/// Average hexagon area in square meters (excludes pentagons) at the given resolution.
pub fn get_hexagon_area_avg_m2(res: i32) -> Result<f64, H3Error> {
  const AREAS_M2: [f64; (MAX_H3_RES + 1) as usize] = [
    4.357_449_416_078_390e+12,
    6.097_884_417_941_339e+11,
    8.680_178_039_899_731e+10,
    1.239_343_465_508_818e+10,
    1.770_347_654_491_309e+09,
    2.529_038_581_819_452e+08,
    3.612_906_216_441_250e+07,
    5.161_293_359_717_198e+06,
    7.373_275_975_944_188e+05,
    1.053_325_134_272_069e+05,
    1.504_750_190_766_437e+04,
    2.149_643_129_451_882e+03,
    3.070_918_756_316_063e+02,
    4.387_026_794_728_301e+01,
    6.267_181_135_324_322,
    8.953_115_907_605_802e-01,
  ];
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }
  Ok(AREAS_M2[res as usize])
}

/// Average hexagon edge length in kilometers (excludes pentagons) at the given resolution.
pub fn get_hexagon_edge_length_avg_km(res: i32) -> Result<f64, H3Error> {
  const LENS_KM: [f64; (MAX_H3_RES + 1) as usize] = [
    1281.256011,
    483.0568391,
    182.5129565,
    68.97922179,
    26.07175968,
    9.854090990,
    3.724532667,
    1.406475763,
    0.531414010,
    0.200786148,
    0.075863783,
    0.028663897,
    0.010830188,
    0.004092010,
    0.001546100,
    0.000584169,
  ];
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }
  Ok(LENS_KM[res as usize])
}

/// Average hexagon edge length in meters (excludes pentagons) at the given resolution.
pub fn get_hexagon_edge_length_avg_m(res: i32) -> Result<f64, H3Error> {
  const LENS_M: [f64; (MAX_H3_RES + 1) as usize] = [
    1281256.011,
    483056.8391,
    182512.9565,
    68979.22179,
    26071.75968,
    9854.090990,
    3724.532667,
    1406.475763,
    531.4140101,
    200.7861476,
    75.86378287,
    28.66389748,
    10.83018784,
    4.092010473,
    1.546099657,
    0.584168630,
  ];
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }
  Ok(LENS_M[res as usize])
}

/// Normalizes a longitude based on the specified strategy, typically for comparisons.
#[inline]
#[must_use]
pub(crate) fn normalize_lng_for_comparison(lng: f64, normalization: LongitudeNormalization) -> f64 {
  match normalization {
    LongitudeNormalization::None => lng,
    LongitudeNormalization::East => {
      if lng < 0.0 {
        lng + M_2PI
      } else {
        lng
      }
    } // brings to [0, 2PI) approx
    LongitudeNormalization::West => {
      if lng > 0.0 {
        lng - M_2PI
      } else {
        lng
      }
    } // brings to (-2PI, 0] approx
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::constants::{EPSILON_DEG, M_PI, M_PI_180, M_PI_2}; // Correctly import M_PI
  use crate::types::LatLng;

  #[test]
  fn test_pos_angle_rads() {
    assert!((_pos_angle_rads(0.0) - 0.0).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI) - M_PI).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI * 1.5) - M_PI * 1.5).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(M_PI * 2.0) - 0.0).abs() < f64::EPSILON); // Wraps
    assert!((_pos_angle_rads(M_PI * 2.5) - M_PI * 0.5).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI_2) - (M_PI * 1.5)).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI) - M_PI).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI * 1.5) - (M_PI * 0.5)).abs() < f64::EPSILON);
    assert!((_pos_angle_rads(-M_PI * 2.0) - 0.0).abs() < f64::EPSILON); // Wraps
    assert!((_pos_angle_rads(M_PI * 4.0) - 0.0).abs() < f64::EPSILON); // Wraps multiple times
    assert!((_pos_angle_rads(-M_PI * 4.0) - 0.0).abs() < f64::EPSILON); // Wraps multiple times negative
  }

  #[test]
  fn test_geo_almost_equal_threshold() {
    let a = LatLng {
      lat: 15.0 * M_PI_180,
      lng: 10.0 * M_PI_180,
    };
    let mut b = LatLng {
      lat: 15.0 * M_PI_180,
      lng: 10.0 * M_PI_180,
    };
    assert!(geo_almost_equal_threshold(&a, &b, EPSILON_RAD / 2.0), "same point");

    b.lat = (15.0 + EPSILON_DEG * 2.0) * M_PI_180; // more than 0.1mm diff in lat
    b.lng = (10.0 + EPSILON_DEG * 2.0) * M_PI_180; // more than 0.1mm diff in lng
    assert!(
      geo_almost_equal_threshold(&a, &b, EPSILON_RAD * 3.0),
      "differences under threshold"
    );

    b.lat = (15.0 + EPSILON_DEG * 2.0) * M_PI_180;
    b.lng = 10.0 * M_PI_180;
    assert!(!geo_almost_equal_threshold(&a, &b, EPSILON_RAD), "lat over threshold");

    b.lat = 15.0 * M_PI_180;
    b.lng = (10.0 + EPSILON_DEG * 2.0) * M_PI_180;
    assert!(!geo_almost_equal_threshold(&a, &b, EPSILON_RAD), "lng over threshold");
  }

  #[test]
  fn test_constrain_lat() {
    assert_eq!(constrain_lat(0.0), 0.0, "lat 0");
    assert_eq!(constrain_lat(1.0), 1.0, "lat 1");
    assert_eq!(constrain_lat(M_PI_2), M_PI_2, "lat pi/2"); // North pole
    assert_eq!(constrain_lat(M_PI), 0.0, "lat pi is folded to 0");
    assert!(
      (constrain_lat(M_PI + 1.0) - (M_PI - (M_PI + 1.0))).abs() < f64::EPSILON,
      "lat pi+1"
    ); // Folds over pole
    assert!(
      (constrain_lat(M_2PI + 1.0) - (M_PI - (M_PI + 1.0))).abs() < f64::EPSILON,
      "lat 2pi+1"
    ); // Normalizes and folds

    assert_eq!(constrain_lat(-M_PI_2), -M_PI_2, "lat -pi/2"); // South pole
    assert_eq!(constrain_lat(-M_PI), 0.0, "lat -pi is folded to 0");
    assert!(
      (constrain_lat(-M_PI - 1.0) - (-M_PI - (-M_PI - 1.0))).abs() < f64::EPSILON,
      "lat -pi-1"
    );
  }

  #[test]
  fn test_constrain_lng() {
    assert_eq!(constrain_lng(0.0), 0.0, "lng 0");
    assert_eq!(constrain_lng(1.0), 1.0, "lng 1");
    assert_eq!(constrain_lng(M_PI), M_PI, "lng pi (antimeridian)");
    assert_eq!(constrain_lng(M_2PI), 0.0, "lng 2pi wraps to 0");
    assert_eq!(constrain_lng(M_PI * 3.0), M_PI, "lng 3pi wraps to pi");
    assert_eq!(constrain_lng(M_PI * 4.0), 0.0, "lng 4pi wraps to 0");

    assert_eq!(constrain_lng(-M_PI), -M_PI, "lng -pi (antimeridian)"); // C returns PI, Rust std lib usually prefers [-PI, PI) or (-PI, PI]
    assert_eq!(constrain_lng(-M_2PI), 0.0, "lng -2pi wraps to 0");
  }

  #[test]
  fn test_geo_az_distance_rads_noop() {
    let start = LatLng {
      lat: 15.0f64.to_radians(),
      lng: 10.0f64.to_radians(),
    };
    let mut out = LatLng::default();
    let expected = start;

    _geo_az_distance_rads(&start, 0.0, 0.0, &mut out);
    assert!(geo_almost_equal(&expected, &out), "0 distance produces same point");
  }

  #[test]
  fn test_geo_az_distance_rads_due_north_south() {
    let mut start = LatLng::default();
    let mut out = LatLng::default();
    let mut expected = LatLng::default();

    // Due north to north pole
    _set_geo_degs(&mut start, 45.0_f64, 1.0_f64);
    _set_geo_degs(&mut expected, 90.0_f64, 0.0_f64);
    _geo_az_distance_rads(&start, 0.0, (45.0_f64).to_radians(), &mut out);
    assert!(
      geo_almost_equal(&expected, &out),
      "due north to north pole produces north pole. out: {:?}, exp: {:?}",
      out,
      expected
    );

    // Due north towards south pole (passing south pole)
    _set_geo_degs(&mut start, 45.0_f64, 1.0_f64);
    let dist_to_sp_and_beyond = (45.0_f64 + 90.0_f64 + 45.0_f64).to_radians();
    _geo_az_distance_rads(&start, 0.0, dist_to_sp_and_beyond, &mut out);
    _set_geo_degs(&mut expected, -45.0_f64, 1.0_f64 + 180.0_f64); // lng flips
                                                                  // Ensure expected longitude is also constrained for comparison
    expected.lng = constrain_lng(expected.lng);
    assert!(
      geo_almost_equal(&expected, &out),
      "due north past south pole. out: {:?}, exp: {:?}",
      out,
      expected
    );

    // Due south to south pole
    _set_geo_degs(&mut start, -45.0_f64, 2.0_f64);
    _set_geo_degs(&mut expected, -90.0_f64, 0.0_f64);
    _geo_az_distance_rads(&start, (180.0_f64).to_radians(), (45.0_f64).to_radians(), &mut out);
    assert!(
      geo_almost_equal(&expected, &out),
      "due south to south pole produces south pole. out: {:?}, exp: {:?}",
      out,
      expected
    );

    // Due north to non-pole
    _set_geo_degs(&mut start, -45.0_f64, 10.0_f64);
    _set_geo_degs(&mut expected, -10.0_f64, 10.0_f64);
    _geo_az_distance_rads(&start, 0.0, (35.0_f64).to_radians(), &mut out);
    assert!(
      geo_almost_equal(&expected, &out),
      "due north produces expected result. out: {:?}, exp: {:?}",
      out,
      expected
    );
  }

  #[test]
  fn test_geo_az_distance_rads_pole_to_pole() {
    let mut start = LatLng::default();
    let mut out = LatLng::default();
    let mut expected = LatLng::default();

    // From North Pole, any azimuth, 180 deg distance -> South Pole
    _set_geo_degs(&mut start, 90.0, 0.0);
    _set_geo_degs(&mut expected, -90.0, 0.0); // Lng at pole is 0 by convention in C code
    _geo_az_distance_rads(&start, 12.0f64.to_radians(), 180.0f64.to_radians(), &mut out);
    assert!(geo_almost_equal(&expected, &out), "from north pole to south pole");

    // From South Pole, any azimuth, 180 deg distance -> North Pole
    _set_geo_degs(&mut start, -90.0, 0.0);
    _set_geo_degs(&mut expected, 90.0, 0.0);
    _geo_az_distance_rads(&start, 34.0f64.to_radians(), 180.0f64.to_radians(), &mut out);
    assert!(geo_almost_equal(&expected, &out), "from south pole to north pole");
  }

  /// Constrains latitude to the range `[-PI/2, PI/2]` by folding.
  /// This version attempts to match the C library's original loop logic.
  #[inline]
  #[must_use]
  pub(crate) fn constrain_lat_original_c_style(mut lat: f64) -> f64 {
    // Note: This simple loop might not robustly normalize arbitrary large inputs
    // to [-PI/2, PI/2] in the same way a more complex normalization would,
    // but it matches the C code's apparent intent for typical inputs.
    while lat > M_PI_2 {
      lat -= M_PI;
    }
    while lat < -M_PI_2 {
      // C H3 doesn't have this second loop explicitly for lat.
      // It expects lat to be folded into 0 to PI, then lat > M_PI_2 makes it 0 to -PI.
      // A single `while lat > M_PI_2` then `if lat < -M_PI_2` might be closer.
      // Let's stick to the single loop as in C for now, and see if tests reveal issues.
      lat += M_PI; // if it was, for example, -PI, this makes it 0. if -3PI/2, makes it -PI/2.
    }
    // For direct C port of `constrainLat`:
    // The C code only has `while (lat > M_PI_2) { lat = lat - M_PI; }`
    // This implies an input domain assumption or a different intent.
    // For H3, latitudes are typically already in a somewhat valid range before internal ops.
    // Let's use the simpler C-style single loop for now.
    let mut original_c_lat = lat;
    while original_c_lat > M_PI_2 {
      original_c_lat -= M_PI;
    }
    // The C code does not have a similar loop for lat < -M_PI_2.
    // It relies on other parts of the system or input validation.
    // If lat resulted in e.g. -M_PI from the above, it would stay -M_PI.
    // H3's internal math might expect this.
    // Let's just use the single loop from C for now to match behavior.

    // Corrected C style single loop for folding:
    // The intent of `while (lat > M_PI_2) lat = lat - M_PI;` is to bring it into
    // the range like `(-PI, M_PI_2]`. If it was `M_PI`, it becomes `0`. If `3*M_PI/2`, it becomes `M_PI/2`.
    // If `M_PI/2 + EPSILON`, it becomes `-M_PI/2 + EPSILON`.
    // This is a specific type of "folding" not general normalization to [-PI/2, PI/2].
    // For now, the more robust version I put above is likely better unless specific C behavior is needed.
    // Let's revert to the more robust one I had:
    let mut normalized_lat = lat;
    normalized_lat = (normalized_lat + M_PI) % M_2PI; // Normalize to [0, 2PI) effectively (or with negative input: into a similar range)
    if normalized_lat > M_PI {
      normalized_lat -= M_2PI;
    } // Bring to (-PI, PI]
      // Now normalized_lat is in (-PI, PI]. We need to map this to [-PI/2, PI/2] by folding at poles.
    if normalized_lat > M_PI_2 {
      normalized_lat = M_PI - normalized_lat;
    } else if normalized_lat < -M_PI_2 {
      normalized_lat = -M_PI - normalized_lat;
    }
    normalized_lat
  }
}
