// src/math/vec3d.rs

use crate::types::{LatLng, Vec3d};
use std::f64;

/// Square of a number.
#[inline]
fn _square(x: f64) -> f64 {
  x * x
}

/// Calculate the square of the Euclidean distance between two 3D coordinates.
#[inline]
#[must_use]
pub(crate) fn _point_square_dist(v1: &Vec3d, v2: &Vec3d) -> f64 {
  _square(v1.x - v2.x) + _square(v1.y - v2.y) + _square(v1.z - v2.z)
}

/// Calculate the 3D Cartesian coordinate on a unit sphere from latitude and longitude.
///
/// # Arguments
///
/// * `geo` - The latitude and longitude of the point (in radians).
/// * `point` - Output: The 3D Cartesian coordinate.
#[inline]
pub(crate) fn _geo_to_vec3d(geo: &LatLng, point: &mut Vec3d) {
  let r = geo.lat.cos();

  point.z = geo.lat.sin();
  point.x = geo.lng.cos() * r;
  point.y = geo.lng.sin() * r;
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::constants::{EPSILON_RAD, M_PI_2}; // For M_PI, M_PI_2
  use crate::types::{LatLng, Vec3d};

  #[test]
  fn test_point_square_dist() {
    let v1 = Vec3d { x: 0.0, y: 0.0, z: 0.0 };
    let v2 = Vec3d { x: 1.0, y: 0.0, z: 0.0 };
    let v3 = Vec3d { x: 0.0, y: 1.0, z: 1.0 };
    let v4 = Vec3d { x: 1.0, y: 1.0, z: 1.0 };
    let v5 = Vec3d { x: 1.0, y: 1.0, z: 2.0 };

    assert!(
      (_point_square_dist(&v1, &v1)).abs() < f64::EPSILON,
      "distance to self is 0"
    );
    assert!(
      (_point_square_dist(&v1, &v2) - 1.0).abs() < f64::EPSILON,
      "distance to <1,0,0> is 1"
    );
    assert!(
      (_point_square_dist(&v1, &v3) - 2.0).abs() < f64::EPSILON,
      "distance to <0,1,1> is 2"
    );
    assert!(
      (_point_square_dist(&v1, &v4) - 3.0).abs() < f64::EPSILON,
      "distance to <1,1,1> is 3"
    );
    assert!(
      (_point_square_dist(&v1, &v5) - 6.0).abs() < f64::EPSILON,
      "distance to <1,1,2> is 6"
    );
  }

  #[test]
  fn test_geo_to_vec3d() {
    let origin = Vec3d::default();
    let mut p1 = Vec3d::default();
    let mut p2 = Vec3d::default();
    let mut p3 = Vec3d::default();

    let c1 = LatLng { lat: 0.0, lng: 0.0 };
    _geo_to_vec3d(&c1, &mut p1);
    // Check that point is on unit sphere (squared distance from origin is 1)
    assert!(
      (_point_square_dist(&origin, &p1) - 1.0).abs() < EPSILON_RAD,
      "Geo point (0,0) is on the unit sphere"
    );
    // Expected: x=1, y=0, z=0 for (lat=0, lng=0)
    assert!((p1.x - 1.0).abs() < f64::EPSILON);
    assert!((p1.y - 0.0).abs() < f64::EPSILON);
    assert!((p1.z - 0.0).abs() < f64::EPSILON);

    let c2 = LatLng { lat: M_PI_2, lng: 0.0 }; // North Pole
    _geo_to_vec3d(&c2, &mut p2);
    assert!(
      (_point_square_dist(&origin, &p2) - 1.0).abs() < EPSILON_RAD,
      "Geo point (North Pole) is on the unit sphere"
    );
    // Expected: x=0, y=0, z=1 for North Pole
    assert!((p2.x - 0.0).abs() < f64::EPSILON);
    assert!((p2.y - 0.0).abs() < f64::EPSILON);
    assert!((p2.z - 1.0).abs() < f64::EPSILON);

    // Check distance between (0,0) and North Pole (pi/2,0)
    // Squared distance should be (1-0)^2 + (0-0)^2 + (0-1)^2 = 1+0+1 = 2
    assert!(
      (_point_square_dist(&p1, &p2) - 2.0).abs() < EPSILON_RAD,
      "Squared distance between (0,0) and North Pole is 2"
    );

    let c3 = LatLng {
      lat: crate::constants::M_PI,
      lng: 0.0,
    }; // This is not a valid lat, but tests math
       // M_PI radians = 180 degrees lat. cos(pi) = -1, sin(pi) = 0.
       // This effectively mirrors the (0,0) point through the Earth's center
       // resulting in x = -1, y=0, z=0 (if lat was allowed to be PI)
       // However, LatLng implies normalized coordinates usually.
       // The C code doesn't normalize before this, so we test its math.
       // cos(PI) = -1. So r = -1.
       // z = sin(PI) = 0.
       // x = cos(0) * -1 = -1.
       // y = sin(0) * -1 = 0.
       // So Vec3d should be (-1, 0, 0).
    _geo_to_vec3d(&c3, &mut p3);
    assert!(
      (_point_square_dist(&origin, &p3) - 1.0).abs() < EPSILON_RAD,
      "Geo point (PI,0) is on the unit sphere (mathematically)"
    );
    // Expected: x=-1, y=0, z=0
    assert!((p3.x - (-1.0)).abs() < f64::EPSILON);
    assert!((p3.y - 0.0).abs() < f64::EPSILON);
    assert!((p3.z - 0.0).abs() < f64::EPSILON);

    // Squared distance between (0,0) [p1: (1,0,0)] and (PI,0) [p3: (-1,0,0)] should be (1 - (-1))^2 = 2^2 = 4
    assert!(
      (_point_square_dist(&p1, &p3) - 4.0).abs() < EPSILON_RAD,
      "Squared distance between (0,0) and (PI,0) [antipodal on equator] is 4"
    );
  }
}
