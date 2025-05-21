// src/math/vec2d.rs

use crate::types::Vec2d;
use std::f64;

/// Calculates the magnitude of a 2D Cartesian vector.
#[inline]
#[must_use]
pub(crate) fn _v2d_mag(v: &Vec2d) -> f64 {
  (v.x * v.x + v.y * v.y).sqrt()
}

/// Finds the intersection between two lines.
///
/// Assumes that the lines intersect and that the intersection is not at an
/// endpoint of either line. (This matches the C library's assumptions for this internal function).
///
/// # Arguments
///
/// * `p0` - The first endpoint of the first line.
/// * `p1` - The second endpoint of the first line.
/// * `p2` - The first endpoint of the second line.
/// * `p3` - The second endpoint of the second line.
/// * `inter` - Output: The intersection point.
#[inline]
pub(crate) fn _v2d_intersect(p0: &Vec2d, p1: &Vec2d, p2: &Vec2d, p3: &Vec2d, inter: &mut Vec2d) {
  let s1x = p1.x - p0.x;
  let s1y = p1.y - p0.y;
  let s2x = p3.x - p2.x;
  let s2y = p3.y - p2.y;

  let denominator = -s2x * s1y + s1x * s2y;

  // Note: The C version does not check for parallel lines (denominator == 0).
  // For a robust library, this should be handled, e.g., by returning a Result.
  // For a direct port matching internal C assumptions, we might omit the check.
  // Let's assume the C caller ensures lines are not parallel.
  // if denominator.abs() < f64::EPSILON { /* handle parallel lines */ }

  let t = (s2x * (p0.y - p2.y) - s2y * (p0.x - p2.x)) / denominator;

  inter.x = p0.x + (t * s1x);
  inter.y = p0.y + (t * s1y);
}

/// Checks if two 2D vectors are "almost" equal, within a small epsilon.
/// The C version uses `FLT_EPSILON` which is for `float`. Rust's `f64::EPSILON` is for `double`.
#[inline]
#[must_use]
pub(crate) fn _v2d_almost_equals(v1: &Vec2d, v2: &Vec2d) -> bool {
  (v1.x - v2.x).abs() < f64::EPSILON && (v1.y - v2.y).abs() < f64::EPSILON
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::types::Vec2d; // Ensure Vec2d is in scope

  #[test]
  fn test_v2d_mag() {
    let v = Vec2d { x: 3.0, y: 4.0 };
    let expected = 5.0;
    let mag = _v2d_mag(&v);
    assert!((mag - expected).abs() < f64::EPSILON, "magnitude as expected");
  }

  #[test]
  fn test_v2d_intersect() {
    let p0 = Vec2d { x: 2.0, y: 2.0 };
    let p1 = Vec2d { x: 6.0, y: 6.0 };
    let p2 = Vec2d { x: 0.0, y: 4.0 };
    let p3 = Vec2d { x: 10.0, y: 4.0 };
    let mut intersection = Vec2d::default();

    _v2d_intersect(&p0, &p1, &p2, &p3, &mut intersection);

    let expected_x = 4.0;
    let expected_y = 4.0;

    assert!(
      (intersection.x - expected_x).abs() < f64::EPSILON,
      "X coord as expected"
    );
    assert!(
      (intersection.y - expected_y).abs() < f64::EPSILON,
      "Y coord as expected"
    );
  }

  #[test]
  fn test_v2d_almost_equals() {
    let v1 = Vec2d { x: 3.0, y: 4.0 };
    let v2 = Vec2d { x: 3.0, y: 4.0 };
    let v3 = Vec2d { x: 3.5, y: 4.0 };
    let v4 = Vec2d { x: 3.0, y: 4.5 };
    let v5 = Vec2d {
      x: 3.0 + f64::EPSILON / 2.0, // Closer than EPSILON
      y: 4.0 - f64::EPSILON / 2.0,
    };
    let v6 = Vec2d {
      x: 3.0 + f64::EPSILON * 2.0, // Farther than EPSILON
      y: 4.0,
    };

    assert!(_v2d_almost_equals(&v1, &v2), "true for equal vectors");
    assert!(!_v2d_almost_equals(&v1, &v3), "false for different x");
    assert!(!_v2d_almost_equals(&v1, &v4), "false for different y");
    assert!(_v2d_almost_equals(&v1, &v5), "true for almost equal within epsilon");
    assert!(!_v2d_almost_equals(&v1, &v6), "false for almost equal outside epsilon");
  }
}
