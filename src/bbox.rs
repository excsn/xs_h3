use crate::constants::{EPSILON, EPSILON_RAD, MAX_H3_RES, M_2PI, M_PI, M_PI_2};
use crate::latlng::{constrain_lng, great_circle_distance_km};
use crate::types::BBox;
use crate::{CellBoundary, GeoLoop, GeoPolygon, H3Error, LatLng, MAX_CELL_BNDRY_VERTS};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LongitudeNormalization {
  None,
  East,
  West,
}

// Other bbox functions like bbox_width_rads, bbox_height_rads
pub(crate) fn bbox_width_rads(bbox: &BBox) -> f64 {
  if bbox_is_transmeridian(bbox) {
    bbox.east - bbox.west + M_2PI
  } else {
    bbox.east - bbox.west
  }
}

pub(crate) fn bbox_height_rads(bbox: &BBox) -> f64 {
  bbox.north - bbox.south
}

pub(crate) fn bbox_is_transmeridian(bbox: &BBox) -> bool {
  bbox.east < bbox.west
}

/// Get the center of a bounding box.
#[inline]
pub(crate) fn bbox_center(bbox: &BBox, center: &mut LatLng) {
  center.lat = (bbox.north + bbox.south) * 0.5;
  let east_for_center = if bbox_is_transmeridian(bbox) {
    bbox.east + M_2PI
  } else {
    bbox.east
  };
  center.lng = constrain_lng((east_for_center + bbox.west) * 0.5); // constrain_lng from latlng.rs
}

/// Whether the bounding box contains a given point.
#[inline]
#[must_use]
pub(crate) fn bbox_contains_point(bbox: &BBox, point: &LatLng) -> bool {
  if point.lat < bbox.south - EPSILON_RAD || point.lat > bbox.north + EPSILON_RAD {
    // Add epsilon for strict float compare
    return false;
  }
  if bbox_is_transmeridian(bbox) {
    // transmeridian case
    point.lng >= bbox.west - EPSILON_RAD || point.lng <= bbox.east + EPSILON_RAD
  } else {
    // standard case
    point.lng >= bbox.west - EPSILON_RAD && point.lng <= bbox.east + EPSILON_RAD
  }
}

/// Internal helper from C: _hexRadiusKm
/// This needs `get_pentagons` if we want to exactly match C logic of using a pentagon for max distortion.
/// For now, let's use a simpler approach or assume it's provided if needed by estimators.
/// If we port get_pentagons, it belongs in h3_index::inspection or similar.
/// For now, a placeholder for estimator dependency:
fn _hex_radius_km_approx(res: i32) -> Result<f64, H3Error> {
  // A very rough approximation based on average edge length.
  // This is NOT how C H3 does it and might not be suitable for accurate estimates.
  // The C version gets an actual pentagon at that resolution and calculates its radius.
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }
  Ok(crate::latlng::get_hexagon_edge_length_avg_km(res)? * 0.8) // Arbitrary factor
}

/// Determine the longitude normalization scheme for two bounding boxes.
pub(crate) fn bbox_normalization(
  a: &BBox,
  b: &BBox,
  a_normalization: &mut LongitudeNormalization,
  b_normalization: &mut LongitudeNormalization,
) {
  let a_is_transmeridian = bbox_is_transmeridian(a);
  let b_is_transmeridian = bbox_is_transmeridian(b);

  // If a trends east to overlap b, based on C logic: a->west - b->east < b->west - a->east
  // This means the distance from a's west edge to b's east edge (going east)
  // is less than distance from b's west edge to a's east edge (going east).
  // This check is trying to determine if normalizing `a` eastward and `b` westward
  // (or vice versa) would result in a smaller overall span or more direct overlap.
  let a_to_b_trends_east = (a.west - b.east).abs() < (b.west - a.east).abs(); // Simplified heuristic

  *a_normalization = if !a_is_transmeridian {
    LongitudeNormalization::None
  } else if b_is_transmeridian {
    LongitudeNormalization::East // Both transmeridian, normalize East by convention
  } else if a_to_b_trends_east {
    LongitudeNormalization::East
  } else {
    LongitudeNormalization::West
  };

  *b_normalization = if !b_is_transmeridian {
    LongitudeNormalization::None
  } else if a_is_transmeridian {
    LongitudeNormalization::East
  } else if a_to_b_trends_east {
    // If a trends east to b, b trends west to a for overlap
    LongitudeNormalization::West
  } else {
    LongitudeNormalization::East
  };
}

/// Helper to normalize a longitude based on the determined scheme.
fn normalize_lng_for_comparison(lng: f64, normalization: LongitudeNormalization) -> f64 {
  match normalization {
    LongitudeNormalization::None => lng,
    LongitudeNormalization::East => {
      if lng < 0.0 {
        lng + M_2PI
      } else {
        lng
      }
    }
    LongitudeNormalization::West => {
      if lng > 0.0 {
        lng - M_2PI
      } else {
        lng
      }
    }
  }
}

/// Whether one bounding box (a) contains another bounding box (b).
#[inline]
#[must_use]
pub(crate) fn bbox_contains_bbox(a: &BBox, b: &BBox) -> bool {
  if a.north < b.north || a.south > b.south {
    return false;
  }

  let mut a_norm = LongitudeNormalization::None;
  let mut b_norm = LongitudeNormalization::None;
  bbox_normalization(a, b, &mut a_norm, &mut b_norm);

  normalize_lng_for_comparison(a.west, a_norm) <= normalize_lng_for_comparison(b.west, b_norm)
    && normalize_lng_for_comparison(a.east, a_norm) >= normalize_lng_for_comparison(b.east, b_norm)
}

/// Whether two bounding boxes overlap.
#[inline]
#[must_use]
pub(crate) fn bbox_overlaps_bbox(a: &BBox, b: &BBox) -> bool {
  if a.north < b.south || a.south > b.north {
    // Corrected: use strict less/greater for non-overlap
    return false;
  }

  let mut a_norm = LongitudeNormalization::None;
  let mut b_norm = LongitudeNormalization::None;
  bbox_normalization(a, b, &mut a_norm, &mut b_norm);

  // If normalized east of A is west of normalized west of B, no overlap.
  // Or, if normalized west of A is east of normalized east of B, no overlap.
  if normalize_lng_for_comparison(a.east, a_norm) < normalize_lng_for_comparison(b.west, b_norm)
    || normalize_lng_for_comparison(a.west, a_norm) > normalize_lng_for_comparison(b.east, b_norm)
  {
    return false;
  }
  true
}

/// Whether two bounding boxes are strictly equal.
#[inline]
#[must_use]
pub(crate) fn bbox_equals(b1: &BBox, b2: &BBox) -> bool {
  // Using almost_equal for float comparisons
  (b1.north - b2.north).abs() < EPSILON_RAD &&
    (b1.south - b2.south).abs() < EPSILON_RAD &&
    (b1.east - b2.east).abs() < EPSILON_RAD && // This won't work well for transmeridian or normalized lng
    (b1.west - b2.west).abs() < EPSILON_RAD
  // A more robust check for longitude would involve normalizing them or checking width/center.
  // For strict H3 port, direct comparison like C is often used, assuming inputs are canonical.
  // However, for general utility, normalizing before comparison is better.
  // For now, direct float comparison (with epsilon).
}

/// Converts a `BBox` to a `CellBoundary` (4 vertices, CCW).
#[inline]
#[must_use]
pub(crate) fn bbox_to_cell_boundary(bbox: &BBox) -> CellBoundary {
  let mut vertices = [LatLng::default(); MAX_CELL_BNDRY_VERTS];

  // Define the 4 significant vertices of the bounding box.
  // Order them to be counter-clockwise (CCW) for consistency with H3 cell boundaries.
  // Standard CCW order for a rectangle: BL, BR, TR, TL
  vertices[0] = LatLng {
    lat: bbox.south,
    lng: bbox.west,
  }; // Bottom-Left (SW)
  vertices[1] = LatLng {
    lat: bbox.south,
    lng: bbox.east,
  }; // Bottom-Right (SE)
  vertices[2] = LatLng {
    lat: bbox.north,
    lng: bbox.east,
  }; // Top-Right (NE)
  vertices[3] = LatLng {
    lat: bbox.north,
    lng: bbox.west,
  }; // Top-Left (NW)

  CellBoundary {
    num_verts: 4,
    verts: vertices,
  }
}

/// Estimates the number of H3 cells needed to cover the given `BBox` at `res`.
pub(crate) fn bbox_hex_estimate(bbox: &BBox, res: i32) -> Result<i64, H3Error> {
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }

  // The C H3 logic uses a pentagon's radius for the most distorted case.
  // This requires `get_pentagons` to be ported.
  // For now, a simpler estimation using average area, which is less accurate
  // but allows us to proceed.
  // This is NOT a direct port of C's bboxHexEstimate if it relies on pentagon radius.

  let avg_cell_area_km2 = crate::latlng::get_hexagon_area_avg_km2(res)?;
  if avg_cell_area_km2 < EPSILON {
    return Ok(1);
  } // Avoid division by zero

  // Estimate bbox area roughly. This is very crude for spherical rectangles.
  let height_km = great_circle_distance_km(
    &LatLng {
      lat: bbox.north,
      lng: bbox.west,
    },
    &LatLng {
      lat: bbox.south,
      lng: bbox.west,
    },
  );
  let width_km = great_circle_distance_km(
    &LatLng {
      lat: bbox.south,
      lng: bbox.east,
    }, // Use south for width to get approx mid-latitude width
    &LatLng {
      lat: bbox.south,
      lng: bbox.west,
    },
  );

  let bbox_area_km2_approx = height_km * width_km;
  if bbox_area_km2_approx < EPSILON {
    return Ok(1);
  } // Avoid issues with tiny/degenerate bboxes

  let estimate = (bbox_area_km2_approx / avg_cell_area_km2).ceil();
  if !estimate.is_finite() {
    return Err(H3Error::Failed); // Or ResDomain if area calculation failed due to res
  }

  Ok(estimate.max(1.0) as i64) // Ensure at least 1
}

/// Estimates the number of H3 cells needed to trace the line from `origin` to `destination`.
pub(crate) fn line_hex_estimate(origin: &LatLng, destination: &LatLng, res: i32) -> Result<i64, H3Error> {
  if res < 0 || res > MAX_H3_RES {
    return Err(H3Error::ResDomain);
  }

  // Similar to bboxHexEstimate, C H3 uses pentagon radius.
  // Using average edge length for now.
  let avg_edge_len_km = crate::latlng::get_hexagon_edge_length_avg_km(res)?;
  if avg_edge_len_km < EPSILON {
    return Ok(1);
  }

  let line_dist_km = great_circle_distance_km(origin, destination);

  // Number of edges to span the distance. Diameter approx 2 * edge_len.
  let estimate = (line_dist_km / (2.0 * avg_edge_len_km)).ceil();
  if !estimate.is_finite() {
    return Err(H3Error::Failed);
  }
  Ok(estimate.max(1.0) as i64)
}

/// Create `BBox`es from a `GeoPolygon` (one for outer, one for each hole).
/// `bboxes_out` must be pre-allocated to `1 + polygon.num_holes`.
/// Create a bounding box from a GeoLoop.
/// This version is a more direct port of the single-pass C H3 logic.
pub(crate) fn bbox_from_geoloop(geoloop: &GeoLoop, bbox: &mut BBox) {
  if geoloop.num_verts == 0 {
    *bbox = BBox::default();
    return;
  }

  bbox.south = f64::MAX;
  bbox.west = f64::MAX; // Smallest positive longitude seen
  bbox.north = -f64::MAX;
  bbox.east = -f64::MAX; // Largest negative longitude seen (used if transmeridian)

  let mut min_lng_regular = f64::MAX; // Smallest longitude in normal [-PI, PI] sense
  let mut max_lng_regular = -f64::MAX; // Largest longitude in normal [-PI, PI] sense
  let mut crosses_antimeridian_arc = false;

  for i in 0..geoloop.num_verts {
    let p = geoloop.verts[i];

    if p.lat < bbox.south {
      bbox.south = p.lat;
    }
    if p.lat > bbox.north {
      bbox.north = p.lat;
    }

    // Track min/max longitude in the standard way first
    if p.lng < min_lng_regular {
      min_lng_regular = p.lng;
    }
    if p.lng > max_lng_regular {
      max_lng_regular = p.lng;
    }

    if i < geoloop.num_verts - 1 {
      // Check arc to next point
      let next_p = geoloop.verts[i + 1];
      if (p.lng - next_p.lng).abs() > M_PI {
        crosses_antimeridian_arc = true;
      }
    } else {
      // Check arc from last point to first point
      if geoloop.num_verts > 1 {
        // Avoid if only one point
        let first_p = geoloop.verts[0];
        if (p.lng - first_p.lng).abs() > M_PI {
          crosses_antimeridian_arc = true;
        }
      }
    }
  }

  if crosses_antimeridian_arc {
    // If any arc crosses, the simple min/max longitude is not sufficient.
    // We need to find the "tightest" box considering the antimeridian.
    // The west edge becomes the smallest positive longitude.
    // The east edge becomes the largest negative longitude.
    bbox.west = f64::MAX; // Re-init for positive longitudes
    bbox.east = -f64::MAX; // Re-init for negative longitudes
    let mut has_pos_lng = false;
    let mut has_neg_lng = false;

    for i in 0..geoloop.num_verts {
      let p_lng = geoloop.verts[i].lng;
      if p_lng > 0.0 {
        // Check positive hemisphere
        if p_lng < bbox.west {
          bbox.west = p_lng;
        }
        has_pos_lng = true;
      }
      if p_lng < 0.0 {
        // Check negative hemisphere
        if p_lng > bbox.east {
          bbox.east = p_lng;
        }
        has_neg_lng = true;
      }
      // If p_lng is exactly 0.0, it could be part of either segment
      // depending on its neighbors, but it won't change min_pos_lng or max_neg_lng.
    }
    // If all points were on one side (e.g. all positive, but an arc crossed),
    // one of these might not have been updated.
    if !has_pos_lng {
      bbox.west = bbox.east;
    } // Effectively, all points were negative or zero.
    if !has_neg_lng {
      bbox.east = bbox.west;
    } // Effectively, all points were positive or zero.
  } else {
    // No arc crossed the antimeridian. Use simple min/max.
    bbox.west = min_lng_regular;
    bbox.east = max_lng_regular;
  }
}

/// Scale a given bounding box by some factor.
pub(crate) fn scale_bbox(bbox: &mut BBox, scale: f64) {
  let width = bbox_width_rads(bbox);
  let height = bbox_height_rads(bbox);

  let width_buffer = (width * scale - width) * 0.5;
  let height_buffer = (height * scale - height) * 0.5;

  bbox.north += height_buffer;
  if bbox.north > M_PI_2 {
    bbox.north = M_PI_2;
  }
  bbox.south -= height_buffer;
  if bbox.south < -M_PI_2 {
    bbox.south = -M_PI_2;
  }

  bbox.east += width_buffer;
  bbox.east = constrain_lng(bbox.east);

  bbox.west -= width_buffer;
  bbox.west = constrain_lng(bbox.west);
}

/// Create `BBox`es from a `GeoPolygon` (one for outer, one for each hole).
/// `bboxes_out` must be pre-allocated to `1 + polygon.num_holes`.
pub(crate) fn bboxes_from_geo_polygon(polygon: &GeoPolygon, bboxes_out: &mut [BBox]) {
  if bboxes_out.is_empty() {
    return; // Nothing to do if output slice is empty
  }

  // Calculate BBox for the outer loop
  bbox_from_geoloop(&polygon.geoloop, &mut bboxes_out[0]);

  // Calculate BBox for each hole
  for i in 0..polygon.num_holes {
    // Ensure we don't write out of bounds of the provided bboxes_out slice
    if (i + 1) < bboxes_out.len() {
      bbox_from_geoloop(&polygon.holes[i], &mut bboxes_out[i + 1]);
    } else {
      // Not enough space in the output slice for all hole bboxes.
      // This indicates a mismatch between allocation and polygon.num_holes.
      // For robustness, we stop here. The C API often relies on the caller
      // to provide sufficiently sized buffers.
      break;
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::constants::{EPSILON_RAD, M_PI};
  use crate::latlng::{_set_geo_degs, geo_almost_equal};
  use crate::LatLng;

  #[test]
  fn test_bbox_width_height_rads() {
    let bbox = BBox {
      north: (1.0_f64).to_radians(),
      south: 0.0,
      east: (1.0_f64).to_radians(),
      west: 0.0,
    };
    assert!((bbox_width_rads(&bbox) - (1.0_f64).to_radians()).abs() < EPSILON_RAD);
    assert!((bbox_height_rads(&bbox) - (1.0_f64).to_radians()).abs() < EPSILON_RAD);

    let transmeridian_bbox = BBox {
      north: 0.1,
      south: -0.1,
      east: -M_PI + 0.2,
      west: M_PI - 0.2,
    }; // Spans 0.4 rads across antimeridian
    assert!(bbox_is_transmeridian(&transmeridian_bbox));
    assert!(
      (bbox_width_rads(&transmeridian_bbox) - (0.4)).abs() < EPSILON_RAD,
      "Transmeridian width. Got: {}",
      bbox_width_rads(&transmeridian_bbox)
    );
  }

  #[test]
  fn test_bbox_center() {
    let mut center = LatLng::default();
    let bbox1 = BBox {
      north: 1.0,
      south: 0.8,
      east: 1.0,
      west: 0.8,
    }; // Not radians, just for math
    let expected1 = LatLng { lat: 0.9, lng: 0.9 };
    bbox_center(&bbox1, &mut center);
    assert!(geo_almost_equal(&center, &expected1));

    // Transmeridian
    let bbox_tm = BBox {
      north: 0.1,
      south: -0.1,
      east: -M_PI + 0.1,
      west: M_PI - 0.1,
    }; // Center on antimeridian
    let expected_tm = LatLng { lat: 0.0, lng: -M_PI };
    bbox_center(&bbox_tm, &mut center);
    assert!(
      geo_almost_equal(&center, &expected_tm),
      "Center transmeridian. Got: {:?}, Expected: {:?}",
      center,
      expected_tm
    );
  }

  #[test]
  fn test_bbox_contains_point() {
    let bbox = BBox {
      north: 0.1,
      south: -0.1,
      east: 0.2,
      west: -0.2,
    };
    let inside = LatLng { lat: 0.0, lng: 0.0 };
    let outside_lat = LatLng { lat: 0.5, lng: 0.0 };
    let outside_lng = LatLng { lat: 0.0, lng: 0.5 };
    assert!(bbox_contains_point(&bbox, &inside));
    assert!(!bbox_contains_point(&bbox, &outside_lat));
    assert!(!bbox_contains_point(&bbox, &outside_lng));

    let trans_bbox = BBox {
      north: 0.1,
      south: -0.1,
      east: -M_PI + 0.1,
      west: M_PI - 0.1,
    }; // Crosses antimeridian
    let inside_tm_east = LatLng {
      lat: 0.0,
      lng: -M_PI + 0.05,
    };
    let inside_tm_west = LatLng {
      lat: 0.0,
      lng: M_PI - 0.05,
    };
    let outside_tm = LatLng { lat: 0.0, lng: 0.0 }; // Now outside this TM bbox
    assert!(bbox_contains_point(&trans_bbox, &inside_tm_east));
    assert!(bbox_contains_point(&trans_bbox, &inside_tm_west));
    assert!(!bbox_contains_point(&trans_bbox, &outside_tm));
  }

  // TODO: Add tests for bbox_contains_bbox, bbox_overlaps_bbox (with normal and transmeridian cases),
  // bbox_equals (with precision), bbox_to_cell_boundary, bbox_hex_estimate, line_hex_estimate,
  // bboxes_from_geo_polygon, scale_bbox.
  // These require careful setup of test cases.
}
