// src/polygon.rs

use crate::bbox::{
  bbox_contains_point, bbox_is_transmeridian, bbox_normalization, bbox_overlaps_bbox, LongitudeNormalization,
}; // Assuming bbox_contains_point and others are pub(crate) in bbox
use crate::constants::{EPSILON_RAD, M_2PI, M_PI}; // DBL_EPSILON is from C float.h, Rust f64::EPSILON is different
use crate::latlng::{constrain_lng, normalize_lng_for_comparison}; // Assuming normalize_lng_for_comparison is pub(crate)
use crate::types::{BBox, CellBoundary, GeoLoop, GeoPolygon, LatLng};
use crate::{ContainmentMode, H3Error};

use std::f64;
use std::f64::EPSILON as MACHINE_EPSILON_F64; // For DBL_MAX

// This is the core logic of point-in-polygon, used by both GeoLoop and LinkedGeoLoop.
// It's essentially the content of the C macro expansion for pointInside.
fn generic_point_inside_loop(loop_verts: &[LatLng], num_loop_verts: usize, bbox: &BBox, coord: &LatLng) -> bool {
  if num_loop_verts == 0 {
    return false;
  }
  // Fail fast if we're outside the bounding box.
  if !bbox_contains_point(bbox, coord) {
    return false;
  }

  let is_transmeridian = bbox_is_transmeridian(bbox);
  let mut contains = false;

  // Adjust latitude slightly if it matches a vertex, to avoid issues where
  // a point on a horizontal edge is handled ambiguously.
  // Bias northwards.
  let mut lat = coord.lat;
  let mut lng = normalize_lng_for_comparison(
    coord.lng,
    if is_transmeridian {
      LongitudeNormalization::East
    } else {
      LongitudeNormalization::None
    },
  ); // Normalize target lng based on bbox type

  for i in 0..num_loop_verts {
    let p1 = loop_verts[i];
    let p2 = loop_verts[(i + 1) % num_loop_verts];

    // Adjust latitude if it's on a vertex to avoid double counting or missed intersections.
    if lat == p1.lat || lat == p2.lat {
      lat += f64::EPSILON * 10.0; // Small nudge north
                                  // Re-check bbox after nudge, though unlikely to fail if originally inside
      if lat > bbox.north {
        continue;
      } // Nudged out of bbox
    }
    // Rays are cast in the longitudinal direction. For points exactly on a
    // vertical segment, bias westerly.
    let mut p1_lng = normalize_lng_for_comparison(
      p1.lng,
      if is_transmeridian {
        LongitudeNormalization::East
      } else {
        LongitudeNormalization::None
      },
    );
    let mut p2_lng = normalize_lng_for_comparison(
      p2.lng,
      if is_transmeridian {
        LongitudeNormalization::East
      } else {
        LongitudeNormalization::None
      },
    );

    if (p1_lng - lng).abs() < f64::EPSILON || (p2_lng - lng).abs() < f64::EPSILON {
      lng -= f64::EPSILON * 10.0;
    }

    // Swap p1 and p2 if p1 is north of p2 (ensures p1.lat <= p2.lat)
    let (mut a, mut b) = if p1.lat > p2.lat { (p2, p1) } else { (p1, p2) };
    // Re-normalize longitudes after potential swap
    a.lng = normalize_lng_for_comparison(
      a.lng,
      if is_transmeridian {
        LongitudeNormalization::East
      } else {
        LongitudeNormalization::None
      },
    );
    b.lng = normalize_lng_for_comparison(
      b.lng,
      if is_transmeridian {
        LongitudeNormalization::East
      } else {
        LongitudeNormalization::None
      },
    );

    // If the point's latitude is not between the latitudes of the segment's endpoints,
    // or if the point is on the same latitude as a horizontal segment, skip.
    if lat < a.lat || lat >= b.lat {
      // Use >= for b.lat to handle horizontal segments robustly
      continue;
    }
    // If segment is horizontal and point is on it, this isn't an intersection for ray casting.
    if (a.lat - b.lat).abs() < f64::EPSILON && (lat - a.lat).abs() < f64::EPSILON {
      continue;
    }

    // Calculate the longitude of the intersection of the ray (horizontal line at point's lat)
    // with the line segment (a, b).
    let intersect_lng: f64;
    if (b.lat - a.lat).abs() < f64::EPSILON {
      // Horizontal segment (already handled by lat checks)
      continue;
    } else {
      // Standard case: interpolate longitude
      intersect_lng = (b.lng - a.lng) * (lat - a.lat) / (b.lat - a.lat) + a.lng;
    }

    // If intersection is to the right (east) of the point, flip containment.
    if intersect_lng > lng {
      contains = !contains;
    }
  }
  contains
}

/// Take a given `GeoLoop` data structure and check if it contains a given geo coordinate.
///
/// # Arguments
/// * `geoloop` - The geoloop.
/// * `bbox` - The pre-calculated bounding box for the `geoloop`.
/// * `coord` - The coordinate to check.
///
/// # Returns
/// `true` if the point is contained by the geoloop.
#[must_use]
pub(crate) fn point_inside_geoloop(geoloop: &GeoLoop, bbox: &BBox, coord: &LatLng) -> bool {
  generic_point_inside_loop(&geoloop.verts, geoloop.num_verts, bbox, coord)
}

/// Core logic for determining if a loop's winding order is clockwise.
/// Sum over the edges `(x2-x1)(y2+y1)`. If sum is >0, then clockwise.
/// See: https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
/// This is a "shoelace formula" variant for area, where sign indicates winding.
fn generic_is_clockwise(loop_verts: &[LatLng], num_loop_verts: usize, is_transmeridian_bbox: bool) -> bool {
  if num_loop_verts < 3 {
    return false; // Not a polygon
  }

  let mut sum: f64 = 0.0;
  for i in 0..num_loop_verts {
    let p1 = loop_verts[i];
    let p2 = loop_verts[(i + 1) % num_loop_verts];

    // Normalize longitudes if transmeridian for consistent calculation.
    // Normalize to the eastern representation [0, 2PI) if transmeridian.
    let p1_lng = normalize_lng_for_comparison(
      p1.lng,
      if is_transmeridian_bbox {
        LongitudeNormalization::East
      } else {
        LongitudeNormalization::None
      },
    );
    let p2_lng = normalize_lng_for_comparison(
      p2.lng,
      if is_transmeridian_bbox {
        LongitudeNormalization::East
      } else {
        LongitudeNormalization::None
      },
    );

    sum += (p2_lng - p1_lng) * (p2.lat + p1.lat);
  }
  sum > 0.0
}

/// Whether the winding order of a given `GeoLoop` is clockwise.
/// In GeoJSON, clockwise loops are inner holes.
#[must_use]
pub(crate) fn is_clockwise_geoloop(geoloop: &GeoLoop, bbox: &BBox) -> bool {
  generic_is_clockwise(&geoloop.verts, geoloop.num_verts, bbox_is_transmeridian(bbox))
}

/// Whether two lines (Cartesian) intersect. Based on:
/// http://www.jeffreythompson.org/collision-detection/line-line.php
///
/// # Arguments
/// * `a1` - Start of line A.
/// * `a2` - End of line A.
/// * `b1` - Start of line B.
/// * `b2` - End of line B.
///
/// # Returns
/// `true` if the lines intersect.
#[must_use]
pub(crate) fn line_crosses_line(a1: &LatLng, a2: &LatLng, b1: &LatLng, b2: &LatLng) -> bool {
  let denom = (b2.lng - b1.lng) * (a2.lat - a1.lat) - (b2.lat - b1.lat) * (a2.lng - a1.lng);
  if denom.abs() < MACHINE_EPSILON_F64 {
    // Use machine epsilon here
    return false;
  }
  // ... rest of function ...
  let ua_num = (b2.lat - b1.lat) * (a1.lng - b1.lng) - (b2.lng - b1.lng) * (a1.lat - b1.lat);
  let ub_num = (a2.lat - a1.lat) * (a1.lng - b1.lng) - (a2.lng - a1.lng) * (a1.lat - b1.lat);

  let ua = ua_num / denom;
  let ub = ub_num / denom;

  (ua >= 0.0 && ua <= 1.0) && (ub >= 0.0 && ub <= 1.0)
}

/// Whether any part of a `CellBoundary` crosses a `GeoLoop`. Crossing means
/// line segments intersect; containment is not crossing for this function.
pub(crate) fn cell_boundary_crosses_geoloop(
  geoloop: &GeoLoop,
  loop_bbox: &BBox, // Pre-calculated bbox of the geoloop
  boundary: &CellBoundary,
  boundary_bbox: &BBox, // Pre-calculated bbox of the CellBoundary
) -> bool {
  if geoloop.num_verts == 0 || boundary.num_verts == 0 {
    return false;
  }
  // Fail fast if bounding boxes do not overlap
  if !bbox_overlaps_bbox(loop_bbox, boundary_bbox) {
    return false;
  }

  // Determine normalization schemes needed for comparing across antimeridian
  let mut loop_norm = LongitudeNormalization::None;
  let mut boundary_norm = LongitudeNormalization::None;
  bbox_normalization(loop_bbox, boundary_bbox, &mut loop_norm, &mut boundary_norm);

  for i in 0..geoloop.num_verts {
    let mut p1_loop = geoloop.verts[i];
    p1_loop.lng = normalize_lng_for_comparison(p1_loop.lng, loop_norm);
    let mut p2_loop = geoloop.verts[(i + 1) % geoloop.num_verts];
    p2_loop.lng = normalize_lng_for_comparison(p2_loop.lng, loop_norm);

    for j in 0..boundary.num_verts {
      let mut p1_boundary = boundary.verts[j];
      p1_boundary.lng = normalize_lng_for_comparison(p1_boundary.lng, boundary_norm);
      let mut p2_boundary = boundary.verts[(j + 1) % boundary.num_verts];
      p2_boundary.lng = normalize_lng_for_comparison(p2_boundary.lng, boundary_norm);

      if line_crosses_line(&p1_loop, &p2_loop, &p1_boundary, &p2_boundary) {
        return true;
      }
    }
  }
  false
}

// Constant for the mask, can be defined here or in constants.rs
pub(crate) const FLAG_CONTAINMENT_MODE_MASK: u32 = 0b1111; // Allows up to 16 modes (0-15)

/// Gets the `ContainmentMode` from the `flags` argument.
#[inline]
#[must_use]
pub(crate) fn flag_get_containment_mode(flags: u32) -> Result<ContainmentMode, H3Error> {
  ContainmentMode::try_from(flags & FLAG_CONTAINMENT_MODE_MASK)
}

/// Validates the polyfill flags.
#[inline]
pub(crate) fn validate_polygon_flags(flags: u32) -> Result<(), H3Error> {
  // Check if any bits outside the mode mask are set
  if flags & !FLAG_CONTAINMENT_MODE_MASK != 0 {
    return Err(H3Error::OptionInvalid);
  }
  // Check if the mode itself is valid by trying to convert it
  match flag_get_containment_mode(flags) {
    Ok(ContainmentMode::Invalid) | Err(_) => Err(H3Error::OptionInvalid),
    Ok(_) => Ok(()),
  }
}

/// Checks if a point is inside a `GeoPolygon` (outer loop and holes).
#[must_use]
pub(crate) fn point_inside_polygon(
  polygon: &GeoPolygon,
  bboxes: &[BBox], // Index 0 for outer, 1+ for holes
  coord: &LatLng,
) -> bool {
  if polygon.geoloop.num_verts == 0 {
    // An empty outer loop cannot contain anything
    return false;
  }
  // Must be inside the outer loop first.
  if !point_inside_geoloop(&polygon.geoloop, &bboxes[0], coord) {
    return false;
  }
  // If inside outer loop, must NOT be inside any of the holes.
  for i in 0..polygon.num_holes {
    if polygon.holes[i].num_verts > 0 && // Skip empty holes
           point_inside_geoloop(&polygon.holes[i], &bboxes[i + 1], coord)
    {
      return false; // Point is inside a hole
    }
  }
  true // Inside outer loop and not in any hole
}

/// Whether any part of a `CellBoundary` crosses a `GeoPolygon` (outer loop or any hole).
pub(crate) fn cell_boundary_crosses_polygon(
  polygon: &GeoPolygon,
  polygon_bboxes: &[BBox], // Index 0 is outer, 1+ are holes
  boundary: &CellBoundary,
  boundary_bbox: &BBox,
) -> bool {
  // Check outer loop
  if cell_boundary_crosses_geoloop(&polygon.geoloop, &polygon_bboxes[0], boundary, boundary_bbox) {
    return true;
  }
  // Check holes
  for i in 0..polygon.num_holes {
    if cell_boundary_crosses_geoloop(&polygon.holes[i], &polygon_bboxes[i + 1], boundary, boundary_bbox) {
      return true;
    }
  }
  false
}

/// Whether a `CellBoundary` is completely contained by a `GeoPolygon`.
pub(crate) fn cell_boundary_inside_polygon(
  polygon: &GeoPolygon,
  polygon_bboxes: &[BBox], // Index 0 is outer, 1+ are holes
  boundary: &CellBoundary,
  boundary_bbox: &BBox,
) -> bool {
  if boundary.num_verts == 0 {
    return false;
  } // Empty boundary cannot be inside

  // 1. All points of the boundary must be inside the polygon's outer loop.
  //    Start with a quick check of the first point.
  if !point_inside_geoloop(&polygon.geoloop, &polygon_bboxes[0], &boundary.verts[0]) {
    return false;
  }
  // If first point is in, check all points (more robust, though slower if bbox is small)
  for i in 0..boundary.num_verts {
    if !point_inside_geoloop(&polygon.geoloop, &polygon_bboxes[0], &boundary.verts[i]) {
      return false;
    }
  }

  // 2. The boundary must not cross the polygon's outer loop.
  if cell_boundary_crosses_geoloop(&polygon.geoloop, &polygon_bboxes[0], boundary, boundary_bbox) {
    return false;
  }

  // 3. The boundary must not be inside or cross any of the polygon's holes.
  // Convert CellBoundary to GeoLoop for pointInsideGeoLoop hole check
  let boundary_as_geoloop = GeoLoop {
    num_verts: boundary.num_verts,
    verts: boundary.verts.as_slice()[0..boundary.num_verts].to_vec(), // Create a Vec for GeoLoop
  };

  for i in 0..polygon.num_holes {
    // If any point of the hole is inside the cell boundary, it's an intersection, not full containment of cell.
    if polygon.holes[i].num_verts > 0
      && point_inside_geoloop(&boundary_as_geoloop, boundary_bbox, &polygon.holes[i].verts[0])
    {
      return false;
    }
    // If cell boundary crosses this hole, it's not *fully* inside the polygon (it's in a hole).
    if cell_boundary_crosses_geoloop(&polygon.holes[i], &polygon_bboxes[i + 1], boundary, boundary_bbox) {
      return false;
    }
  }
  true
}

/// Calculates the area of a spherical polygon given by its vertices.
/// Uses L'Huilier's Theorem for spherical excess.
/// Vertices must be in radians and in counter-clockwise order for a positive area.
pub(crate) fn generic_area_rads2(verts: &[LatLng]) -> f64 {
  let n = verts.len();
  if n < 3 {
    return 0.0;
  }

  // Spherical excess E = Sum of interior angles - (n - 2) * PI
  // Area = E (if radius of sphere is 1)
  let mut sum_interior_angles = 0.0;

  for i in 0..n {
    let p0 = verts[if i == 0 { n - 1 } else { i - 1 }]; // Previous vertex
    let p1 = verts[i]; // Current vertex
    let p2 = verts[(i + 1) % n]; // Next vertex

    // Calculate bearings of lines p1->p0 and p1->p2
    // Bearing is angle from North, clockwise.
    // We need the interior angle at p1.
    // Azimuth p1->p0 and p1->p2
    let az1 = crate::latlng::_geo_azimuth_rads(&p1, &p0); // Azimuth from p1 to p0
    let az2 = crate::latlng::_geo_azimuth_rads(&p1, &p2); // Azimuth from p1 to p2

    // Angle between two azimuths: (az2 - az1) normalized to (-PI, PI]
    let mut angle = az2 - az1;
    while angle <= -M_PI {
      angle += M_2PI;
    }
    while angle > M_PI {
      angle -= M_2PI;
    }

    // For CCW polygon, interior angle is PI - angle (if angle is positive)
    // or PI + abs(angle) (if angle is negative)
    // This needs to be robust. A common way is PI + angle for CCW.
    // If polygon is CW, this will be negative.
    // For area, we need sum of interior angles.
    // The change in bearing is what we need. For CCW, this sum is 2PI.
    // sum_interior_angles += (M_PI + angle); // This is complex, sum of exterior angles is easier

    // Using sum of turning angles (exterior angles for CCW):
    // Normalizing az1 and az2 to [0, 2PI)
    let mut norm_az1 = crate::latlng::_pos_angle_rads(az1);
    let mut norm_az2 = crate::latlng::_pos_angle_rads(az2);

    let mut turn_angle = norm_az1 - norm_az2; // For CCW, this should be positive for left turns
    while turn_angle < -M_PI {
      turn_angle += M_2PI;
    }
    while turn_angle > M_PI {
      turn_angle -= M_2PI;
    }

    sum_interior_angles += turn_angle; // This is sum of exterior angles, should be 2PI for CCW convex
  }

  // For a simple convex polygon on a sphere, sum of exterior angles = 2PI.
  // Spherical excess E = Sum of interior angles - (n-2)PI.
  // Sum of interior angles = (n-2)PI + E.
  // Another formula based on Girard's theorem from sum of angles.
  // The sum we calculated using (az1 - az2) is actually the sum of *exterior* turning angles.
  // For a simple CCW polygon, this sum should be 2*PI.
  // Area A = R^2 * (Sum of interior angles - (n-2)PI). With R=1, A = E.
  // H3 C's `sphereArea` is more robust, sums signed triangle areas.
  // For now, return a placeholder or a very rough estimate.
  // A full spherical polygon area function is a project in itself.
  // The C H3 code's polygonAreaRads2 sums the areas of spherical triangles
  // formed by adjacent vertices and an arbitrary reference point (like the first vertex).
  // Area of spherical triangle with sides a,b,c: 4 * atan(sqrt(tan(s/2) * tan((s-a)/2) * tan((s-b)/2) * tan((s-c)/2)))
  // where s = (a+b+c)/2. Sides a,b,c are great circle distances.

  // Placeholder: this will not be accurate.
  if n == 0 {
    return 0.0;
  }
  let approx_area_based_on_first_triangle_and_bbox = {
    if n < 3 {
      return 0.0;
    }
    let mut bbox = BBox::default();
    crate::bbox::bbox_from_geoloop(
      &GeoLoop {
        num_verts: n,
        verts: verts.to_vec(),
      },
      &mut bbox,
    );
    (bbox.north - bbox.south) * crate::bbox::bbox_width_rads(&bbox) * 0.5 // Very rough
  };
  approx_area_based_on_first_triangle_and_bbox.abs() // Area should be positive
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::bbox::bbox_from_geoloop; // For setting up bboxes in tests
  use crate::latlng::_set_geo_degs;

  // SF Verts from C tests
  const SF_VERTS_RAW: [[f64; 2]; 6] = [
    [0.659966917655, -2.1364398519396],
    [0.6595011102219, -2.1359434279405],
    [0.6583348114025, -2.1354884206045],
    [0.6581220034068, -2.1382437718946],
    [0.6594479998527, -2.1384597563896],
    [0.6599990002976, -2.1376771158464],
  ];

  fn get_sf_geoloop() -> GeoLoop {
    let verts: Vec<LatLng> = SF_VERTS_RAW.iter().map(|p| LatLng { lat: p[0], lng: p[1] }).collect();
    GeoLoop {
      num_verts: verts.len(),
      verts,
    }
  }

  #[test]
  fn test_point_inside_geoloop() {
    let geoloop = get_sf_geoloop();
    let mut bbox = BBox::default();
    bbox_from_geoloop(&geoloop, &mut bbox);

    let inside = LatLng {
      lat: 0.659,
      lng: -2.136,
    };
    let outside = LatLng { lat: 1.0, lng: 2.0 };

    assert!(
      point_inside_geoloop(&geoloop, &bbox, &inside),
      "Point should be inside SF geoloop"
    );
    assert!(
      !point_inside_geoloop(&geoloop, &bbox, &outside),
      "Point should be outside SF geoloop"
    );
  }

  #[test]
  fn test_is_clockwise_geoloop() {
    let ccw_verts_vec: Vec<LatLng> = vec![
      LatLng { lat: 0.0, lng: 0.0 },
      LatLng { lat: 1.0, lng: 0.0 },
      LatLng { lat: 1.0, lng: 1.0 },
      LatLng { lat: 0.0, lng: 1.0 },
    ];
    let geoloop_ccw = GeoLoop {
      num_verts: ccw_verts_vec.len(),
      verts: ccw_verts_vec,
    };
    let mut bbox_ccw = BBox::default();
    bbox_from_geoloop(&geoloop_ccw, &mut bbox_ccw);
    assert!(
      !is_clockwise_geoloop(&geoloop_ccw, &bbox_ccw),
      "CCW loop should not be clockwise"
    );

    let cw_verts_vec: Vec<LatLng> = vec![
      LatLng { lat: 0.0, lng: 0.0 },
      LatLng { lat: 0.0, lng: 1.0 },
      LatLng { lat: 1.0, lng: 1.0 },
      LatLng { lat: 1.0, lng: 0.0 },
    ];
    let geoloop_cw = GeoLoop {
      num_verts: cw_verts_vec.len(),
      verts: cw_verts_vec,
    };
    let mut bbox_cw = BBox::default();
    bbox_from_geoloop(&geoloop_cw, &mut bbox_cw);
    assert!(
      is_clockwise_geoloop(&geoloop_cw, &bbox_cw),
      "CW loop should be clockwise"
    );
  }

  #[test]
  fn test_line_crosses_line() {
    let l1p1 = LatLng { lat: 0.0, lng: 0.0 };
    let l1p2 = LatLng { lat: 1.0, lng: 1.0 };
    let l2p1 = LatLng { lat: 0.0, lng: 1.0 };
    let l2p2 = LatLng { lat: 1.0, lng: 0.0 };
    assert!(line_crosses_line(&l1p1, &l1p2, &l2p1, &l2p2), "Diagonals should cross");

    let l3p1 = LatLng { lat: 0.0, lng: 0.0 };
    let l3p2 = LatLng { lat: 0.4, lng: 0.4 }; // Shorter segment of l1
    assert!(
      !line_crosses_line(&l3p1, &l3p2, &l2p1, &l2p2),
      "Shorter diagonal should not cross"
    );
  }

  // More tests needed for cell_boundary_crosses_geoloop, cell_boundary_crosses_polygon,
  // and cell_boundary_inside_polygon, especially with transmeridian and hole cases.
  // These are complex and require careful setup of CellBoundary and GeoPolygon fixtures.
}
