// src/polygon.rs

use crate::bbox::{
  bbox_contains_point, bbox_is_transmeridian, bbox_normalization, bbox_overlaps_bbox, LongitudeNormalization,
}; // Assuming bbox_contains_point and others are pub(crate) in bbox
use crate::constants::{EPSILON_RAD, M_2PI, M_PI}; // DBL_EPSILON is from C float.h, Rust f64::EPSILON is different
use crate::latlng::{constrain_lng, normalize_lng_for_comparison}; use crate::math::vec3d::_geo_to_vec3d;
// Assuming normalize_lng_for_comparison is pub(crate)
use crate::types::{BBox, CellBoundary, GeoLoop, GeoPolygon, LatLng};
use crate::{ContainmentMode, H3Error, Vec3d};

use std::f64;
use std::f64::EPSILON as MACHINE_EPSILON_F64; // For DBL_MAX

// Helper for 3D vector cross product
#[inline]
fn v3d_cross(v1: &Vec3d, v2: &Vec3d, out: &mut Vec3d) {
  out.x = v1.y * v2.z - v1.z * v2.y;
  out.y = v1.z * v2.x - v1.x * v2.z;
  out.z = v1.x * v2.y - v1.y * v2.x;
}

// Helper for 3D vector dot product
#[inline]
fn v3d_dot(v1: &Vec3d, v2: &Vec3d) -> f64 {
  v1.x * v2.x + v1.y * v2.y + v1.z * v2.z
}

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
pub(crate) fn generic_is_clockwise(loop_verts: &[LatLng], num_loop_verts: usize, is_transmeridian_bbox: bool) -> bool {
    if num_loop_verts < 3 {
        return false; // Not a polygon
    }

    let mut sum: f64 = 0.0;
    // Optional: Remove or comment out the println! statements after debugging
    // println!("--- generic_is_clockwise trace ---");
    // println!("is_transmeridian_bbox: {}", is_transmeridian_bbox);
    for i in 0..num_loop_verts {
        let p1_orig = loop_verts[i];
        let p2_orig = loop_verts[(i + 1) % num_loop_verts];

        let p1_lng = normalize_lng_for_comparison(
            p1_orig.lng,
            if is_transmeridian_bbox {
                LongitudeNormalization::East
            } else {
                LongitudeNormalization::None
            },
        );
        let p2_lng = normalize_lng_for_comparison(
            p2_orig.lng,
            if is_transmeridian_bbox {
                LongitudeNormalization::East
            } else {
                LongitudeNormalization::None
            },
        );
        
        let p1_lat = p1_orig.lat;
        let p2_lat = p2_orig.lat;

        let term = (p2_lng - p1_lng) * (p2_lat + p1_lat);
        sum += term;
        // Optional: Remove or comment out these println!
        // println!(
        //     "Iter {}: p1=({:.3},{:.3}), p2=({:.3},{:.3}) -> p1_lng={:.3}, p2_lng={:.3}, p1_lat={:.3}, p2_lat={:.3}",
        //     i, p1_orig.lat, p1_orig.lng, p2_orig.lat, p2_orig.lng, p1_lng, p2_lng, p1_lat, p2_lat
        // );
        // println!("  (p2_lng - p1_lng) = {:.3}", p2_lng - p1_lng);
        // println!("  (p2_lat + p1_lat) = {:.3}", p2_lat + p1_lat);
        // println!("  term = {:.3}, current sum = {:.3}", term, sum);
    }
    // Optional: Remove or comment out this println!
    // println!("Final sum: {:.3}", sum);
    
    // For the surveyor's formula variant (x2-x1)(y2+y1), a positive sum typically indicates CCW
    // traversal if y increases upwards and x increases to the right.
    // If lat=y and lng=x, our CCW test case yields sum=2.0.
    // Therefore, clockwise should yield a negative sum.
    sum < 0.0
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
/// Vertices must be in radians. Area is in radians^2.
/// Based on H3 C's `sphereArea` which uses summation of signed spherical triangle areas.
pub(crate) fn generic_area_rads2(verts_ll: &[LatLng]) -> f64 {
  let num_verts = verts_ll.len();
  if num_verts < 3 {
    return 0.0;
  }

  let mut total_angle_sum: f64 = 0.0;

  // Convert first vertex to 3D, use as anchor.
  let mut anchor_v3d = Vec3d::default();
  _geo_to_vec3d(&verts_ll[0], &mut anchor_v3d);

  // Iterate through triangles formed by (anchor, verts[i], verts[i+1])
  for i in 1..(num_verts - 1) {
    // Triangles are (v0, v1, v2), (v0, v2, v3), ...
    let mut v1_v3d = Vec3d::default();
    _geo_to_vec3d(&verts_ll[i], &mut v1_v3d);

    let mut v2_v3d = Vec3d::default();
    _geo_to_vec3d(&verts_ll[i + 1], &mut v2_v3d);

    // Calculate V = (anchor_v3d x v1_v3d) . v2_v3d  (scalar triple product)
    let mut cross_tmp = Vec3d::default();
    v3d_cross(&anchor_v3d, &v1_v3d, &mut cross_tmp);
    let v_val = v3d_dot(&cross_tmp, &v2_v3d);

    // Calculate S = 1 + (anchor . v1) + (v1 . v2) + (v2 . anchor)
    let s_val = 1.0 + v3d_dot(&anchor_v3d, &v1_v3d) + v3d_dot(&v1_v3d, &v2_v3d) + v3d_dot(&v2_v3d, &anchor_v3d);

    // atan2(V, S) gives half the signed area of the spherical triangle (on unit sphere)
    // or rather, it gives an angle related to the spherical excess.
    // The sum of these angles is related to the total spherical excess of the polygon.
    total_angle_sum += v_val.atan2(s_val);
  }

  // The total area is abs(total_angle_sum * 2.0)
  // The factor of 2 comes from the formula relating sum of atan2 to spherical excess.
  // The absolute value ensures positive area, as winding order might affect sign of sum.
  // H3 C's sphereArea returns fabs(total * 2.0).
  (total_angle_sum * 2.0).abs()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::bbox::bbox_from_geoloop; use crate::constants::M_PI_2;

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
    // println!("geoloop {:?}", geoloop_ccw);
    // println!("bbox {:?}", bbox_ccw);
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

  #[test]
  fn test_generic_area_rads2_simple_triangle() {
    // Define a simple spherical triangle, e.g., one octant of the sphere.
    // Vertices: (0,0), (PI/2, 0), (0, PI/2)
    let p1 = LatLng { lat: 0.0, lng: 0.0 };
    let p2 = LatLng { lat: M_PI_2, lng: 0.0 }; // North Pole on lng 0
    let p3 = LatLng { lat: 0.0, lng: M_PI_2 }; // Equator at 90 deg east

    let verts = [p1, p2, p3];
    let area = generic_area_rads2(&verts);

    // This triangle is 1/8 of the sphere's surface.
    // Sphere area is 4*PI. So, area should be (4*PI)/8 = PI/2.
    let expected_area = M_PI / 2.0;
    assert!(
      (area - expected_area).abs() < 1e-9,
      "Area of octant triangle. Got: {}, Expected: {}",
      area,
      expected_area
    );

    // Test with CCW square on equator
    let sq1 = LatLng { lat: 0.0, lng: 0.0 };
    let sq2 = LatLng { lat: 0.0, lng: M_PI_2 }; // 90 deg east
    let sq3 = LatLng {
      lat: M_PI_2,
      lng: M_PI_2,
    }; // No, this is wrong for a square
       // A square on a sphere is harder to define with equal area vs. angles
       // Let's use the H3 test case for cell 85283473fffffff if we can get its boundary easily
       // The CLI test value for 85283473fffffff (Res 5 hex) is 0.0000065310 rads^2
    let cell_for_area_test = crate::H3Index(0x85283473fffffff);
    let boundary = crate::indexing::cell_to_boundary(cell_for_area_test).unwrap();
    let area_from_boundary = generic_area_rads2(&boundary.verts[0..boundary.num_verts]);
    let expected_cli_area = 0.0000065310;
    assert!(
      (area_from_boundary - expected_cli_area).abs() < expected_cli_area * 0.01, // 1% tolerance for now
      "Area for cell 85283473fffffff. Got: {}, Expected from CLI: {}",
      area_from_boundary,
      expected_cli_area
    );
  }

  // More tests needed for cell_boundary_crosses_geoloop, cell_boundary_crosses_polygon,
  // and cell_boundary_inside_polygon, especially with transmeridian and hole cases.
  // These are complex and require careful setup of CellBoundary and GeoPolygon fixtures.
}
