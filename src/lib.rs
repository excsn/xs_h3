#![deny(clippy::all)] // Enforce clippy lints
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)] // Often a matter of taste
#![allow(clippy::missing_errors_doc)] // TODO: Add error docs later
#![allow(clippy::cast_possible_truncation)] // Common in C ports, review carefully
#![allow(clippy::cast_precision_loss)] // Common in C ports, review carefully
#![allow(clippy::cast_sign_loss)] // Common in C ports, review carefully
#![allow(clippy::must_use_candidate)] // For functions where side effects are intended
#![allow(clippy::unreadable_literal)] // Sometimes C constants are clearer this way
#![allow(clippy::similar_names)] // Can be common in math-heavy code
#![allow(clippy::wildcard_imports)] // Allow for re-exporting from modules

//! `xs-h3` is a Rust implementation of Uber's H3 geospatial indexing library.
//!
//! This library aims to provide a safe and performant H3 experience in Rust.

// Declare modules
pub mod base_cells;
pub mod bbox;
pub mod constants;
pub mod coords;
pub mod h3_index;
pub mod hierarchy;
pub mod indexing;
pub mod iterators;
pub mod latlng;
pub mod local_ij;
pub mod math;
pub mod measures;
pub mod polygon;
pub mod regions;
pub mod traversal;
pub mod types;

// Re-export key public types and functions for easier use
pub use constants::MAX_CELL_BNDRY_VERTS;
pub use latlng::{
  degs_to_rads, get_hexagon_area_avg_km2, get_hexagon_area_avg_m2, get_hexagon_edge_length_avg_km,
  get_hexagon_edge_length_avg_m, great_circle_distance_km, great_circle_distance_m, great_circle_distance_rads,
  rads_to_degs,
};
pub use types::{
  BBox, CellBoundary, ContainmentMode, CoordIJ, CoordIJK, FaceIJK, GeoLoop, GeoPolygon, H3Error, H3Index, LatLng,
  Vec2d, Vec3d, H3_NULL,
};

pub use h3_index::inspection::{
  get_base_cell_number,
  get_icosahedron_faces,
  get_num_cells,
  get_pentagons,
  get_res0_cells,
  is_pentagon,
  is_valid_cell,
  max_face_count,
  pentagon_count,
};
pub use h3_index::string_conv::{h3_to_string, h3_to_string_alloc, string_to_h3};
pub use hierarchy::{
  cell_to_center_child, cell_to_child_pos, cell_to_children, cell_to_children_size, cell_to_parent, child_pos_to_cell,
  compact_cells, uncompact_cells, uncompact_cells_size,
};
pub use indexing::{cell_to_boundary, cell_to_lat_lng, lat_lng_to_cell};
pub use measures::{
  cell_area_km2,
  cell_area_m2,
  cell_area_rads2,
  exact_edge_length_km,
  exact_edge_length_m,
  exact_edge_length_rads,
};
pub use regions::{cells_to_multi_polygon, max_polygon_to_cells_size, polygon_to_cells};
pub use traversal::{
  are_neighbor_cells, grid_disk, grid_disk_distances, grid_distance, grid_path_cells, grid_path_cells_size,
  grid_ring_unsafe, max_grid_disk_size,
};
