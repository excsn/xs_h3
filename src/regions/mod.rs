// src/regions/mod.rs
#![allow(clippy::module_name_repetitions)]

pub mod polyfill;
pub mod to_polygon; // For cellsToMultiPolygon and its helpers

// Re-export public API functions
pub use polyfill::{
  max_polygon_to_cells_size,
  polygon_to_cells, // These would be the experimental versions
                    // Or, if keeping C API names:
                    // maxPolygonToCellsSize, polygonToCells,
                    // maxPolygonToCellsSizeExperimental, polygonToCellsExperimental
};
pub use to_polygon::{cells_to_multi_polygon, MultiPolygonRust, PolygonRust};
