// src/traversal/mod.rs
#![allow(clippy::module_name_repetitions)]

// Declare submodules as we create them
pub mod grid_disk;
pub mod grid_path;
pub mod distance;
pub mod neighbors; // For h3NeighborRotations, directionForNeighbor, are_neighbor_cells

// Re-export public API functions

pub use grid_disk::{grid_disk, grid_disk_distances, max_grid_disk_size, grid_ring_unsafe}; // grid_ring_unsafe might be internal only
pub use neighbors::are_neighbor_cells;
pub use distance::grid_distance;
pub use grid_path::{grid_path_cells, grid_path_cells_size};