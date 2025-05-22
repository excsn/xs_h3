# XS-H3 Usage Guide

This guide provides a detailed overview of how to use the `xs-h3` library, including core concepts, common API usage, and examples for various functionalities.

## Table of Contents

*   [Core Concepts](#core-concepts)
*   [Quick Start Examples](#quick-start-examples)
*   [Main API Sections / Functional Areas](#main-api-sections--functional-areas)
    *   [Indexing: Coordinates and H3 Cells](#indexing-coordinates-and-h3-cells)
    *   [Inspection: Understanding H3 Indexes](#inspection-understanding-h3-indexes)
    *   [Hierarchy: Parents, Children, and Compaction](#hierarchy-parents-children-and-compaction)
    *   [Traversal: Navigating the H3 Grid](#traversal-navigating-the-h3-grid)
    *   [Regions: Working with Polygons (Polyfill)](#regions-working-with-polygons-polyfill)
    *   [Measurements: Cell and Edge Metrics](#measurements-cell-and-edge-metrics)
    *   [Utilities: Geometric and String Conversions](#utilities-geometric-and-string-conversions)
*   [Data Types](#data-types)
*   [Error Handling](#error-handling)
*   [Optional Features](#optional-features)
    *   [Serde Support](#serde-support)

## Core Concepts

Before diving into the API, understanding these core H3 concepts will be helpful:

*   **`H3Index`**: This is the fundamental type in `xs-h3`, representing a unique H3 cell (or sometimes a directed edge) on the Earth's surface. It's internally a 64-bit unsigned integer. All H3 operations primarily revolve around manipulating these indexes.
*   **`LatLng`**: A simple struct holding `latitude` and `longitude` values, always expressed in **radians** for consistency with the H3 C library and mathematical formulas. Helper functions are provided to convert between degrees and radians.
*   **Resolution**: H3 cells are organized into 16 discrete resolution levels (0 to 15). Resolution 0 represents the coarsest cells (122 base cells covering the globe), while resolution 15 represents the finest (sub-meter precision). Most functions that operate on cells require a resolution parameter.
*   **Hexagonal Grid**: H3 uses a primarily hexagonal grid. This provides more uniform adjacencies and distances compared to square grids.
*   **Pentagons**: Due to the icosahedral base of the H3 system, each resolution level includes exactly 12 pentagonal cells. These pentagons have slightly different properties (e.g., 5 neighbors instead of 6 for their children) and require special handling in some algorithms.
*   **Error Handling**: Most functions in `xs-h3` that can encounter issues (e.g., invalid input, domain errors) return a `Result<T, H3Error>`. The `H3Error` enum provides specific error codes.

## Quick Start Examples

### Example 1: Convert Lat/Lng to H3 Cell and Back

This example shows the most common use case: finding the H3 cell for a given geographic coordinate and then finding the center of that cell.

```rust
use xs_h3::{
    lat_lng_to_cell, cell_to_lat_lng, get_resolution, 
    H3Index, LatLng, degs_to_rads, rads_to_degs, H3Error
};

fn main() -> Result<(), H3Error> {
    // Define a point (e.g., somewhere in San Francisco)
    let lat_deg = 37.7749;
    let lng_deg = -122.4194;

    let geo_point = LatLng {
        lat: degs_to_rads(lat_deg),
        lng: degs_to_rads(lng_deg),
    };

    // Choose a resolution (0-15)
    let resolution = 9;

    // Get the H3 cell index
    let h3_cell: H3Index = lat_lng_to_cell(&geo_point, resolution)?;
    println!("H3 Cell for ({}, {}) at res {}: {:x} (hex: {})", 
             lat_deg, lng_deg, resolution, h3_cell.0, xs_h3::h3_to_string_alloc(h3_cell));

    // Get the center coordinates of this H3 cell
    let center_coords_rad: LatLng = cell_to_lat_lng(h3_cell)?;
    println!(
        "Center of H3 cell {:x}: (lat: {:.6}, lng: {:.6}) degrees",
        h3_cell.0,
        rads_to_degs(center_coords_rad.lat),
        rads_to_degs(center_coords_rad.lng)
    );

    // Get the boundary of the H3 cell
    let boundary: xs_h3::CellBoundary = cell_to_boundary(h3_cell)?;
    println!("Boundary for cell {:x} has {} vertices:", h3_cell.0, boundary.num_verts);
    for i in 0..boundary.num_verts {
        println!(
            "  Vertex {}: (lat: {:.6}, lng: {:.6}) degrees",
            i,
            rads_to_degs(boundary.verts[i].lat),
            rads_to_degs(boundary.verts[i].lng)
        );
    }
    Ok(())
}
```

### Example 2: Get Neighbors (k-Ring)

This example demonstrates how to find all H3 cells within a certain grid distance (k-ring) of an origin cell.

```rust
use xs_h3::{lat_lng_to_cell, grid_disk, max_grid_disk_size, H3Index, LatLng, degs_to_rads, H3Error, H3_NULL};

fn main() -> Result<(), H3Error> {
    let origin_geo = LatLng { lat: degs_to_rads(40.6892), lng: degs_to_rads(-74.0445) }; // Statue of Liberty
    let resolution = 8;
    let origin_cell = lat_lng_to_cell(&origin_geo, resolution)?;

    let k = 1; // Get immediate neighbors (k-ring of 1)
    
    // Determine the maximum possible size for the k-ring output
    let max_size = max_grid_disk_size(k)? as usize;
    let mut neighbors = vec![H3_NULL; max_size]; // Initialize with H3_NULL

    grid_disk(origin_cell, k, &mut neighbors)?;

    println!("Neighbors (k={}) of cell {:x}:", k, origin_cell.0);
    for neighbor_cell in neighbors.iter().filter(|&&h| h != H3_NULL) {
        println!("  - {:x}", neighbor_cell.0);
    }
    Ok(())
}
```

## Main API Sections / Functional Areas

### Indexing: Coordinates and H3 Cells

These functions handle the fundamental conversions between geographic coordinates and H3 cell indexes.

*   `pub fn lat_lng_to_cell(geo: &LatLng, res: i32) -> Result<H3Index, H3Error>`
    *   Finds the H3 cell index containing the given `LatLng` point at the specified resolution.
*   `pub fn cell_to_lat_lng(cell: H3Index) -> Result<LatLng, H3Error>`
    *   Calculates the geographic coordinates (latitude/longitude in radians) of the center of the given H3 cell.
*   `pub fn cell_to_boundary(cell: H3Index) -> Result<CellBoundary, H3Error>`
    *   Determines the geographic boundary of an H3 cell, returning a `CellBoundary` struct containing the vertices.

### Inspection: Understanding H3 Indexes

Functions to get properties of H3 indexes or validate them.

*   `pub fn get_resolution(h: H3Index) -> i32`
    *   Returns the resolution (0-15) of the given H3 index.
*   `pub fn get_base_cell_number(h: H3Index) -> i32`
    *   Returns the base cell number (0-121) of the H3 index.
*   `pub fn is_valid_cell(h: H3Index) -> bool`
    *   Checks if the given `H3Index` is a valid H3 cell index.
*   `pub fn is_pentagon(h: H3Index) -> bool`
    *   Determines if the given H3 cell is one of the 12 pentagons at its resolution.
*   `pub fn is_res_class_iii(h: H3Index) -> bool`
    *   Checks if the resolution of the H3 cell is a Class III resolution (odd resolutions).
*   `pub fn get_num_cells(res: i32) -> Result<i64, H3Error>`
    *   Calculates the total number of unique H3 cells at the given resolution.
*   `pub fn pentagon_count() -> i32`
    *   Returns the number of pentagons at any H3 resolution (always 12).
*   `pub fn get_pentagons(res: i32, out_pentagons: &mut [H3Index; 12]) -> Result<(), H3Error>`
    *   Fills the `out_pentagons` array with the H3 indexes of all 12 pentagons at the specified resolution.
*   `pub fn get_res0_cells(out_res0_cells: &mut [H3Index; xs_h3::constants::NUM_BASE_CELLS])`
    *   Fills `out_res0_cells` with all 122 resolution 0 base cell H3 indexes.
*   `pub fn get_icosahedron_faces(h: H3Index, out_faces: &mut [i32]) -> Result<usize, H3Error>`
    *   Finds all icosahedron faces (0-19) intersected by a given H3 cell. The `out_faces` slice should typically be sized to 2, as a cell can cross at most one icosahedron edge.
*   `pub fn max_face_count(_h: H3Index) -> Result<usize, H3Error>`
    *   Returns the maximum number of icosahedron faces an H3 cell's boundary may cross (typically 2).

### Hierarchy: Parents, Children, and Compaction

Functions for navigating the H3 hierarchy and optimizing sets of cells.

*   `pub fn cell_to_parent(h: H3Index, parent_res: i32) -> Result<H3Index, H3Error>`
    *   Returns the H3 parent index of cell `h` at the specified `parent_res`. `parent_res` must be coarser than or equal to `h`'s resolution.
*   `pub fn cell_to_children_size(h: H3Index, child_res: i32) -> Result<i64, H3Error>`
    *   Calculates the number of children cell `h` has at `child_res`. `child_res` must be finer than or equal to `h`'s resolution.
*   `pub fn cell_to_children(h: H3Index, child_res: i32, children: &mut [H3Index]) -> Result<(), H3Error>`
    *   Populates the `children` slice with all children of cell `h` at `child_res`. The slice must be pre-sized correctly using `cell_to_children_size`.
*   `pub fn cell_to_center_child(h: H3Index, child_res: i32) -> Result<H3Index, H3Error>`
    *   Returns the center child of cell `h` at `child_res`.
*   `pub fn cell_to_child_pos(child: H3Index, parent_res: i32) -> Result<i64, H3Error>`
    *   Returns the position (0-indexed) of the `child` cell within an ordered list of all children of its conceptual parent at `parent_res`.
*   `pub fn child_pos_to_cell(child_pos: i64, parent: H3Index, child_res: i32) -> Result<H3Index, H3Error>`
    *   Returns the H3 child cell at `child_pos` relative to `parent`, at `child_res`.
*   `pub fn compact_cells(h3_set: &mut [H3Index], out_compacted_set: &mut [H3Index]) -> Result<usize, H3Error>`
    *   Compacts a set of H3 cells of the same resolution into a smaller set of cells covering the same area. `h3_set` may be sorted in place. Returns the number of cells written to `out_compacted_set`.
*   `pub fn uncompact_cells_size(compacted_set: &[H3Index], res: i32) -> Result<i64, H3Error>`
    *   Calculates the number of cells that will result from uncompacting `compacted_set` to resolution `res`.
*   `pub fn uncompact_cells(compacted_set: &[H3Index], res: i32, out_set: &mut [H3Index]) -> Result<(), H3Error>`
    *   Uncompacts `compacted_set` to resolution `res`, filling `out_set`.

### Traversal: Navigating the H3 Grid

Functions for finding neighbors, paths, and disks of cells.

*   `pub fn grid_distance(origin: H3Index, destination: H3Index) -> Result<i64, H3Error>`
    *   Calculates the grid distance (number of H3 cell steps) between `origin` and `destination`. They must be at the same resolution.
*   `pub fn are_neighbor_cells(origin: H3Index, destination: H3Index) -> Result<bool, H3Error>`
    *   Determines if two H3 cells are immediate neighbors.
*   `pub fn max_grid_disk_size(k: i32) -> Result<i64, H3Error>`
    *   Calculates the maximum number of cells that can be in a k-ring (hexagonal disk) of radius `k`.
*   `pub fn grid_disk(origin: H3Index, k: i32, out_cells: &mut [H3Index]) -> Result<(), H3Error>`
    *   Populates `out_cells` with all H3 cells within grid distance `k` of `origin` (inclusive). Handles pentagons correctly.
*   `pub fn grid_disk_distances(origin: H3Index, k: i32, out_cells: &mut [H3Index], out_distances_opt: Option<&mut [i32]>) -> Result<(), H3Error>`
    *   Similar to `grid_disk`, but also populates `out_distances_opt` with the grid distance of each cell from the origin.
*   `pub fn grid_ring_unsafe(origin: H3Index, k: i32, out_cells: &mut [H3Index]) -> Result<(), H3Error>`
    *   Produces the "hollow" ring of cells exactly `k` distance from `origin`. Marked "unsafe" as its behavior is undefined if a pentagon is encountered during traversal.
*   `pub fn grid_path_cells_size(start: H3Index, end: H3Index) -> Result<i64, H3Error>`
    *   Calculates the number of cells in the line from `start` to `end` (inclusive).
*   `pub fn grid_path_cells(start: H3Index, end: H3Index, out_path: &mut [H3Index]) -> Result<(), H3Error>`
    *   Populates `out_path` with the H3 cells forming a line between `start` and `end` (inclusive).

### Regions: Working with Polygons (Polyfill)

These functions are used to find H3 cells that cover a given polygon. This functionality is considered **experimental**.

*   `pub fn max_polygon_to_cells_size(polygon: &GeoPolygon, res: i32, flags: u32) -> Result<i64, H3Error>`
    *   Estimates the maximum number of H3 cells required to cover the given `polygon` at resolution `res` using the specified `flags`.
*   `pub fn polygon_to_cells(polygon: &GeoPolygon, res: i32, flags: u32, out: &mut [H3Index]) -> Result<(), H3Error>`
    *   Fills the `out` slice with H3 cells at resolution `res` that cover the `polygon`. The `flags` argument controls the containment mode (e.g., `ContainmentMode::Center`).
*   `pub fn cells_to_multi_polygon(h3_set: &[H3Index]) -> Result<MultiPolygonRust, H3Error>`
    *   Converts a set of H3 cells (presumably contiguous or forming distinct groups) into a `MultiPolygonRust` structure, which is a `Vec` of `PolygonRust` objects. Each `PolygonRust` has an outer boundary and a `Vec` of hole boundaries. Hole detection and full normalization are currently basic.

### Measurements: Cell and Edge Metrics

Functions for calculating geometric properties of H3 cells and edges.

*   `pub fn cell_area_rads2(cell: H3Index) -> Result<f64, H3Error>`
    *   Calculates the area of the H3 cell in square radians.
*   `pub fn cell_area_km2(cell: H3Index) -> Result<f64, H3Error>`
    *   Calculates the area of the H3 cell in square kilometers.
*   `pub fn cell_area_m2(cell: H3Index) -> Result<f64, H3Error>`
    *   Calculates the area of the H3 cell in square meters.
*   `pub fn get_hexagon_area_avg_km2(res: i32) -> Result<f64, H3Error>`
    *   Returns the average area of an H3 hexagon (excluding pentagons) at the given resolution in square kilometers.
*   `pub fn get_hexagon_area_avg_m2(res: i32) -> Result<f64, H3Error>`
    *   Returns the average area of an H3 hexagon at the given resolution in square meters.
*   `pub fn get_hexagon_edge_length_avg_km(res: i32) -> Result<f64, H3Error>`
    *   Returns the average edge length of an H3 hexagon at the given resolution in kilometers.
*   `pub fn get_hexagon_edge_length_avg_m(res: i32) -> Result<f64, H3Error>`
    *   Returns the average edge length of an H3 hexagon at the given resolution in meters.
*   `pub fn exact_edge_length_rads(edge: H3Index) -> Result<f64, H3Error>`
    *   Calculates the exact great circle distance of a directed H3 edge in radians. (Currently a placeholder, requires directed edge support).
*   `pub fn exact_edge_length_km(edge: H3Index) -> Result<f64, H3Error>`
    *   Calculates the exact great circle distance of a directed H3 edge in kilometers. (Currently a placeholder).
*   `pub fn exact_edge_length_m(edge: H3Index) -> Result<f64, H3Error>`
    *   Calculates the exact great circle distance of a directed H3 edge in meters. (Currently a placeholder).

### Utilities: Geometric and String Conversions

General utility functions.

*   `pub fn degs_to_rads(degrees: f64) -> f64`
    *   Converts an angle from degrees to radians.
*   `pub fn rads_to_degs(radians: f64) -> f64`
    *   Converts an angle from radians to degrees.
*   `pub fn great_circle_distance_rads(a: &LatLng, b: &LatLng) -> f64`
    *   Calculates the great circle (haversine) distance between two `LatLng` points in radians.
*   `pub fn great_circle_distance_km(a: &LatLng, b: &LatLng) -> f64`
    *   Calculates the great circle distance in kilometers.
*   `pub fn great_circle_distance_m(a: &LatLng, b: &LatLng) -> f64`
    *   Calculates the great circle distance in meters.
*   `pub fn string_to_h3(s: &str) -> Result<H3Index, H3Error>`
    *   Converts a hexadecimal string representation to an `H3Index`.
*   `pub fn h3_to_string(h: H3Index, buffer: &mut [u8]) -> Result<(), H3Error>`
    *   Converts an `H3Index` to its hexadecimal string representation, storing it in the provided `buffer`. The buffer should be at least 17 bytes long for a 16-character hex string + null terminator.
*   `pub fn h3_to_string_alloc(h: H3Index) -> String`
    *   Converts an `H3Index` to an allocated `String`.

## Data Types

The library uses several key structs to represent H3 concepts:

*   **`xs_h3::H3Index`**: The core 64-bit unsigned integer for H3 cell/edge indexes.
*   **`xs_h3::LatLng`**: `struct { lat: f64, lng: f64 }` (coordinates in radians).
*   **`xs_h3::CellBoundary`**: `struct { num_verts: usize, verts: [LatLng; MAX_CELL_BNDRY_VERTS] }`.
*   **`xs_h3::GeoLoop`**: `struct { num_verts: usize, verts: Vec<LatLng> }`.
*   **`xs_h3::GeoPolygon`**: `struct { geoloop: GeoLoop, num_holes: usize, holes: Vec<GeoLoop> }`.
*   **`xs_h3::PolygonRust`**: `struct { outer: Vec<LatLng>, holes: Vec<Vec<LatLng>> }` (used by `cells_to_multi_polygon`).
*   **`xs_h3::BBox`**: `struct { north: f64, south: f64, east: f64, west: f64 }`.
*   **`xs_h3::Direction`**: Enum for H3 digit directions (0-7).
*   **`xs_h3::H3Error`**: Enum for error codes.

Refer to the [API Reference documentation](https://docs.rs/xs-h3) for detailed struct field descriptions.

## Error Handling

Most functions that can encounter an error return a `Result<T, xs_h3::H3Error>`. The `H3Error` enum lists possible error conditions:

*   `Success`: No error.
*   `Failed`: Generic failure.
*   `Domain`: Argument out of its valid domain.
*   `LatLngDomain`: Latitude/longitude out of range.
*   `ResDomain`: Resolution out of range (0-15).
*   `CellInvalid`: Provided H3 cell index is invalid.
*   `DirEdgeInvalid`: Provided H3 directed edge index is invalid.
*   `UndirEdgeInvalid`: Provided H3 undirected edge index is invalid.
*   `VertexInvalid`: Provided H3 vertex index is invalid.
*   `Pentagon`: Algorithm encountered a pentagon requiring special handling that was not supported or led to an issue.
*   `DuplicateInput`: Duplicate H3 indexes provided where not allowed.
*   `NotNeighbors`: Cells provided are not neighbors as expected.
*   `ResMismatch`: Cells provided have incompatible resolutions.
*   `MemoryAlloc`: Failed to allocate necessary memory.
*   `MemoryBounds`: Provided output buffer was too small.
*   `OptionInvalid`: Invalid mode or flags argument.

## Optional Features

### Serde Support

To enable serialization and deserialization for `xs-h3` types like `H3Index`, `LatLng`, `H3Error`, etc., enable the `serde` feature in your `Cargo.toml`:

```toml
[dependencies]
xs-h3 = { version = "0.1.0", features = ["serde"] }
```

This will allow you to use Serde's `Serialize` and `Deserialize` traits with these types. For enums like `H3Error` and `Direction` that have a C-style integer representation, `serde_repr` is used to ensure they serialize/deserialize to/from their integer values.