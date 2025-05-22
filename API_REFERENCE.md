# xs-h3 API Reference

This document provides an API reference for the `xs-h3` crate, a Rust implementation of Uber's H3 geospatial indexing system.

## 1. Core Concepts

The `xs-h3` library revolves around a few key concepts:

*   **`H3Index`**: The primary type representing an H3 cell index or a directed H3 edge index. It's a 64-bit unsigned integer (`u64`). Most operations in the library consume or produce `H3Index` values.
*   **`LatLng`**: Represents a geographic coordinate with latitude and longitude, stored in radians.
*   **Resolutions**: H3 cells exist at 16 different resolutions (0 to 15), with 0 being the coarsest and 15 the finest. Many functions require a resolution parameter.
*   **Cells**: H3 partitions the globe into hexagonal cells. At each resolution, there are 12 pentagonal cells located at the vertices of the underlying icosahedron.
*   **Error Handling**: Functions that can fail return a `Result<T, H3Error>`, where `H3Error` is an enum detailing the type of error.

The main entry points for interacting with the library are the free-standing functions, often categorized by functionality (indexing, traversal, hierarchy, etc.).

## 2. Main Types and Their Public Methods

The library primarily uses structs to represent data. These structs generally do not have public methods for direct manipulation beyond what's provided by the free-standing functions. The core data types are:

### `xs_h3::H3Index`

A transparent wrapper around a `u64` representing an H3 cell or edge index.

*   **Public Fields**: `pub u64` (the raw index value)

### `xs_h3::LatLng`

Represents a geographic coordinate.

*   **Public Fields**:
    *   `lat: f64` (Latitude in radians)
    *   `lng: f64` (Longitude in radians)

### `xs_h3::CellBoundary`

Represents the boundary of an H3 cell.

*   **Public Fields**:
    *   `num_verts: usize` (Number of vertices in the boundary)
    *   `verts: [LatLng; MAX_CELL_BNDRY_VERTS]` (Array of `LatLng` vertices. `MAX_CELL_BNDRY_VERTS` is typically 10.)

### `xs_h3::GeoLoop`

Represents a single closed loop of geographic coordinates.

*   **Public Fields**:
    *   `num_verts: usize` (Number of vertices in the loop)
    *   `verts: Vec<LatLng>` (Vertices forming the loop)

### `xs_h3::GeoPolygon`

Represents a polygon with an outer loop and optional inner hole loops.

*   **Public Fields**:
    *   `geoloop: GeoLoop` (The outer loop of the polygon)
    *   `num_holes: usize` (Number of hole loops)
    *   `holes: Vec<GeoLoop>` (Array of hole loops)

### `xs_h3::BBox`

Geographic bounding box with coordinates defined in radians.

*   **Public Fields**:
    *   `north: f64` (North latitude in radians)
    *   `south: f64` (South latitude in radians)
    *   `east: f64` (East longitude in radians)
    *   `west: f64` (West longitude in radians)

### `xs_h3::CoordIJ`

IJ hexagon coordinates.

*   **Public Fields**:
    *   `i: i32`
    *   `j: i32`

### `xs_h3::CoordIJK`

IJK hexagon coordinates.

*   **Public Fields**:
    *   `i: i32`
    *   `j: i32`
    *   `k: i32`

### `xs_h3::FaceIJK`

Face number and IJK coordinates on that face-centered coordinate system.

*   **Public Fields**:
    *   `face: i32` (Icosahedron face number, 0-19)
    *   `coord: CoordIJK` (IJK coordinates on that face)

### `xs_h3::Vec2d`

2D floating-point vector.

*   **Public Fields**:
    *   `x: f64`
    *   `y: f64`

### `xs_h3::Vec3d`

3D floating-point vector.

*   **Public Fields**:
    *   `x: f64`
    *   `y: f64`
    *   `z: f64`

### `xs_h3::PolygonRust`
(from `xs_h3::regions::to_polygon`, re-exported via `xs_h3::regions`)
Represents a single polygon with an outer loop and holes, using `Vec<LatLng>`.

*   **Public Fields**:
    *   `outer: Vec<LatLng>`
    *   `holes: Vec<Vec<LatLng>>`

## 3. Public Enums (Non-Config)

### `xs_h3::Direction`
Represents H3 digit axes directions. `repr(u8)`.
*   `Center = 0`
*   `KAxes = 1`
*   `JAxes = 2`
*   `JkAxes = 3`
*   `IAxes = 4`
*   `IkAxes = 5`
*   `IjAxes = 6`
*   `InvalidDigit = 7`

### `xs_h3::ContainmentMode`
(from `xs_h3::types`, used in `xs_h3::regions`)
Values representing polyfill containment modes. `repr(u32)`.
*   `Center = 0`: Cell center is contained in the shape.
*   `Full = 1`: Cell is fully contained in the shape.
*   `Overlapping = 2`: Cell overlaps the shape at any point.
*   `OverlappingBbox = 3`: Cell bounding box overlaps shape.
*   `Invalid = 4`: Invalid mode, should not be used.

## 4. Public Functions (Free-standing)

All functions return `Result<T, H3Error>` if they can fail.

### Module `xs_h3::indexing`

*   `pub fn lat_lng_to_cell(geo: &LatLng, res: i32) -> Result<H3Index, H3Error>`
    *   Finds the H3 cell containing the given `LatLng` point at the specified resolution.
*   `pub fn cell_to_lat_lng(cell: H3Index) -> Result<LatLng, H3Error>`
    *   Finds the center `LatLng` point of the given H3 cell.
*   `pub fn cell_to_boundary(cell: H3Index) -> Result<CellBoundary, H3Error>`
    *   Finds the boundary of the given H3 cell.

### Module `xs_h3::h3_index::inspection` (Re-exported at `xs_h3::`)

*   `pub fn get_resolution(h: H3Index) -> i32`
    *   Gets the resolution of the H3 index. (Note: In `xs_h3::h3_index` internally, then re-exported)
*   `pub fn get_base_cell_number(h: H3Index) -> i32`
    *   Returns the base cell number (0-121) for an H3 cell index.
*   `pub fn is_valid_cell(h: H3Index) -> bool`
    *   Validates an H3 cell index.
*   `pub fn is_pentagon(h: H3Index) -> bool`
    *   Determines if an H3 cell is a pentagon.
*   `pub fn is_res_class_iii(h: H3Index) -> bool`
    *   Determines if an H3 cell's resolution is Class III (odd).
*   `pub fn get_num_cells(res: i32) -> Result<i64, H3Error>`
    *   Number of unique H3 cells at the given resolution.
*   `pub fn pentagon_count() -> i32`
    *   Number of H3 pentagons at any given resolution (always 12).
*   `pub fn get_pentagons(res: i32, out_pentagons: &mut [H3Index; 12]) -> Result<(), H3Error>`
    *   Get all H3 pentagon indexes at the specified resolution.
*   `pub fn get_res0_cells(out_res0_cells: &mut [H3Index; NUM_BASE_CELLS])`
    *   Get all 122 H3 resolution 0 indexes. (`NUM_BASE_CELLS` is typically 122)
*   `pub fn max_face_count(_h: H3Index) -> Result<usize, H3Error>`
    *   Maximum number of icosahedron faces an H3 cell's boundary may cross (typically 2).
*   `pub fn get_icosahedron_faces(h: H3Index, out_faces: &mut [i32]) -> Result<usize, H3Error>`
    *   Finds all icosahedron faces intersected by a given H3 cell. `out_faces` should be sized by `max_face_count`.

### Module `xs_h3::h3_index::string_conv` (Re-exported at `xs_h3::`)

*   `pub fn string_to_h3(s: &str) -> Result<H3Index, H3Error>`
    *   Converts a string representation (hexadecimal) of an H3 index into an `H3Index`.
*   `pub fn h3_to_string(h: H3Index, buffer: &mut [u8]) -> Result<(), H3Error>`
    *   Converts an `H3Index` into its string representation, writing to the provided buffer.
*   `pub fn h3_to_string_alloc(h: H3Index) -> String`
    *   Converts an `H3Index` into an allocated `String`.

### Module `xs_h3::hierarchy`

*   `pub fn cell_to_parent(h: H3Index, parent_res: i32) -> Result<H3Index, H3Error>`
    *   Produces the parent H3 index of `h` at `parent_res`.
*   `pub fn cell_to_children_size(h: H3Index, child_res: i32) -> Result<i64, H3Error>`
    *   Determines the exact number of children for a cell at `child_res`.
*   `pub fn cell_to_center_child(h: H3Index, child_res: i32) -> Result<H3Index, H3Error>`
    *   Returns the center child of `h` at `child_res`.
*   `pub fn cell_to_children(h: H3Index, child_res: i32, children: &mut [H3Index]) -> Result<(), H3Error>`
    *   Fills `children` with all H3 cells that are children of `h` at `child_res`.
*   `pub fn cell_to_child_pos(child: H3Index, parent_res: i32) -> Result<i64, H3Error>`
    *   Returns the position of `child` within an ordered list of all children of its parent at `parent_res`.
*   `pub fn child_pos_to_cell(child_pos: i64, parent: H3Index, child_res: i32) -> Result<H3Index, H3Error>`
    *   Returns the child cell at `child_pos` of `parent` at `child_res`.
*   `pub fn compact_cells(h3_set: &mut [H3Index], out_compacted_set: &mut [H3Index]) -> Result<usize, H3Error>`
    *   Compacts a set of H3 cells of the same resolution. `h3_set` may be modified (sorted).
*   `pub fn uncompact_cells_size(compacted_set: &[H3Index], res: i32) -> Result<i64, H3Error>`
    *   Calculates the number of cells resulting from uncompacting `compacted_set` to `res`.
*   `pub fn uncompact_cells(compacted_set: &[H3Index], res: i32, out_set: &mut [H3Index]) -> Result<(), H3Error>`
    *   Uncompacts `compacted_set` to `res`, filling `out_set`.

### Module `xs_h3::traversal`

*   `pub fn grid_distance(origin: H3Index, destination: H3Index) -> Result<i64, H3Error>`
    *   Produces the grid distance between two H3 cells.
*   `pub fn max_grid_disk_size(k: i32) -> Result<i64, H3Error>`
    *   Maximum number of H3 cells in a k-ring disk of radius `k`.
*   `pub fn grid_disk(origin: H3Index, k: i32, out_cells: &mut [H3Index]) -> Result<(), H3Error>`
    *   Produces H3 cells within `k` distance of `origin`. Handles pentagons.
*   `pub fn grid_disk_distances(origin: H3Index, k: i32, out_cells: &mut [H3Index], out_distances_opt: Option<&mut [i32]>) -> Result<(), H3Error>`
    *   Produces H3 cells and their distances from `origin`, up to distance `k`.
*   `pub fn grid_ring_unsafe(origin: H3Index, k: i32, out_cells: &mut [H3Index]) -> Result<(), H3Error>`
    *   Returns the "hollow" ring of cells at distance `k`. Behavior undefined if a pentagon is encountered.
*   `pub fn are_neighbor_cells(origin: H3Index, destination: H3Index) -> Result<bool, H3Error>`
    *   Returns whether or not the provided H3 cells are neighbors.
*   `pub fn grid_path_cells_size(start: H3Index, end: H3Index) -> Result<i64, H3Error>`
    *   Number of H3 cells in a line from `start` to `end`.
*   `pub fn grid_path_cells(start: H3Index, end: H3Index, out_path: &mut [H3Index]) -> Result<(), H3Error>`
    *   Given two H3 cells, returns the line of H3 cells between them (inclusive).

### Module `xs_h3::regions`

*   `pub fn max_polygon_to_cells_size(polygon: &GeoPolygon, res: i32, flags: u32) -> Result<i64, H3Error>`
    *   Estimates the maximum number of H3 cells needed to polyfill `polygon` at `res`.
*   `pub fn polygon_to_cells(polygon: &GeoPolygon, res: i32, flags: u32, out: &mut [H3Index]) -> Result<(), H3Error>`
    *   Fills a geopolitical `polygon` with H3 cells at resolution `res`. `flags` control containment mode.
*   `pub fn cells_to_multi_polygon(h3_set: &[H3Index]) -> Result<MultiPolygonRust, H3Error>`
    *   Converts a set of H3 cells to a `MultiPolygonRust` (list of polygons with outer/hole loops).

### Module `xs_h3::measures`

*   `pub fn cell_area_rads2(cell: H3Index) -> Result<f64, H3Error>`
    *   Area of H3 cell in radians^2.
*   `pub fn cell_area_km2(cell: H3Index) -> Result<f64, H3Error>`
    *   Area of H3 cell in kilometers^2.
*   `pub fn cell_area_m2(cell: H3Index) -> Result<f64, H3Error>`
    *   Area of H3 cell in meters^2.
*   `pub fn exact_edge_length_rads(edge: H3Index) -> Result<f64, H3Error>`
    *   Exact length of a specific H3 directed edge in radians. (Currently placeholder)
*   `pub fn exact_edge_length_km(edge: H3Index) -> Result<f64, H3Error>`
    *   Exact length of a specific H3 directed edge in kilometers. (Currently placeholder)
*   `pub fn exact_edge_length_m(edge: H3Index) -> Result<f64, H3Error>`
    *   Exact length of a specific H3 directed edge in meters. (Currently placeholder)

### Module `xs_h3::latlng` (Geometric utilities)

*   `pub fn degs_to_rads(degrees: f64) -> f64`
*   `pub fn rads_to_degs(radians: f64) -> f64`
*   `pub fn great_circle_distance_rads(a: &LatLng, b: &LatLng) -> f64`
*   `pub fn great_circle_distance_km(a: &LatLng, b: &LatLng) -> f64`
*   `pub fn great_circle_distance_m(a: &LatLng, b: &LatLng) -> f64`
*   `pub fn get_hexagon_area_avg_km2(res: i32) -> Result<f64, H3Error>`
*   `pub fn get_hexagon_area_avg_m2(res: i32) -> Result<f64, H3Error>`
*   `pub fn get_hexagon_edge_length_avg_km(res: i32) -> Result<f64, H3Error>`
*   `pub fn get_hexagon_edge_length_avg_m(res: i32) -> Result<f64, H3Error>`

## 5. Public Type Aliases

### `xs_h3::regions::MultiPolygonRust`
`pub type MultiPolygonRust = Vec<PolygonRust>;`
    *   Represents a collection of polygons, where each `PolygonRust` has an outer boundary and a list of hole boundaries.

## 6. Public Constants

### `xs_h3::H3_NULL`
`pub const H3_NULL: H3Index = H3Index(0);`
    *   Represents an invalid or null H3 index.

### `xs_h3::MAX_CELL_BNDRY_VERTS`
`pub const MAX_CELL_BNDRY_VERTS: usize = 10;`
    *   Maximum number of vertices a cell boundary can have (pentagons can have up to 10 after distortion).

## 7. Error Handling

The primary error type is `xs_h3::H3Error`. Most functions that can fail return `Result<T, H3Error>`.

### Enum `xs_h3::H3Error`
`repr(u32)`
*   `Success = 0`: Success (no error).
*   `Failed = 1`: The operation failed but a more specific error is not available.
*   `Domain = 2`: Argument was outside of acceptable range.
*   `LatLngDomain = 3`: Latitude or longitude arguments were outside of acceptable range.
*   `ResDomain = 4`: Resolution argument was outside of acceptable range.
*   `CellInvalid = 5`: `H3Index` cell argument was not valid.
*   `DirEdgeInvalid = 6`: `H3Index` directed edge argument was not valid.
*   `UndirEdgeInvalid = 7`: `H3Index` undirected edge argument was not valid.
*   `VertexInvalid = 8`: `H3Index` vertex argument was not valid.
*   `Pentagon = 9`: Pentagon distortion was encountered which the algorithm could not handle.
*   `DuplicateInput = 10`: Duplicate input was encountered in the arguments and the algorithm could not handle it.
*   `NotNeighbors = 11`: `H3Index` cell arguments were not neighbors.
*   `ResMismatch = 12`: `H3Index` cell arguments had incompatible resolutions.
*   `MemoryAlloc = 13`: Necessary memory allocation failed.
*   `MemoryBounds = 14`: Bounds of provided memory were not large enough.
*   `OptionInvalid = 15`: Mode or flags argument was not valid.