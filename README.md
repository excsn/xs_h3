# XS-H3: A High-Performance H3 Geospatial Indexing Library

[![Crates.io](https://img.shields.io/crates/v/xs_h3.svg)](https://crates.io/crates/xs_h3)
[![Docs.rs](https://docs.rs/xs_h3/badge.svg)](https://docs.rs/xs_h3)
[![License: MPL-2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)

XS-H3 is Excerion Sun's pure Rust implementation of Uber's H3, a discrete global grid system. It provides a set of functions for converting between geographic coordinates (latitude/longitude) and H3 cell indexes, as well as for performing grid traversal, hierarchical operations, and region-based analysis.

This library is designed with a primary focus on **high performance** and **memory safety**, leveraging the strengths of Rust. It aims to serve as a robust and efficient alternative or complement to existing H3 libraries by closely following the C H3 library's API and behavior where appropriate.

The "xs" in `xs-h3` can be thought of as "Cross-Safe" or "Extra Safe," highlighting the emphasis on Rust's safety guarantees.

---

**⚠️ Experimental Status ⚠️**

**`xs-h3` (version 0.1.0) is currently in an early, experimental stage of development.**

*   **API Instability:** The public API is subject to change without notice as the library evolves and matures.
*   **Incomplete Features:** While many core H3 functionalities are present, not all features from the C H3 library are fully implemented or optimized. Some existing features may be partial or have known limitations.
*   **Limited Testing:** Comprehensive testing across all edge cases, resolutions, and direct behavioral comparisons with the C H3 library is ongoing.
*   **Potential Bugs:** As with any early-stage software, there may be bugs or correctness issues.

**Use in production environments is NOT recommended at this stage.** We encourage experimentation, feedback, and contributions to help stabilize and complete the library.

---

## What is H3?

H3 is a geospatial indexing system that partitions the world into hexagonal cells across multiple resolutions. This hexagonal grid offers several advantages for spatial analysis, including more uniform adjacency and distance calculations compared to traditional square grids.

Key features of the H3 system include:

*   **Hierarchical Indexing**: Cells at finer resolutions are perfectly contained within coarser resolution cells.
*   **Global Coverage**: Covers the entire globe using an icosahedron-based projection.
*   **Hexagonal Grid**: Provides approximately uniform cell sizes and shapes at each resolution (with 12 pentagonal cells per resolution at icosahedron vertices).
*   **Efficient Algorithms**: For neighborhood lookup, region filling, and other common geospatial operations.

## XS-H3 Features (Current & Planned)

XS-H3 aims to provide a comprehensive implementation of the H3 specification. Current capabilities include:

*   **Indexing**:
    *   Conversion from latitude/longitude to H3 cell index (`lat_lng_to_cell`). (Largely implemented)
    *   Conversion from H3 cell index to its center latitude/longitude (`cell_to_lat_lng`). (Implemented)
    *   Conversion from H3 cell index to its boundary coordinates (`cell_to_boundary`). (Implemented)
*   **Inspection**:
    *   Determining the resolution, base cell, or high bit of an H3 index. (Implemented)
    *   Validating H3 indexes (`is_valid_cell`). (Implemented)
    *   Checking if a cell is a pentagon (`is_pentagon`) or Class III (`is_res_class_iii`). (Implemented)
    *   Getting lists of all resolution 0 cells or all pentagons at a resolution. (Implemented)
*   **Traversal**:
    *   Finding cells within a k-ring distance (`grid_disk`, `grid_disk_distances`). (Safe BFS version implemented)
    *   Finding cells forming a "hollow" ring (`grid_ring_unsafe`). (Implemented, unsafe due to pentagons)
    *   Finding cells forming a line between two H3 indexes (`grid_path_cells`). (Implemented, some path deviations from C noted)
    *   Calculating grid distance between two cells (`grid_distance`). (Implemented, some distance deviations from C noted for specific pairs due to underlying local IJK differences)
    *   Identifying neighboring cells (`are_neighbor_cells`). (Implemented)
*   **Hierarchy**:
    *   Finding the parent or children of an H3 cell (`cell_to_parent`, `cell_to_children`, `cell_to_center_child`). (Implemented)
    *   Determining child position relative to parent (`cell_to_child_pos`, `child_pos_to_cell`). (Implemented)
    *   Compacting a set of H3 cells into a minimal representation (`compact_cells`). (Basic implementation)
    *   Uncompacting a set of cells to a specific resolution (`uncompact_cells`). (Implemented)
*   **Regions (Polyfill)**:
    *   Filling a geographic polygon with H3 cells (`polygon_to_cells`). (Experimental, basic implementation)
    *   Estimating the maximum number of cells needed for a polyfill (`max_polygon_to_cells_size`). (Experimental)
    *   Converting a set of H3 cells into a multi-polygon representation (`cells_to_multi_polygon`). (Experimental, hole handling and full normalization are basic)
*   **Measurements & Utilities**:
    *   Average cell area and edge length calculations. (Implemented)
    *   Exact cell area calculations (`cell_area_km2`, etc.). (Implemented)
    *   Conversion between H3 index and string representations. (Implemented)
    *   Optional Serde support via the `serde` feature. (Implemented for core types)

**Planned / In Progress for future versions:**
*   Full Local IJ Coordinate System support and related APIs.
*   Complete Directed Edge functions (creation, inspection, traversal).
*   Vertex functions.
*   Robust and complete `cells_to_multi_polygon` (including full GeoJSON-compliant hole handling and winding order normalization).
*   Comprehensive benchmarking against C H3 and performance optimization.
*   More extensive test coverage, including property-based testing and systematic comparison with C H3 outputs.
*   Detailed error handling and potentially more specific `H3Error` variants.
*   Exact edge length functions for directed edges.

## Performance Goals

The primary engineering goal of XS-H3 is to achieve **high performance**, while providing Rust's strong memory safety guarantees. This is pursued through:

*   Careful, idiomatic Rust algorithm implementation.
*   Leveraging Rust's zero-cost abstractions and control over memory layout.
*   Focus on minimizing allocations and optimizing critical code paths.

Benchmarks against reference implementations are planned as a key part of the development and validation process.

## Use Cases

H3 (and therefore XS-H3) is suitable for a wide range of geospatial applications, including:

*   **Data Aggregation**: Grouping point data into discrete grid cells for analysis and visualization.
*   **Spatial Indexing**: Efficiently querying and joining spatial datasets.
*   **Proximity Analysis**: Finding nearby points of interest or features.
*   **Coverage Analysis**: Determining the H3 cells that cover a given area.
*   **Large-Scale Geospatial Processing**: Handling and analyzing massive datasets.

## Getting Started

Add `xs-h3` to your `Cargo.toml`:

```toml
[dependencies]
xs-h3 = "0.1.0"
```

If you need Serde support for `H3Index`, `LatLng`, etc.:
```toml
[dependencies]
xs-h3 = { version = "0.1.0", features = ["serde"] }
```

**Given the experimental nature of version 0.1.0, please lock to this specific version if you decide to try it, as future 0.1.x versions may introduce breaking API changes.**

### Basic Usage Example

```rust
use xs_h3::{lat_lng_to_cell, cell_to_lat_lng, get_resolution, H3Index, LatLng, degs_to_rads, rads_to_degs, H3Error};

fn main() -> Result<(), H3Error> {
    let sf_lat_deg = 37.7749;
    let sf_lng_deg = -122.4194;

    let geo_point = LatLng {
        lat: degs_to_rads(sf_lat_deg),
        lng: degs_to_rads(sf_lng_deg),
    };

    let resolution = 9;
    let h3_cell: H3Index = lat_lng_to_cell(&geo_point, resolution)?;
    println!("H3 Cell for ({}, {}) at res {}: {:x}", sf_lat_deg, sf_lng_deg, resolution, h3_cell.0);

    let cell_res = get_resolution(h3_cell);
    println!("Resolution of cell {:x}: {}", h3_cell.0, cell_res);

    let center_coords_rad: LatLng = cell_to_lat_lng(h3_cell)?;
    println!(
        "Center of cell {:x}: (lat: {:.6}, lng: {:.6}) degrees",
        h3_cell.0,
        rads_to_degs(center_coords_rad.lat),
        rads_to_degs(center_coords_rad.lng)
    );
    Ok(())
}
```
*(For more examples, please see the `examples/`, `tests/` and `benches/` directories in the repository.)*

## Contributing

Contributions are highly welcome, especially given the library's experimental state! Please feel free to open an issue to discuss potential changes, report bugs, or submit a pull request for fixes or new features.
When contributing, please try to:
*   Add tests for new functionality or bug fixes.
*   Follow Rust's idiomatic coding style (`cargo fmt` and `cargo clippy`).
*   Benchmark performance-sensitive changes.
*   Clearly document any deviations from the H3 C library's behavior if they are intentional, or note areas where parity is still being worked towards.

## License

XS-H3 is licensed under the Mozilla Public License Version 2.0 (MPL-2.0). See the [LICENSE](LICENSE) for details.

## Acknowledgements

This library is an implementation of the H3 global grid system, originally developed by Uber Technologies, Inc. We are grateful for their innovative and great pioneering work and for making the H3 specification and reference implementation publicly available for free.