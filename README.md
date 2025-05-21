# XS-H3: A High-Performance H3 Geospatial Indexing Library

**EXPERIMENTAL LIBRARY**

XS-H3 is a [Rust-based] implementation of Uber's H3, a discrete global grid system. It provides a set of functions for converting between geographic coordinates (latitude/longitude) and H3 cell indexes, as well as for performing grid traversal, hierarchical operations, and region-based analysis.

This library is designed with a primary focus on **high performance** and **memory safety**, leveraging the strengths of its underlying implementation language [Rust]. It aims to serve as a robust and efficient alternative or complement to existing H3 libraries.

## What is H3?

H3 is a geospatial indexing system that partitions the world into hexagonal cells across multiple resolutions. This hexagonal grid offers several advantages for spatial analysis, including more uniform adjacency and distance calculations compared to traditional square grids.

Key features of the H3 system include:

*   **Hierarchical Indexing**: Cells at finer resolutions are perfectly contained within coarser resolution cells.
*   **Global Coverage**: Covers the entire globe using an icosahedron-based projection.
*   **Hexagonal Grid**: Provides approximately uniform cell sizes and shapes at each resolution (with 12 pentagonal cells per resolution at icosahedron vertices).
*   **Efficient Algorithms**: For neighborhood lookup, region filling, and other common geospatial operations.

## XS-H3 Features

XS-H3 aims to provide a comprehensive implementation of the H3 specification, including:

*   **Indexing**:
    *   Conversion from latitude/longitude to H3 cell index at a given resolution (`latLngToCell`).
    *   Conversion from H3 cell index to its center latitude/longitude (`cellToLatLng`).
    *   Conversion from H3 cell index to its boundary coordinates (`cellToBoundary`).
*   **Inspection**:
    *   Determining the resolution or base cell of an H3 index.
    *   Validating H3 indexes (`isValidCell`).
    *   Checking if a cell is a pentagon (`isPentagon`).
*   **Traversal**:
    *   Finding cells within a k-ring distance (`gridDisk`, `gridDiskDistances`).
    *   Finding cells forming a line between two H3 indexes (`gridPathCells`).
    *   Calculating grid distance between two cells (`gridDistance`).
    *   Identifying neighboring cells (`areNeighborCells`).
*   **Hierarchy**:
    *   Finding the parent or children of an H3 cell (`cellToParent`, `cellToChildren`).
    *   Compacting a set of H3 cells into a minimal representation (`compactCells`).
    *   Uncompacting a set of cells to a specific resolution (`uncompactCells`).
*   **Regions**:
    *   Filling a geographic polygon with H3 cells (`polygonToCells`).
    *   Converting a set of H3 cells into a multi-polygon representation (`cellsToMultiPolygon`).
*   **Miscellaneous Utilities**:
    *   Average cell area and edge length calculations.
    *   Conversion between H3 index and string representations.

## Performance Goals

The primary engineering goal of XS-H3 is to achieve **high performance**, rivaling or exceeding that of the original C implementation where possible, while providing strong memory safety guarantees. This is achieved through:

*   Careful algorithm implementation.
*   Leveraging [Rust's] zero-cost abstractions and control over memory layout.
*   Focus on minimizing allocations and optimizing critical code paths.

Benchmarks against reference implementations are a key part of the development and validation process.

## Use Cases

H3 (and therefore XS-H3) is suitable for a wide range of geospatial applications, including:

*   **Data Aggregation**: Grouping point data into discrete grid cells for analysis and visualization.
*   **Spatial Indexing**: Efficiently querying and joining spatial datasets.
*   **Proximity Analysis**: Finding nearby points of interest or features.
*   **Coverage Analysis**: Determining the H3 cells that cover a given area.
*   **Large-Scale Geospatial Processing**: Handling and analyzing massive datasets.

## License

XS-H3 is licensed under the MPL 2.0. See `LICENSE` file for details.

## Acknowledgements

This library is an implementation of the H3 global grid system, originally developed by Uber Technologies, Inc. We are grateful for their pioneering work and for making the H3 specification and reference implementation publicly available.