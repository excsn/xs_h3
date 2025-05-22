//! Core H3 data structures.

use crate::constants::MAX_CELL_BNDRY_VERTS;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
#[cfg(feature = "serde")]
use serde_repr::{Deserialize_repr, Serialize_repr};

/// Represents an H3 cell index or a directed H3 edge index.
/// This is a 64-bit unsigned integer.
#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct H3Index(pub u64);

/// Invalid H3 index, often used to signify an error or missing data.
pub const H3_NULL: H3Index = H3Index(0);

/// Latitude/longitude coordinates in radians.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct LatLng {
  /// Latitude in radians.
  pub lat: f64,
  /// Longitude in radians.
  pub lng: f64,
}

/// Represents the boundary of an H3 cell.
///
/// Contains the number of vertices and an array of `LatLng` coordinates
/// forming the cell boundary in counter-clockwise order.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CellBoundary {
  /// Number of vertices in the boundary.
  pub num_verts: usize,
  /// Array of vertices. The maximum number of vertices is defined by `MAX_CELL_BNDRY_VERTS`.
  /// Unused vertices at the end of the array are not significant.
  pub verts: [LatLng; MAX_CELL_BNDRY_VERTS],
}

impl Default for CellBoundary {
  fn default() -> Self {
    Self {
      num_verts: 0,
      verts: [LatLng::default(); MAX_CELL_BNDRY_VERTS],
    }
  }
}

/// Represents a single closed loop of geographic coordinates.
/// The last vertex is not implicitly connected to the first.
#[derive(Debug, Clone, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GeoLoop {
  /// Number of vertices in the loop.
  pub num_verts: usize,
  /// Vertices forming the loop.
  pub verts: Vec<LatLng>,
}

/// Represents a polygon with an outer loop and zero or more inner hole loops.
#[derive(Debug, Clone, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GeoPolygon {
  /// The outer loop of the polygon.
  pub geoloop: GeoLoop,
  /// Number of hole loops.
  pub num_holes: usize,
  /// Array of hole loops.
  pub holes: Vec<GeoLoop>,
}

// TODO: GeoMultiPolygon if needed later

/// Represents an H3 error code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)]
#[cfg_attr(feature = "serde", derive(Serialize_repr, Deserialize_repr))]
#[allow(clippy::enum_variant_names)] // To match C naming
pub enum H3Error {
  /// Success (no error).
  Success = 0,
  /// The operation failed but a more specific error is not available.
  Failed = 1,
  /// Argument was outside of acceptable range.
  Domain = 2,
  /// Latitude or longitude arguments were outside of acceptable range.
  LatLngDomain = 3,
  /// Resolution argument was outside of acceptable range.
  ResDomain = 4,
  /// `H3Index` cell argument was not valid.
  CellInvalid = 5,
  /// `H3Index` directed edge argument was not valid.
  DirEdgeInvalid = 6,
  /// `H3Index` undirected edge argument was not valid.
  UndirEdgeInvalid = 7,
  /// `H3Index` vertex argument was not valid.
  VertexInvalid = 8,
  /// Pentagon distortion was encountered which the algorithm could not handle.
  Pentagon = 9,
  /// Duplicate input was encountered in the arguments and the algorithm could not handle it.
  DuplicateInput = 10,
  /// `H3Index` cell arguments were not neighbors.
  NotNeighbors = 11,
  /// `H3Index` cell arguments had incompatible resolutions.
  ResMismatch = 12,
  /// Necessary memory allocation failed.
  MemoryAlloc = 13,
  /// Bounds of provided memory were not large enough.
  MemoryBounds = 14,
  /// Mode or flags argument was not valid.
  OptionInvalid = 15,
}

/// IJ hexagon coordinates. Each axis is spaced 120 degrees apart.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CoordIJ {
  /// I component.
  pub i: i32,
  /// J component.
  pub j: i32,
}

/// IJK hexagon coordinates. Each axis is spaced 120 degrees apart.
/// The K component is derived from `i` and `j` (`k = -i - j`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct CoordIJK {
  /// I component.
  pub i: i32,
  /// J component.
  pub j: i32,
  /// K component.
  pub k: i32,
}

/// Face number and IJK coordinates on that face-centered coordinate system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct FaceIJK {
  /// Icosahedron face number (0-19).
  pub face: i32,
  /// IJK coordinates on that face.
  pub coord: CoordIJK,
}

/// 2D floating-point vector.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Vec2d {
  /// X component.
  pub x: f64,
  /// Y component.
  pub y: f64,
}

/// 3D floating-point vector.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Vec3d {
  /// X component.
  pub x: f64,
  /// Y component.
  pub y: f64,
  /// Z component.
  pub z: f64,
}

/// H3 digit representing IJK+ axes direction (0-6), or invalid (7).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Hash)]
#[repr(u8)]
#[cfg_attr(feature = "serde", derive(Serialize_repr, Deserialize_repr))]
pub enum Direction {
  /// H3 digit in center.
  Center = 0,
  /// H3 digit in k-axes direction.
  KAxes = 1,
  /// H3 digit in j-axes direction.
  JAxes = 2,
  /// H3 digit in j == k direction.
  JkAxes = 3, // J_AXES_DIGIT | K_AXES_DIGIT
  /// H3 digit in i-axes direction.
  IAxes = 4,
  /// H3 digit in i == k direction.
  IkAxes = 5, // I_AXES_DIGIT | K_AXES_DIGIT
  /// H3 digit in i == j direction.
  IjAxes = 6, // I_AXES_DIGIT | J_AXES_DIGIT
  /// H3 digit in the invalid direction.
  InvalidDigit = 7,
  // NUM_DIGITS = INVALID_DIGIT (from C)
  // PENTAGON_SKIPPED_DIGIT = K_AXES_DIGIT (from C)
}

impl Default for Direction {
  fn default() -> Self {
    Direction::Center
  }
}

impl TryFrom<u8> for Direction {
  type Error = H3Error; // Or a more specific error type

  fn try_from(value: u8) -> Result<Self, Self::Error> {
    match value {
      0 => Ok(Direction::Center),
      1 => Ok(Direction::KAxes),
      2 => Ok(Direction::JAxes),
      3 => Ok(Direction::JkAxes),
      4 => Ok(Direction::IAxes),
      5 => Ok(Direction::IkAxes),
      6 => Ok(Direction::IjAxes),
      7 => Ok(Direction::InvalidDigit), // Technically valid enum variant, but invalid as H3 digit
      _ => Err(H3Error::Domain),        // Value out of range for Direction
    }
  }
}

/// Geographic bounding box with coordinates defined in radians.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct BBox {
  /// North latitude in radians.
  pub north: f64,
  /// South latitude in radians.
  pub south: f64,
  /// East longitude in radians.
  pub east: f64,
  /// West longitude in radians.
  pub west: f64,
}

/// Values representing polyfill containment modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u32)] // To match C flags if needed, though Rust usually uses u8 or usize for enums like this
#[cfg_attr(feature = "serde", derive(Serialize_repr, Deserialize_repr))]
pub enum ContainmentMode {
  /// Cell center is contained in the shape.
  Center = 0,
  /// Cell is fully contained in the shape.
  Full = 1,
  /// Cell overlaps the shape at any point.
  Overlapping = 2,
  /// Cell bounding box overlaps shape.
  OverlappingBbox = 3,
  /// This mode is invalid and should not be used (used for bounds checking).
  Invalid = 4,
}

impl Default for ContainmentMode {
  fn default() -> Self {
    ContainmentMode::Center // Or whatever default makes sense
  }
}

// Optional: if you need to convert from u32 flags easily
impl TryFrom<u32> for ContainmentMode {
  type Error = H3Error;
  fn try_from(value: u32) -> Result<Self, Self::Error> {
    match value {
      0 => Ok(ContainmentMode::Center),
      1 => Ok(ContainmentMode::Full),
      2 => Ok(ContainmentMode::Overlapping),
      3 => Ok(ContainmentMode::OverlappingBbox),
      _ => Err(H3Error::OptionInvalid), // Or H3Error::Domain
    }
  }
}
