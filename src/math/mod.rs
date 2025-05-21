// src/math/mod.rs
#![allow(clippy::module_name_repetitions)] // Common in math modules

pub mod extensions;
pub mod vec2d;
pub mod vec3d;

// Re-export if needed for internal crate use, or they can be used as `math::vec2d::foo`
// pub use extensions::*;
// pub use vec2d::*;
// pub use vec3d::*;