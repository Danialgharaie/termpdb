//! 3D Software Rendering Pipeline for TermPDB.
//!
//! Provides truecolor half-block framebuffer rendering, floating-point Z-buffering,
//! orbit camera controls, directional lighting with depth fog, color schemes,
//! and software rasterizers for spheres, cylinders, lines, and triangles.

pub mod buffer;
pub mod camera;
pub mod color;
pub mod lighting;
pub mod rasterizer;
pub mod representations;

pub use buffer::{Framebuffer, PixelColor};
pub use camera::Camera;
pub use color::{ColorScheme, color_for_atom};
pub use lighting::Lighting;
pub use rasterizer::{draw_cylinder, draw_line_3d, draw_sphere, draw_triangle_3d};
pub use representations::{RenderMode, render_structure};
