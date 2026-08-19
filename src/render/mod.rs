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

use crate::model::Structure;

/// Renders a molecular structure headlessly into a standalone ANSI truecolor string.
///
/// `width` specifies terminal columns and `height` specifies terminal rows.
pub fn export_ansi(
    structure: &Structure,
    mode: RenderMode,
    color_scheme: ColorScheme,
    width: u16,
    height: u16,
) -> String {
    if width == 0 || height == 0 {
        return String::new();
    }

    let pixel_width = width as usize;
    let pixel_height = (height as usize) * 2;
    let mut buffer = Framebuffer::new(pixel_width, pixel_height);
    let mut camera = Camera::new();
    camera.aspect = (pixel_width as f32) / (pixel_height as f32);
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);
    let lighting = Lighting::default();

    render_structure(
        structure,
        mode,
        color_scheme,
        &camera,
        &mut buffer,
        &lighting,
    );
    buffer.to_ansi()
}
