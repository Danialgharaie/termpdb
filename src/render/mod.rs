//! 3D Software Rendering Pipeline for TermPDB.
//!
//! Provides truecolor half-block framebuffer rendering, floating-point Z-buffering,
//! orbit camera controls, directional lighting with depth fog, color schemes,
//! and software rasterizers for spheres, cylinders, lines, and triangles.

pub mod braille;
pub mod buffer;
pub mod camera;
pub mod color;
pub mod export;
pub mod lighting;
pub mod postprocess;
pub mod rasterizer;
pub mod representations;

pub use braille::BrailleBuffer;
pub use buffer::{Framebuffer, FramebufferBand, PixelColor};
pub use camera::Camera;
pub use color::{ColorScheme, color_for_atom};
pub use export::{downsample_rgba, export_mp4, render_supersampled, render_svg, write_png};
pub use lighting::Lighting;
pub use postprocess::{PostProcessConfig, apply_postprocessing};
pub use rasterizer::{
    draw_cylinder, draw_dashed_line_3d, draw_line_3d, draw_overlay_line, draw_overlay_sphere,
    draw_sphere, draw_sphere_band, draw_triangle_3d,
};
pub use representations::{
    LodMode, RenderContext, RenderMode, RibbonPrimitive, Visibility, build_render_cache,
    build_ribbon_geometry, render_structure, render_structure_ctx,
};

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
    export_ansi_with_visibility(
        structure,
        mode,
        color_scheme,
        width,
        height,
        Visibility::default(),
        LodMode::Auto,
    )
}

/// Headless ANSI export with an explicit atom-visibility filter.
#[allow(clippy::too_many_arguments)]
pub fn export_ansi_with_visibility(
    structure: &Structure,
    mode: RenderMode,
    color_scheme: ColorScheme,
    width: u16,
    height: u16,
    visibility: Visibility,
    lod: LodMode,
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
        visibility,
        lod,
    );
    apply_postprocessing(&mut buffer, &PostProcessConfig::default());
    buffer.to_ansi()
}

/// Draws selected atoms (and a distance line) on top of the current framebuffer.
pub fn draw_selection_markers(
    structure: &Structure,
    camera: &Camera,
    buffer: &mut Framebuffer,
    indices: &[usize],
) {
    if indices.is_empty() {
        return;
    }

    let colors = [(255, 220, 40), (80, 220, 255)];
    let mut screens: Vec<(f32, f32, f32)> = Vec::new();

    for (n, &idx) in indices.iter().enumerate() {
        let Some(atom) = structure.atoms().get(idx) else {
            continue;
        };
        let Some(pt) = camera.world_to_screen(atom.pos, buffer.width, buffer.height) else {
            continue;
        };
        let r = representations::project_radius(0.55, pt.2, camera.fov, buffer.height).max(2.5);
        draw_overlay_sphere(buffer, pt, r, colors[n.min(1)]);
        screens.push(pt);
    }

    if screens.len() == 2 {
        draw_overlay_line(buffer, screens[0], screens[1], (255, 255, 255));
    }
}
