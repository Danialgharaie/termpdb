//! Van der Waals (VDW) space-filling representation.
//!
//! Renders every atom in the molecular structure as a full-size analytical sphere
//! with radius corresponding to its element Van der Waals radius (1.2 - 2.2 A).

use crate::render::buffer::Framebuffer;
use crate::render::rasterizer::draw_sphere;
use crate::render::representations::{LodLevel, RenderContext, project_radius};

/// Renders the molecular structure in VDW space-filling mode.
pub fn render_vdw(ctx: &RenderContext, buffer: &mut Framebuffer) {
    let structure = ctx.structure;
    let camera = ctx.camera;
    let mats = ctx.mats;
    let lighting = ctx.lighting;
    let lod = ctx.lod;
    let colors = ctx.colors;
    let visible = ctx.visible;

    // Coarse whole-scene cull: if even the largest atom projects to sub-pixel
    // size at the structure's nearest depth, nothing can be drawn -- skip the
    // per-atom loop entirely. (Safe for VDW: everything drawn is an atom sphere.)
    if let Some((_, _, depth)) = camera.project(&mats, ctx.com, buffer.width, buffer.height) {
        let near_depth = (depth - ctx.radius).max(0.1);
        if project_radius(ctx.max_vdw, near_depth, camera.fov, buffer.height) < 0.4 {
            return;
        }
    }

    let atoms = structure.atoms();
    let mut projected: Vec<((f32, f32, f32), f32, crate::render::buffer::PixelColor)> =
        Vec::with_capacity(atoms.len());

    for atom in atoms {
        if !visible[atom.index] {
            continue;
        }
        if let Some(pt) = camera.project(&mats, atom.pos, buffer.width, buffer.height) {
            let r_world = atom.vdw_radius();
            let r_world = if r_world > 0.1 { r_world } else { 1.5 };
            let sphere_r = project_radius(r_world, pt.2, camera.fov, buffer.height);
            if sphere_r < 0.4 {
                continue;
            }
            let sphere_r = if lod == LodLevel::Full {
                sphere_r.max(1.0)
            } else {
                sphere_r
            };

            let color = colors[atom.index];
            projected.push((pt, sphere_r, color));
        }
    }

    let band_height = 16;
    let mut bands = buffer.par_bands_mut(band_height);
    if bands.len() > 1 && projected.len() > 50 {
        use rayon::prelude::*;
        bands.par_iter_mut().for_each(|band| {
            for &(pt, r, color) in &projected {
                crate::render::rasterizer::draw_sphere_band(band, pt, r, color, lighting);
            }
        });
    } else {
        for (pt, r, color) in projected {
            draw_sphere(buffer, pt, r, color, lighting);
        }
    }
}
