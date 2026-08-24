//! All-atom Ball and Stick representation.
//!
//! Renders every atom as an analytical sphere (R ~ 0.35 A) and every covalent
//! bond as a thick shaded cylinder split proportionally by atom colors.

use crate::model::Bond;
use crate::model::bond::BondDetector;
use crate::render::buffer::Framebuffer;
use crate::render::rasterizer::{draw_cylinder, draw_sphere};
use crate::render::representations::{LodLevel, RenderContext, project_radius};

/// Renders the molecular structure in Ball and Stick mode.
pub fn render_ball_stick(ctx: &RenderContext, buffer: &mut Framebuffer) {
    let structure = ctx.structure;
    let camera = ctx.camera;
    let mats = ctx.mats;
    let lighting = ctx.lighting;
    let lod = ctx.lod;
    let colors = ctx.colors;
    let visible = ctx.visible;

    // 1. Render all atoms as analytical spheres
    let atoms = structure.atoms();
    let mut projected: Vec<((f32, f32, f32), f32, crate::render::buffer::PixelColor)> =
        Vec::with_capacity(atoms.len());

    for atom in atoms {
        if !visible[atom.index] {
            continue;
        }
        // project_sphere keeps spheres whose center is behind the near plane
        // but whose surface still crosses it (conservative approximation).
        if let Some((pt, sphere_r)) =
            camera.project_sphere(&mats, atom.pos, 0.38, buffer.width, buffer.height)
        {
            if sphere_r < 0.4 {
                continue;
            }
            let sphere_r = if lod == LodLevel::Full {
                sphere_r.max(0.8)
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

    // 2. Render all covalent bonds as cylinders
    let detected_bonds;
    let bonds: &[Bond] = if !structure.bonds().is_empty() {
        structure.bonds()
    } else {
        detected_bonds = BondDetector::detect_bonds(structure.atoms());
        &detected_bonds
    };

    for bond in bonds {
        if bond.atom1_idx < structure.atoms().len() && bond.atom2_idx < structure.atoms().len() {
            let atom1 = &structure.atoms()[bond.atom1_idx];
            let atom2 = &structure.atoms()[bond.atom2_idx];
            if !visible[atom1.index] || !visible[atom2.index] {
                continue;
            }

            // Segment clip against the near plane keeps half-behind bonds drawn.
            if let Some((p1, p2)) =
                camera.project_segment(&mats, atom1.pos, atom2.pos, buffer.width, buffer.height)
            {
                let avg_depth = (p1.2 + p2.2) * 0.5;
                let bond_r = project_radius(0.18, avg_depth, camera.fov, buffer.height).max(0.5);

                let c1 = colors[atom1.index];
                let c2 = colors[atom2.index];

                if c1 == c2 {
                    draw_cylinder(buffer, p1, p2, bond_r, c1, lighting);
                } else {
                    let pmid = ((p1.0 + p2.0) * 0.5, (p1.1 + p2.1) * 0.5, avg_depth);
                    draw_cylinder(buffer, p1, pmid, bond_r, c1, lighting);
                    draw_cylinder(buffer, pmid, p2, bond_r, c2, lighting);
                }
            }
        }
    }
}
