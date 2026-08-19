//! Van der Waals (VDW) space-filling representation.
//!
//! Renders every atom in the molecular structure as a full-size analytical sphere
//! with radius corresponding to its element Van der Waals radius (1.2 - 2.2 Å).

use crate::model::Structure;
use crate::render::buffer::Framebuffer;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, color_for_atom};
use crate::render::lighting::Lighting;
use crate::render::rasterizer::draw_sphere;
use crate::render::representations::{build_atom_residue_map, project_radius};

/// Renders the molecular structure in VDW space-filling mode.
pub fn render_vdw(
    structure: &Structure,
    color_scheme: ColorScheme,
    camera: &Camera,
    buffer: &mut Framebuffer,
    lighting: &Lighting,
) {
    let residue_map = build_atom_residue_map(structure);

    for atom in &structure.atoms {
        if let Some(pt) = camera.world_to_screen(atom.pos, buffer.width, buffer.height) {
            let r_world = atom.vdw_radius();
            let r_world = if r_world > 0.1 { r_world } else { 1.5 };
            let sphere_r = project_radius(r_world, pt.2, camera.fov, buffer.height).max(1.0);

            let res = residue_map.get(&atom.index).copied();
            let color = color_for_atom(atom, res, structure, color_scheme);
            draw_sphere(buffer, pt, sphere_r, color, lighting);
        }
    }
}
