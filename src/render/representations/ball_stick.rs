//! All-atom Ball and Stick representation.
//!
//! Renders every atom as an analytical sphere ($R \approx 0.35\text{ \AA}$) and every
//! covalent bond as a thick shaded cylinder split proportionally by atom colors.

use crate::model::bond::BondDetector;
use crate::model::{Bond, Structure};
use crate::render::buffer::Framebuffer;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, color_for_atom};
use crate::render::lighting::Lighting;
use crate::render::rasterizer::{draw_cylinder, draw_sphere};
use crate::render::representations::{build_atom_residue_map, project_radius};

/// Renders the molecular structure in Ball and Stick mode.
pub fn render_ball_stick(
    structure: &Structure,
    color_scheme: ColorScheme,
    camera: &Camera,
    buffer: &mut Framebuffer,
    lighting: &Lighting,
) {
    let residue_map = build_atom_residue_map(structure);

    // 1. Render all atoms as analytical spheres
    for atom in &structure.atoms {
        if let Some(pt) = camera.world_to_screen(atom.pos, buffer.width, buffer.height) {
            let r_world = 0.38;
            let sphere_r = project_radius(r_world, pt.2, camera.fov, buffer.height).max(0.8);
            let res = residue_map.get(&atom.index).copied();
            let color = color_for_atom(atom, res, structure, color_scheme);
            draw_sphere(buffer, pt, sphere_r, color, lighting);
        }
    }

    // 2. Render all covalent bonds as cylinders
    let detected_bonds;
    let bonds: &[Bond] = if !structure.bonds.is_empty() {
        &structure.bonds
    } else {
        detected_bonds = BondDetector::detect_bonds(&structure.atoms);
        &detected_bonds
    };

    for bond in bonds {
        if bond.atom1_idx < structure.atoms.len() && bond.atom2_idx < structure.atoms.len() {
            let atom1 = &structure.atoms[bond.atom1_idx];
            let atom2 = &structure.atoms[bond.atom2_idx];

            let p1_opt = camera.world_to_screen(atom1.pos, buffer.width, buffer.height);
            let p2_opt = camera.world_to_screen(atom2.pos, buffer.width, buffer.height);

            if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
                let avg_depth = (p1.2 + p2.2) * 0.5;
                let bond_r = project_radius(0.18, avg_depth, camera.fov, buffer.height).max(0.5);

                let res1 = residue_map.get(&atom1.index).copied();
                let res2 = residue_map.get(&atom2.index).copied();
                let c1 = color_for_atom(atom1, res1, structure, color_scheme);
                let c2 = color_for_atom(atom2, res2, structure, color_scheme);

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
