//! Backbone / C-alpha trace representation.
//!
//! Connects sequential C-alpha (CA) and nucleic acid phosphorus (P) backbone guide atoms
//! with smooth 3D cylinders, and renders heteroatom ligands/cofactors as analytical spheres.

use crate::model::atom::Atom;
use crate::model::residue::Residue;
use crate::render::buffer::Framebuffer;
use crate::render::rasterizer::{draw_cylinder, draw_sphere};
use crate::render::representations::{RenderContext, project_radius};

/// Maximum distance in A between consecutive guide atoms before assuming a chain break.
pub const MAX_TRACE_BOND_DISTANCE: f32 = 8.0;

/// Finds the guide atom for trace rendering in a given residue.
/// Returns C-alpha (CA) for amino acids, Phosphorus (P) for nucleic acids, or first atom.
pub fn find_trace_guide_atom<'a>(residue: &'a Residue, atoms: &'a [Atom]) -> Option<&'a Atom> {
    if let Some(ca) = residue.ca_atom(atoms) {
        return Some(ca);
    }

    if residue.is_nucleic() {
        for &idx in &residue.atom_indices {
            if let Some(atom) = atoms.get(idx) {
                let name = atom.name.trim();
                if name.eq_ignore_ascii_case("P") {
                    return Some(atom);
                }
            }
        }
    }

    for &idx in &residue.atom_indices {
        if let Some(atom) = atoms.get(idx).filter(|a| a.is_backbone()) {
            return Some(atom);
        }
    }

    residue.atom_indices.first().and_then(|&idx| atoms.get(idx))
}

/// Renders the molecular structure in Trace mode.
pub fn render_trace(ctx: &RenderContext, buffer: &mut Framebuffer) {
    let structure = ctx.structure;
    let camera = ctx.camera;
    let mats = ctx.mats;
    let lighting = ctx.lighting;
    let visibility = ctx.visibility;
    let colors = ctx.colors;
    let visible = ctx.visible;

    // 1. Render backbone trace for each chain
    for chain in structure.chains() {
        let mut guide_list: Vec<(&Atom, &Residue)> = Vec::with_capacity(chain.residues.len());

        for res in &chain.residues {
            if !visibility.residue_visible(res) {
                continue;
            }
            if let Some(guide_atom) = find_trace_guide_atom(res, structure.atoms()) {
                guide_list.push((guide_atom, res));
            }
        }

        // Draw connecting cylinders between consecutive guide atoms
        for window in guide_list.windows(2) {
            let (atom1, _) = window[0];
            let (atom2, _) = window[1];

            if atom1.pos.distance(&atom2.pos) > MAX_TRACE_BOND_DISTANCE {
                continue;
            }

            let p1_opt = camera.project(&mats, atom1.pos, buffer.width, buffer.height);
            let p2_opt = camera.project(&mats, atom2.pos, buffer.width, buffer.height);

            if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
                let avg_depth = (p1.2 + p2.2) * 0.5;
                let cyl_r = project_radius(0.35, avg_depth, camera.fov, buffer.height).max(0.6);

                let c1 = colors[atom1.index];
                let c2 = colors[atom2.index];

                if c1 == c2 {
                    draw_cylinder(buffer, p1, p2, cyl_r, c1, lighting);
                } else {
                    let pmid = ((p1.0 + p2.0) * 0.5, (p1.1 + p2.1) * 0.5, avg_depth);
                    draw_cylinder(buffer, p1, pmid, cyl_r, c1, lighting);
                    draw_cylinder(buffer, pmid, p2, cyl_r, c2, lighting);
                }
            }
        }

        // Draw joint spheres at guide atom positions
        for (atom, _) in &guide_list {
            if let Some(pt) = camera.project(&mats, atom.pos, buffer.width, buffer.height) {
                let sphere_r = project_radius(0.40, pt.2, camera.fov, buffer.height).max(0.7);
                let color = colors[atom.index];
                draw_sphere(buffer, pt, sphere_r, color, lighting);
            }
        }
    }

    // 2. Render HETATM / ligand atoms as analytical spheres
    for atom in structure.atoms() {
        if atom.is_hetatm {
            if !visible[atom.index] {
                continue;
            }
            let Some(pt) = camera.project(&mats, atom.pos, buffer.width, buffer.height) else {
                continue;
            };
            let r_world = (atom.vdw_radius() * 0.28).clamp(0.35, 0.65);
            let sphere_r = project_radius(r_world, pt.2, camera.fov, buffer.height).max(0.7);
            let color = colors[atom.index];
            draw_sphere(buffer, pt, sphere_r, color, lighting);
        }
    }

    // 3. Render bonds between HETATMs if present
    for bond in structure.bonds() {
        if bond.atom1_idx < structure.atoms().len() && bond.atom2_idx < structure.atoms().len() {
            let atom1 = &structure.atoms()[bond.atom1_idx];
            let atom2 = &structure.atoms()[bond.atom2_idx];
            // Skip non-ligand bonds early (cheap field check) so large structures
            // don't pay O(bonds) in residue lookup / visibility checks here.
            if !atom1.is_hetatm && !atom2.is_hetatm {
                continue;
            }
            if !visible[atom1.index] || !visible[atom2.index] {
                continue;
            }

            let p1_opt = camera.project(&mats, atom1.pos, buffer.width, buffer.height);
            let p2_opt = camera.project(&mats, atom2.pos, buffer.width, buffer.height);

            if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
                let avg_depth = (p1.2 + p2.2) * 0.5;
                let cyl_r = project_radius(0.18, avg_depth, camera.fov, buffer.height).max(0.5);

                let c1 = colors[atom1.index];
                let c2 = colors[atom2.index];

                if c1 == c2 {
                    draw_cylinder(buffer, p1, p2, cyl_r, c1, lighting);
                } else {
                    let pmid = ((p1.0 + p2.0) * 0.5, (p1.1 + p2.1) * 0.5, avg_depth);
                    draw_cylinder(buffer, p1, pmid, cyl_r, c1, lighting);
                    draw_cylinder(buffer, pmid, p2, cyl_r, c2, lighting);
                }
            }
        }
    }
}
