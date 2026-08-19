//! Secondary structure cartoon ribbon representation.
//!
//! Generates smooth Catmull-Rom spline curves along polymer chains with secondary structure styling:
//! - **Alpha Helix**: Thick shaded cylinder / spiral tube ($R \approx 1.2\text{ \AA}$).
//! - **Beta Sheet**: Flat ribbon quad-strip with arrowhead pointing towards C-terminal end.
//! - **Random Coil / Loop**: Smooth thin tube ($R \approx 0.3\text{ \AA}$).
//! - **Heteroatoms / Ligands**: Rendered in Ball & Stick style.

use crate::math::Vec3;
use crate::math::spline::CatmullRomSpline;
use crate::model::Structure;
use crate::model::atom::Atom;
use crate::model::residue::{Residue, SecondaryStructure};
use crate::render::buffer::Framebuffer;
use crate::render::camera::Camera;
use crate::render::color::{ColorScheme, color_for_atom};
use crate::render::lighting::Lighting;
use crate::render::rasterizer::{draw_cylinder, draw_line_3d, draw_sphere, draw_triangle_3d};
use crate::render::representations::{build_atom_residue_map, lerp_rgb, project_radius};

/// Maximum distance in Å between consecutive guide atoms before splitting chain segments.
const MAX_RIBBON_SEGMENT_DISTANCE: f32 = 8.0;

/// Number of spline interpolation samples per residue step.
const SAMPLES_PER_RESIDUE: usize = 6;

/// Information for a single residue guide point along a polymer chain.
struct ResidueGuide<'a> {
    atom: &'a Atom,
    residue: &'a Residue,
    pos: Vec3,
    normal: Vec3,
    ss: SecondaryStructure,
}

/// Finds the guide atom for ribbon spline extraction in a given residue.
fn find_ribbon_guide_atom<'a>(residue: &'a Residue, atoms: &'a [Atom]) -> Option<&'a Atom> {
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

/// Computes an orientation normal vector for a residue (e.g. from CA to carbonyl O or N).
fn compute_residue_normal(residue: &Residue, ca_pos: Vec3, atoms: &[Atom]) -> Vec3 {
    let mut o_pos = None;
    let mut n_pos = None;

    for &idx in &residue.atom_indices {
        if let Some(atom) = atoms.get(idx) {
            let name = atom.name.trim();
            if name.eq_ignore_ascii_case("O") || name.eq_ignore_ascii_case("OXT") {
                o_pos = Some(atom.pos);
            } else if name.eq_ignore_ascii_case("N") {
                n_pos = Some(atom.pos);
            }
        }
    }

    if let Some(o) = o_pos {
        let v = o - ca_pos;
        if v.norm() > 1e-4 {
            return v.normalize();
        }
    }

    if let Some(n) = n_pos {
        let v = ca_pos - n;
        if v.norm() > 1e-4 {
            return v.normalize();
        }
    }

    Vec3::new(0.0, 1.0, 0.0)
}

/// Returns a stable perpendicular vector for a normalized tangent.
fn fallback_perp(t: Vec3) -> Vec3 {
    let up = if t.y.abs() < 0.9 { Vec3::Y } else { Vec3::X };
    let perp = t.cross(up);
    if perp.norm() > 1e-4 {
        perp.normalize()
    } else {
        Vec3::Z
    }
}

/// Renders the molecular structure in Ribbon (secondary structure cartoon) mode.
pub fn render_ribbon(
    structure: &Structure,
    color_scheme: ColorScheme,
    camera: &Camera,
    buffer: &mut Framebuffer,
    lighting: &Lighting,
) {
    let residue_map = build_atom_residue_map(structure);

    // 1. Render Cartoon Ribbon for each polymer chain
    for chain in &structure.chains {
        let mut all_guides: Vec<ResidueGuide> = Vec::with_capacity(chain.residues.len());

        for res in &chain.residues {
            if let Some(guide_atom) = find_ribbon_guide_atom(res, &structure.atoms) {
                let norm = compute_residue_normal(res, guide_atom.pos, &structure.atoms);
                all_guides.push(ResidueGuide {
                    atom: guide_atom,
                    residue: res,
                    pos: guide_atom.pos,
                    normal: norm,
                    ss: res.secondary_structure,
                });
            }
        }

        if all_guides.is_empty() {
            continue;
        }

        // Split into contiguous segments based on distance gaps
        let mut segments: Vec<Vec<ResidueGuide>> = Vec::new();
        let mut current_segment = Vec::new();

        for guide in all_guides {
            if let Some(last) = current_segment.last() {
                let last_guide: &ResidueGuide = last;
                if last_guide.pos.distance(&guide.pos) > MAX_RIBBON_SEGMENT_DISTANCE {
                    let prev = std::mem::take(&mut current_segment);
                    segments.push(prev);
                }
            }
            current_segment.push(guide);
        }
        if !current_segment.is_empty() {
            segments.push(current_segment);
        }

        // Render each contiguous segment
        for segment in &segments {
            render_ribbon_segment(segment, structure, color_scheme, camera, buffer, lighting);
        }
    }

    // 2. Render HETATM / ligands in Ball & Stick style
    for atom in &structure.atoms {
        if atom.is_hetatm {
            let Some(pt) = camera.world_to_screen(atom.pos, buffer.width, buffer.height) else {
                continue;
            };
            let r_world = (atom.vdw_radius() * 0.28).clamp(0.35, 0.65);
            let sphere_r = project_radius(r_world, pt.2, camera.fov, buffer.height).max(0.8);
            let res = residue_map.get(&atom.index).copied();
            let color = color_for_atom(atom, res, structure, color_scheme);
            draw_sphere(buffer, pt, sphere_r, color, lighting);
        }
    }

    // 3. Render bonds involving HETATMs
    for bond in &structure.bonds {
        if bond.atom1_idx < structure.atoms.len() && bond.atom2_idx < structure.atoms.len() {
            let atom1 = &structure.atoms[bond.atom1_idx];
            let atom2 = &structure.atoms[bond.atom2_idx];

            if atom1.is_hetatm || atom2.is_hetatm {
                let p1_opt = camera.world_to_screen(atom1.pos, buffer.width, buffer.height);
                let p2_opt = camera.world_to_screen(atom2.pos, buffer.width, buffer.height);

                if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
                    let avg_depth = (p1.2 + p2.2) * 0.5;
                    let cyl_r = project_radius(0.18, avg_depth, camera.fov, buffer.height).max(0.5);

                    let res1 = residue_map.get(&atom1.index).copied();
                    let res2 = residue_map.get(&atom2.index).copied();
                    let c1 = color_for_atom(atom1, res1, structure, color_scheme);
                    let c2 = color_for_atom(atom2, res2, structure, color_scheme);

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
}

/// Renders a single contiguous chain segment using Catmull-Rom spline interpolation.
fn render_ribbon_segment(
    guides: &[ResidueGuide],
    structure: &Structure,
    color_scheme: ColorScheme,
    camera: &Camera,
    buffer: &mut Framebuffer,
    lighting: &Lighting,
) {
    let n = guides.len();
    if n == 0 {
        return;
    }

    if n == 1 {
        let g = &guides[0];
        if let Some(pt) = camera.world_to_screen(g.pos, buffer.width, buffer.height) {
            let sphere_r = project_radius(0.5, pt.2, camera.fov, buffer.height).max(0.8);
            let color = color_for_atom(g.atom, Some(g.residue), structure, color_scheme);
            draw_sphere(buffer, pt, sphere_r, color, lighting);
        }
        return;
    }

    let points: Vec<Vec3> = guides.iter().map(|g| g.pos).collect();
    let colors: Vec<_> = guides
        .iter()
        .map(|g| color_for_atom(g.atom, Some(g.residue), structure, color_scheme))
        .collect();

    let num_segments = n - 1;

    for seg in 0..num_segments {
        let g1 = &guides[seg];
        let g2 = &guides[seg + 1];
        let c1 = colors[seg];
        let c2 = colors[seg + 1];

        let p1 = points[seg];
        let p2 = points[seg + 1];
        let p0 = if seg > 0 {
            points[seg - 1]
        } else {
            p1 * 2.0 - p2
        };
        let p3 = if seg + 2 < n {
            points[seg + 2]
        } else {
            p2 * 2.0 - p1
        };

        let is_strand_c_term = g1.ss == SecondaryStructure::Sheet
            && (seg + 2 >= n || guides[seg + 2].ss != SecondaryStructure::Sheet);

        let k_samples = SAMPLES_PER_RESIDUE;

        for k in 0..k_samples {
            let u0 = (k as f32) / (k_samples as f32);
            let u1 = ((k + 1) as f32) / (k_samples as f32);
            let u_mid = (u0 + u1) * 0.5;

            let pt0 = CatmullRomSpline::interpolate_segment(p0, p1, p2, p3, u0);
            let pt1 = CatmullRomSpline::interpolate_segment(p0, p1, p2, p3, u1);

            let tan0 = CatmullRomSpline::tangent_segment(p0, p1, p2, p3, u0);
            let tan1 = CatmullRomSpline::tangent_segment(p0, p1, p2, p3, u1);

            let ss = if u_mid < 0.5 { g1.ss } else { g2.ss };
            let color = lerp_rgb(c1, c2, u_mid);

            match ss {
                SecondaryStructure::Helix => {
                    // Alpha-Helix: thick shaded cylinder
                    let p1_opt = camera.world_to_screen(pt0, buffer.width, buffer.height);
                    let p2_opt = camera.world_to_screen(pt1, buffer.width, buffer.height);

                    if let (Some(s0), Some(s1)) = (p1_opt, p2_opt) {
                        let avg_depth = (s0.2 + s1.2) * 0.5;
                        let cyl_r =
                            project_radius(1.2, avg_depth, camera.fov, buffer.height).max(1.2);
                        draw_cylinder(buffer, s0, s1, cyl_r, color, lighting);
                    }
                }

                SecondaryStructure::Sheet => {
                    // Beta-Sheet: flat ribbon quad-strip with arrowhead
                    let n0 = g1.normal.lerp(g2.normal, u0);
                    let mut w0_vec = tan0.cross(n0);
                    if w0_vec.norm() > 1e-4 {
                        w0_vec = w0_vec.normalize();
                    } else {
                        w0_vec = fallback_perp(tan0);
                    }
                    let ribbon_norm0 = w0_vec.cross(tan0).normalize();

                    let n1 = g1.normal.lerp(g2.normal, u1);
                    let mut w1_vec = tan1.cross(n1);
                    if w1_vec.norm() > 1e-4 {
                        w1_vec = w1_vec.normalize();
                    } else {
                        w1_vec = fallback_perp(tan1);
                    }
                    let ribbon_norm1 = w1_vec.cross(tan1).normalize();

                    // Calculate width factoring in arrowhead
                    let w0 = if is_strand_c_term {
                        if u0 < 0.35 {
                            2.2 + (u0 / 0.35) * 1.0
                        } else {
                            (3.2 * (1.0 - (u0 - 0.35) / 0.65)).max(0.1)
                        }
                    } else {
                        2.2
                    };

                    let w1 = if is_strand_c_term {
                        if u1 < 0.35 {
                            2.2 + (u1 / 0.35) * 1.0
                        } else {
                            (3.2 * (1.0 - (u1 - 0.35) / 0.65)).max(0.1)
                        }
                    } else {
                        2.2
                    };

                    let v0_l = pt0 - w0_vec * (w0 * 0.5);
                    let v0_r = pt0 + w0_vec * (w0 * 0.5);
                    let v1_l = pt1 - w1_vec * (w1 * 0.5);
                    let v1_r = pt1 + w1_vec * (w1 * 0.5);

                    let s0_l = camera.world_to_screen(v0_l, buffer.width, buffer.height);
                    let s0_r = camera.world_to_screen(v0_r, buffer.width, buffer.height);
                    let s1_l = camera.world_to_screen(v1_l, buffer.width, buffer.height);
                    let s1_r = camera.world_to_screen(v1_r, buffer.width, buffer.height);

                    if let (Some(sl0), Some(sr0), Some(sl1), Some(sr1)) = (s0_l, s0_r, s1_l, s1_r) {
                        let n_avg = (ribbon_norm0 + ribbon_norm1).normalize();

                        // Top face
                        draw_triangle_3d(buffer, sl0, sr0, sr1, n_avg, color, lighting);
                        draw_triangle_3d(buffer, sl0, sr1, sl1, n_avg, color, lighting);

                        // Bottom face (two-sided lighting)
                        draw_triangle_3d(buffer, sr0, sl0, sr1, -n_avg, color, lighting);
                        draw_triangle_3d(buffer, sr1, sl0, sl1, -n_avg, color, lighting);

                        // Edge lines for crisp visual definition
                        draw_line_3d(buffer, sl0, sl1, color);
                        draw_line_3d(buffer, sr0, sr1, color);
                    }
                }

                SecondaryStructure::Coil => {
                    // Random coil / loop: thin smooth tube
                    let p1_opt = camera.world_to_screen(pt0, buffer.width, buffer.height);
                    let p2_opt = camera.world_to_screen(pt1, buffer.width, buffer.height);

                    if let (Some(s0), Some(s1)) = (p1_opt, p2_opt) {
                        let avg_depth = (s0.2 + s1.2) * 0.5;
                        let cyl_r =
                            project_radius(0.3, avg_depth, camera.fov, buffer.height).max(0.6);
                        draw_cylinder(buffer, s0, s1, cyl_r, color, lighting);
                    }
                }
            }
        }
    }
}
