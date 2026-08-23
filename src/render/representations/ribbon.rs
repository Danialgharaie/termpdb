//! Secondary structure cartoon ribbon representation.
//!
//! Generates smooth Catmull-Rom spline curves along polymer chains with secondary structure styling:
//! - Alpha Helix: thick shaded cylinder / spiral tube (R ~ 1.2 A).
//! - Beta Sheet: flat ribbon quad-strip with arrowhead pointing towards C-terminal end.
//! - Random Coil / Loop: thin smooth tube (R ~ 0.3 A).
//! - Heteroatoms / Ligands: rendered in Ball & Stick style.
//!
//! The camera-independent tessellation (guide atoms, spline samples, sheet
//! geometry, ligand primitives) is built once into a Vec<RibbonPrimitive> and
//! re-projected each frame, so orbit/spin skip the O(residues) guide building,
//! the Catmull-Rom evaluation, and the O(atoms)/O(bonds) ligand scans.

use crate::math::Vec3;
use crate::math::spline::CatmullRomSpline;
use crate::model::atom::Atom;
use crate::model::residue::{Residue, SecondaryStructure};
use crate::render::buffer::{Framebuffer, PixelColor};
use crate::render::rasterizer::{draw_cylinder, draw_line_3d, draw_sphere, draw_triangle_3d};
use crate::render::representations::{LodLevel, RenderContext, project_radius};

/// Maximum distance in A between consecutive guide atoms before splitting chain segments.
const MAX_RIBBON_SEGMENT_DISTANCE: f32 = 8.0;

/// A camera-independent ribbon primitive in world space. Built once per
/// structure/LOD/color/visibility change; render_ribbon just projects + draws it.
pub enum RibbonPrimitive {
    /// Shaded tube segment (alpha-helix or coil).
    Cylinder {
        a: Vec3,
        b: Vec3,
        r_world: f32,
        min_r: f32,
        color: PixelColor,
    },
    /// Beta-sheet quad strip segment: four world corners + averaged face normal.
    /// Drawn as two top triangles (+normal), two bottom triangles (-normal), and
    /// two edge lines for crisp definition.
    SheetQuad {
        v0l: Vec3,
        v0r: Vec3,
        v1l: Vec3,
        v1r: Vec3,
        normal: Vec3,
        color: PixelColor,
    },
    /// Joint / ligand sphere.
    Sphere {
        c: Vec3,
        r_world: f32,
        min_r: f32,
        color: PixelColor,
    },
}

/// Information for a single residue guide point along a polymer chain.
struct ResidueGuide<'a> {
    atom: &'a Atom,
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

/// Beta-sheet half-width (with C-terminal arrowhead flare) at parametric position u.
fn sheet_half_width(u: f32, is_strand_c_term: bool) -> f32 {
    if is_strand_c_term {
        if u < 0.35 {
            2.2 + (u / 0.35) * 1.0
        } else {
            (3.2 * (1.0 - (u - 0.35) / 0.65)).max(0.1)
        }
    } else {
        2.2
    }
}

/// Builds the camera-independent ribbon geometry (spline + ligand primitives) in
/// one pass. Callers cache this and rebuild only when the structure, color
/// scheme, visibility, or LOD changes.
pub fn build_ribbon_geometry(
    structure: &crate::model::Structure,
    colors: &[PixelColor],
    visible: &[bool],
    visibility: crate::render::Visibility,
    lod: LodLevel,
) -> Vec<RibbonPrimitive> {
    let atoms = structure.atoms();
    let samples = lod.ribbon_samples(structure.residue_count());
    let mut out: Vec<RibbonPrimitive> = Vec::new();

    // 1. Cartoon ribbon for each polymer chain
    for chain in structure.chains() {
        let mut all_guides: Vec<ResidueGuide> = Vec::with_capacity(chain.residues.len());
        for res in &chain.residues {
            if !visibility.residue_visible(res) {
                continue;
            }
            if let Some(guide_atom) = find_ribbon_guide_atom(res, atoms) {
                let norm = compute_residue_normal(res, guide_atom.pos, atoms);
                all_guides.push(ResidueGuide {
                    atom: guide_atom,
                    pos: guide_atom.pos,
                    normal: norm,
                    ss: res.secondary_structure,
                });

                if (res.is_nucleic()
                    || crate::render::representations::nucleic::is_nucleic_residue(&res.name))
                    && let Some(slab) = crate::render::representations::nucleic::build_base_slab(
                        res,
                        atoms,
                        guide_atom.pos,
                    )
                {
                    out.push(slab);
                }
            }
        }
        if all_guides.is_empty() {
            continue;
        }

        let mut segments: Vec<Vec<ResidueGuide>> = Vec::new();
        let mut current_segment = Vec::new();
        for guide in all_guides {
            if let Some(last) = current_segment.last() {
                let last_guide: &ResidueGuide = last;
                if last_guide.pos.distance(&guide.pos) > MAX_RIBBON_SEGMENT_DISTANCE {
                    segments.push(std::mem::take(&mut current_segment));
                }
            }
            current_segment.push(guide);
        }
        if !current_segment.is_empty() {
            segments.push(current_segment);
        }

        for segment in &segments {
            emit_segment_primitives(segment, colors, samples, &mut out);
        }
    }

    // 2. HETATM / ligand spheres (Ball & Stick style)
    for atom in atoms {
        if atom.is_hetatm && visible[atom.index] {
            let r_world = (atom.vdw_radius() * 0.28).clamp(0.35, 0.65);
            out.push(RibbonPrimitive::Sphere {
                c: atom.pos,
                r_world,
                min_r: 0.8,
                color: colors[atom.index],
            });
        }
    }

    // 3. HETATM bonds
    for bond in structure.bonds() {
        if bond.atom1_idx < atoms.len() && bond.atom2_idx < atoms.len() {
            let atom1 = &atoms[bond.atom1_idx];
            let atom2 = &atoms[bond.atom2_idx];
            if !atom1.is_hetatm && !atom2.is_hetatm {
                continue;
            }
            if !visible[atom1.index] || !visible[atom2.index] {
                continue;
            }
            let c1 = colors[atom1.index];
            let c2 = colors[atom2.index];
            if c1 == c2 {
                out.push(RibbonPrimitive::Cylinder {
                    a: atom1.pos,
                    b: atom2.pos,
                    r_world: 0.18,
                    min_r: 0.5,
                    color: c1,
                });
            } else {
                let pmid = (atom1.pos + atom2.pos) * 0.5;
                out.push(RibbonPrimitive::Cylinder {
                    a: atom1.pos,
                    b: pmid,
                    r_world: 0.18,
                    min_r: 0.5,
                    color: c1,
                });
                out.push(RibbonPrimitive::Cylinder {
                    a: pmid,
                    b: atom2.pos,
                    r_world: 0.18,
                    min_r: 0.5,
                    color: c2,
                });
            }
        }
    }

    out
}

/// Tessellates one contiguous chain segment into world-space primitives.
fn emit_segment_primitives(
    guides: &[ResidueGuide],
    colors: &[PixelColor],
    samples: usize,
    out: &mut Vec<RibbonPrimitive>,
) {
    let n = guides.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        let g = &guides[0];
        out.push(RibbonPrimitive::Sphere {
            c: g.pos,
            r_world: 0.5,
            min_r: 0.8,
            color: colors[g.atom.index],
        });
        return;
    }

    let points: Vec<Vec3> = guides.iter().map(|g| g.pos).collect();
    let seg_colors: Vec<PixelColor> = guides.iter().map(|g| colors[g.atom.index]).collect();
    let k_samples = samples.max(1);

    for seg in 0..(n - 1) {
        let g1 = &guides[seg];
        let g2 = &guides[seg + 1];
        let c1 = seg_colors[seg];
        let c2 = seg_colors[seg + 1];

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

        for k in 0..k_samples {
            let u0 = (k as f32) / (k_samples as f32);
            let u1 = ((k + 1) as f32) / (k_samples as f32);
            let u_mid = (u0 + u1) * 0.5;

            let pt0 = CatmullRomSpline::interpolate_segment(p0, p1, p2, p3, u0);
            let pt1 = CatmullRomSpline::interpolate_segment(p0, p1, p2, p3, u1);
            let tan0 = CatmullRomSpline::tangent_segment(p0, p1, p2, p3, u0);
            let tan1 = CatmullRomSpline::tangent_segment(p0, p1, p2, p3, u1);

            let ss = if u_mid < 0.5 { g1.ss } else { g2.ss };
            let color = crate::render::representations::lerp_rgb(c1, c2, u_mid);

            match ss {
                SecondaryStructure::Helix => out.push(RibbonPrimitive::Cylinder {
                    a: pt0,
                    b: pt1,
                    r_world: 1.2,
                    min_r: 1.2,
                    color,
                }),
                SecondaryStructure::Sheet => {
                    let n0 = g1.normal.lerp(g2.normal, u0);
                    let mut w0_vec = tan0.cross(n0);
                    w0_vec = if w0_vec.norm() > 1e-4 {
                        w0_vec.normalize()
                    } else {
                        fallback_perp(tan0)
                    };
                    let ribbon_norm0 = w0_vec.cross(tan0).normalize();

                    let n1 = g1.normal.lerp(g2.normal, u1);
                    let mut w1_vec = tan1.cross(n1);
                    w1_vec = if w1_vec.norm() > 1e-4 {
                        w1_vec.normalize()
                    } else {
                        fallback_perp(tan1)
                    };
                    let ribbon_norm1 = w1_vec.cross(tan1).normalize();

                    let hw0 = sheet_half_width(u0, is_strand_c_term);
                    let hw1 = sheet_half_width(u1, is_strand_c_term);

                    out.push(RibbonPrimitive::SheetQuad {
                        v0l: pt0 - w0_vec * (hw0 * 0.5),
                        v0r: pt0 + w0_vec * (hw0 * 0.5),
                        v1l: pt1 - w1_vec * (hw1 * 0.5),
                        v1r: pt1 + w1_vec * (hw1 * 0.5),
                        normal: (ribbon_norm0 + ribbon_norm1).normalize(),
                        color,
                    });
                }
                SecondaryStructure::Coil => out.push(RibbonPrimitive::Cylinder {
                    a: pt0,
                    b: pt1,
                    r_world: 0.3,
                    min_r: 0.6,
                    color,
                }),
            }
        }
    }
}

/// Renders the molecular structure in Ribbon mode.
///
/// Uses the cached geometry from ctx.ribbon_geometry when available (interactive
/// path); otherwise builds it fresh (one-shot / export path).
pub fn render_ribbon(ctx: &RenderContext, buffer: &mut Framebuffer) {
    let camera = ctx.camera;
    let mats = ctx.mats;
    let lighting = ctx.lighting;

    let fresh;
    let geometry: &[RibbonPrimitive] = match ctx.ribbon_geometry {
        Some(g) => g,
        None => {
            fresh = build_ribbon_geometry(
                ctx.structure,
                ctx.colors,
                ctx.visible,
                ctx.visibility,
                ctx.lod,
            );
            &fresh
        }
    };

    for prim in geometry {
        match prim {
            RibbonPrimitive::Cylinder {
                a,
                b,
                r_world,
                min_r,
                color,
            } => {
                let pa = camera.project(&mats, *a, buffer.width, buffer.height);
                let pb = camera.project(&mats, *b, buffer.width, buffer.height);
                if let (Some(sa), Some(sb)) = (pa, pb) {
                    let avg_depth = (sa.2 + sb.2) * 0.5;
                    let cyl_r =
                        project_radius(*r_world, avg_depth, camera.fov, buffer.height).max(*min_r);
                    draw_cylinder(buffer, sa, sb, cyl_r, *color, lighting);
                }
            }
            RibbonPrimitive::SheetQuad {
                v0l,
                v0r,
                v1l,
                v1r,
                normal,
                color,
            } => {
                let sl0 = camera.project(&mats, *v0l, buffer.width, buffer.height);
                let sr0 = camera.project(&mats, *v0r, buffer.width, buffer.height);
                let sl1 = camera.project(&mats, *v1l, buffer.width, buffer.height);
                let sr1 = camera.project(&mats, *v1r, buffer.width, buffer.height);
                if let (Some(sl0), Some(sr0), Some(sl1), Some(sr1)) = (sl0, sr0, sl1, sr1) {
                    draw_triangle_3d(buffer, sl0, sr0, sr1, *normal, *color, lighting);
                    draw_triangle_3d(buffer, sl0, sr1, sl1, *normal, *color, lighting);
                    draw_triangle_3d(buffer, sr0, sl0, sr1, -*normal, *color, lighting);
                    draw_triangle_3d(buffer, sr1, sl0, sl1, -*normal, *color, lighting);
                    draw_line_3d(buffer, sl0, sl1, *color);
                    draw_line_3d(buffer, sr0, sr1, *color);
                }
            }
            RibbonPrimitive::Sphere {
                c,
                r_world,
                min_r,
                color,
            } => {
                if let Some(pt) = camera.project(&mats, *c, buffer.width, buffer.height) {
                    let sphere_r =
                        project_radius(*r_world, pt.2, camera.fov, buffer.height).max(*min_r);
                    draw_sphere(buffer, pt, sphere_r, *color, lighting);
                }
            }
        }
    }
}
