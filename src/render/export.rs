//! Offline high-quality rendering: supersampled PNG, vector SVG, and MP4 (via ffmpeg).

use std::fs::File;
use std::io::BufWriter;
use std::process::Command;

use crate::error::{Result, TermPdbError};
use crate::math::Vec3;
use crate::model::Structure;
use crate::render::PixelColor;
use crate::render::buffer::Framebuffer;
use crate::render::camera::{Camera, CameraMatrices};
use crate::render::color::ColorScheme;
use crate::render::fit_render_size;
use crate::render::lighting::Lighting;
use crate::render::representations::trace::{MAX_TRACE_BOND_DISTANCE, find_trace_guide_atom};
use crate::render::representations::{
    LodMode, RenderMode, RibbonPrimitive, Visibility, build_render_cache, build_ribbon_geometry,
    project_radius, render_structure,
};

fn fit_camera(structure: &Structure, sw: usize, sh: usize) -> Camera {
    let mut camera = Camera::new();
    camera.aspect = (sw as f32) / (sh as f32);
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);
    camera
}

#[allow(clippy::too_many_arguments)]
pub fn render_supersampled(
    structure: &Structure,
    mode: RenderMode,
    color: ColorScheme,
    width: usize,
    height: usize,
    ssaa: usize,
    visibility: Visibility,
    lod: LodMode,
) -> Vec<u8> {
    let (width, height, ssaa) = fit_render_size(width, height, ssaa);
    let ssaa = ssaa.max(1);
    let sw = width * ssaa;
    let sh = height * ssaa;
    let mut fb = Framebuffer::new(sw, sh);
    let camera = fit_camera(structure, sw, sh);
    let lighting = Lighting::default();
    fb.clear((0, 0, 0));
    render_structure(
        structure, mode, color, &camera, &mut fb, &lighting, visibility, lod,
    );
    downsample_rgba(&fb, width, height, ssaa)
}

#[allow(clippy::manual_checked_ops)]
pub fn downsample_rgba(fb: &Framebuffer, width: usize, height: usize, ssaa: usize) -> Vec<u8> {
    let sw = fb.width;
    let sh = fb.height;
    let total = (ssaa * ssaa) as u32;
    let mut out = vec![0u8; width * height * 4];
    for y in 0..height {
        for x in 0..width {
            let mut r = 0u32;
            let mut g = 0u32;
            let mut b = 0u32;
            let mut drawn = 0u32;
            for dy in 0..ssaa {
                for dx in 0..ssaa {
                    let sx = x * ssaa + dx;
                    let sy = y * ssaa + dy;
                    if sx < sw && sy < sh {
                        let idx = sy * sw + sx;
                        if fb.depth[idx].is_finite() {
                            let p = fb.pixels[idx];
                            r += p.0 as u32;
                            g += p.1 as u32;
                            b += p.2 as u32;
                            drawn += 1;
                        }
                    }
                }
            }
            let o = (y * width + x) * 4;
            if drawn > 0 {
                out[o] = (r / drawn) as u8;
                out[o + 1] = (g / drawn) as u8;
                out[o + 2] = (b / drawn) as u8;
                out[o + 3] = ((drawn * 255) / total).min(255) as u8;
            }
        }
    }
    out
}

pub fn write_png(path: &str, rgba: &[u8], width: u32, height: u32) -> Result<()> {
    let file = File::create(path)?;
    let w = BufWriter::new(file);
    let mut enc = png::Encoder::new(w, width, height);
    enc.set_color(png::ColorType::Rgba);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc
        .write_header()
        .map_err(|e| TermPdbError::RenderError(format!("png encode: {e}")))?;
    writer
        .write_image_data(rgba)
        .map_err(|e| TermPdbError::RenderError(format!("png encode: {e}")))?;
    Ok(())
}

enum SvgPrim {
    Circle {
        cx: f32,
        cy: f32,
        r: f32,
        depth: f32,
        color: PixelColor,
    },
    Line {
        x1: f32,
        y1: f32,
        x2: f32,
        y2: f32,
        w: f32,
        depth: f32,
        color: PixelColor,
    },
    Polygon {
        pts: Vec<(f32, f32)>,
        depth: f32,
        color: PixelColor,
    },
}

impl SvgPrim {
    fn depth(&self) -> f32 {
        match self {
            SvgPrim::Circle { depth, .. }
            | SvgPrim::Line { depth, .. }
            | SvgPrim::Polygon { depth, .. } => *depth,
        }
    }
    fn color(&self) -> PixelColor {
        match self {
            SvgPrim::Circle { color, .. }
            | SvgPrim::Line { color, .. }
            | SvgPrim::Polygon { color, .. } => *color,
        }
    }
    fn emit(&self) -> String {
        let (cr, cg, cb) = self.color();
        let fill = format!("rgb({},{},{})", cr, cg, cb);
        match self {
            SvgPrim::Circle { cx, cy, r, .. } => {
                format!(
                    r#"<circle cx="{:.2}" cy="{:.2}" r="{:.2}" fill="{}"/>"#,
                    cx, cy, r, fill
                )
            }
            SvgPrim::Line {
                x1, y1, x2, y2, w, ..
            } => format!(
                r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}" stroke-width="{:.2}" stroke-linecap="round"/>"#,
                x1, y1, x2, y2, fill, w
            ),
            SvgPrim::Polygon { pts, .. } => {
                let pts: Vec<String> = pts
                    .iter()
                    .map(|(x, y)| format!("{:.2},{:.2}", x, y))
                    .collect();
                format!(r#"<polygon points="{}" fill="{}"/>"#, pts.join(" "), fill)
            }
        }
    }
}

pub fn render_svg(
    structure: &Structure,
    mode: RenderMode,
    color: ColorScheme,
    width: usize,
    height: usize,
    visibility: Visibility,
    lod: LodMode,
) -> String {
    let camera = fit_camera(structure, width, height);
    let lighting = Lighting::default();
    let mats = camera.matrices();
    let (colors, visible, _, _, _) = build_render_cache(structure, color, visibility, lod);
    let atoms = structure.atoms();
    let mut prims: Vec<SvgPrim> = Vec::new();

    let proj = |pos: Vec3| camera.project(&mats, pos, width, height);
    let cam_z = Vec3::new(0.0, 0.0, 1.0);

    match mode {
        RenderMode::Vdw => {
            for atom in atoms {
                if !visible[atom.index] {
                    continue;
                }
                if let Some((sx, sy, depth)) = proj(atom.pos) {
                    let rw = atom.vdw_radius();
                    let rw = if rw > 0.1 { rw } else { 1.5 };
                    let r = project_radius(rw, depth, camera.fov, height).max(1.0);
                    let lit =
                        lighting.shade(cam_z, depth, colors[atom.index], depth - rw, depth + rw);
                    prims.push(SvgPrim::Circle {
                        cx: sx,
                        cy: sy,
                        r,
                        depth,
                        color: lit,
                    });
                }
            }
        }
        RenderMode::BallAndStick => {
            for atom in atoms {
                if !visible[atom.index] {
                    continue;
                }
                if let Some((sx, sy, depth)) = proj(atom.pos) {
                    let r = project_radius(0.38, depth, camera.fov, height).max(0.8);
                    let lit = lighting.shade(
                        cam_z,
                        depth,
                        colors[atom.index],
                        depth - 0.38,
                        depth + 0.38,
                    );
                    prims.push(SvgPrim::Circle {
                        cx: sx,
                        cy: sy,
                        r,
                        depth,
                        color: lit,
                    });
                }
            }
            push_bond_lines(
                structure,
                &colors,
                &visible,
                &mats,
                &camera,
                &lighting,
                width,
                height,
                &mut prims,
                |_| true,
            );
        }
        RenderMode::Trace => {
            for chain in structure.chains() {
                let mut guides: Vec<&crate::model::atom::Atom> =
                    Vec::with_capacity(chain.residues.len());
                for res in &chain.residues {
                    if !visibility.residue_visible(res) {
                        continue;
                    }
                    if let Some(g) = find_trace_guide_atom(res, atoms) {
                        guides.push(g);
                    }
                }
                for win in guides.windows(2) {
                    let a1 = win[0];
                    let a2 = win[1];
                    if a1.pos.distance(&a2.pos) > MAX_TRACE_BOND_DISTANCE {
                        continue;
                    }
                    if let (Some(p1), Some(p2)) = (proj(a1.pos), proj(a2.pos)) {
                        let avg = (p1.2 + p2.2) * 0.5;
                        let w = project_radius(0.35, avg, camera.fov, height).max(0.6) * 2.0;
                        let c =
                            lighting.shade(cam_z, avg, colors[a1.index], avg - 0.35, avg + 0.35);
                        prims.push(SvgPrim::Line {
                            x1: p1.0,
                            y1: p1.1,
                            x2: p2.0,
                            y2: p2.1,
                            w,
                            depth: avg,
                            color: c,
                        });
                    }
                }
                for g in &guides {
                    if let Some((sx, sy, depth)) = proj(g.pos) {
                        let r = project_radius(0.40, depth, camera.fov, height).max(0.7);
                        let lit = lighting.shade(
                            cam_z,
                            depth,
                            colors[g.index],
                            depth - 0.40,
                            depth + 0.40,
                        );
                        prims.push(SvgPrim::Circle {
                            cx: sx,
                            cy: sy,
                            r,
                            depth,
                            color: lit,
                        });
                    }
                }
            }
            for atom in atoms {
                if atom.is_hetatm
                    && visible[atom.index]
                    && let Some((sx, sy, depth)) = proj(atom.pos)
                {
                    let rw = (atom.vdw_radius() * 0.28).clamp(0.35, 0.65);
                    let r = project_radius(rw, depth, camera.fov, height).max(0.7);
                    let lit =
                        lighting.shade(cam_z, depth, colors[atom.index], depth - rw, depth + rw);
                    prims.push(SvgPrim::Circle {
                        cx: sx,
                        cy: sy,
                        r,
                        depth,
                        color: lit,
                    });
                }
            }
            push_bond_lines(
                structure,
                &colors,
                &visible,
                &mats,
                &camera,
                &lighting,
                width,
                height,
                &mut prims,
                |a| a.is_hetatm,
            );
        }
        RenderMode::Ribbon => {
            let level = lod.resolve(structure.atom_count());
            let geom = build_ribbon_geometry(structure, &colors, &visible, visibility, level);
            for p in &geom {
                match p {
                    RibbonPrimitive::Cylinder {
                        a,
                        b,
                        r_world,
                        min_r,
                        color,
                    } => {
                        if let (Some(p1), Some(p2)) = (proj(*a), proj(*b)) {
                            let avg = (p1.2 + p2.2) * 0.5;
                            let w =
                                project_radius(*r_world, avg, camera.fov, height).max(*min_r) * 2.0;
                            let lit =
                                lighting.shade(cam_z, avg, *color, avg - *r_world, avg + *r_world);
                            prims.push(SvgPrim::Line {
                                x1: p1.0,
                                y1: p1.1,
                                x2: p2.0,
                                y2: p2.1,
                                w,
                                depth: avg,
                                color: lit,
                            });
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
                        let s0l = proj(*v0l);
                        let s0r = proj(*v0r);
                        let s1l = proj(*v1l);
                        let s1r = proj(*v1r);
                        if let (Some(a), Some(b), Some(c), Some(d)) = (s0l, s0r, s1l, s1r) {
                            let avg = (a.2 + b.2 + c.2 + d.2) * 0.25;
                            let lit = lighting.shade(*normal, avg, *color, avg - 0.5, avg + 0.5);
                            prims.push(SvgPrim::Polygon {
                                pts: vec![(a.0, a.1), (b.0, b.1), (d.0, d.1), (c.0, c.1)],
                                depth: avg,
                                color: lit,
                            });
                        }
                    }
                    RibbonPrimitive::Sphere {
                        c,
                        r_world,
                        min_r,
                        color,
                    } => {
                        if let Some((sx, sy, depth)) = proj(*c) {
                            let r = project_radius(*r_world, depth, camera.fov, height).max(*min_r);
                            let lit = lighting.shade(
                                cam_z,
                                depth,
                                *color,
                                depth - *r_world,
                                depth + *r_world,
                            );
                            prims.push(SvgPrim::Circle {
                                cx: sx,
                                cy: sy,
                                r,
                                depth,
                                color: lit,
                            });
                        }
                    }
                }
            }
        }
        RenderMode::Wireframe => {
            push_bond_lines(
                structure,
                &colors,
                &visible,
                &mats,
                &camera,
                &lighting,
                width,
                height,
                &mut prims,
                |_| true,
            );
        }
    }

    prims.sort_by(|a, b| {
        b.depth()
            .partial_cmp(&a.depth())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut s = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {} {}" width="{}" height="{}">
<rect width="{}" height="{}" fill="black"/>
"#,
        width, height, width, height, width, height
    );
    for p in &prims {
        s.push_str(&p.emit());
        s.push('\n');
    }
    s.push_str("</svg>");
    s
}

#[allow(clippy::too_many_arguments)]
fn push_bond_lines(
    structure: &Structure,
    colors: &[PixelColor],
    visible: &[bool],
    mats: &CameraMatrices,
    camera: &Camera,
    lighting: &Lighting,
    width: usize,
    height: usize,
    out: &mut Vec<SvgPrim>,
    keep: impl Fn(&crate::model::atom::Atom) -> bool,
) {
    let atoms = structure.atoms();
    for bond in structure.bonds() {
        if bond.atom1_idx >= atoms.len() || bond.atom2_idx >= atoms.len() {
            continue;
        }
        let a1 = &atoms[bond.atom1_idx];
        let a2 = &atoms[bond.atom2_idx];
        if !keep(a1) && !keep(a2) {
            continue;
        }
        if !visible[a1.index] || !visible[a2.index] {
            continue;
        }
        if let (Some(p1), Some(p2)) = (
            camera.project(mats, a1.pos, width, height),
            camera.project(mats, a2.pos, width, height),
        ) {
            let avg = (p1.2 + p2.2) * 0.5;
            let w = project_radius(0.18, avg, camera.fov, height).max(0.5) * 2.0;
            let c1 = lighting.shade(
                Vec3::new(0.0, 0.0, 1.0),
                avg,
                colors[a1.index],
                avg - 0.18,
                avg + 0.18,
            );
            let c2 = lighting.shade(
                Vec3::new(0.0, 0.0, 1.0),
                avg,
                colors[a2.index],
                avg - 0.18,
                avg + 0.18,
            );
            let mid = ((p1.0 + p2.0) * 0.5, (p1.1 + p2.1) * 0.5);
            out.push(SvgPrim::Line {
                x1: p1.0,
                y1: p1.1,
                x2: mid.0,
                y2: mid.1,
                w,
                depth: avg,
                color: c1,
            });
            out.push(SvgPrim::Line {
                x1: mid.0,
                y1: mid.1,
                x2: p2.0,
                y2: p2.1,
                w,
                depth: avg,
                color: c2,
            });
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn export_mp4(
    structure: &Structure,
    mode: RenderMode,
    color: ColorScheme,
    width: usize,
    height: usize,
    ssaa: usize,
    frames: u32,
    fps: u32,
    visibility: Visibility,
    lod: LodMode,
    out_path: &str,
) -> Result<()> {
    let probe = Command::new("ffmpeg").arg("-version").output();
    let ok = matches!(probe, Ok(o) if o.status.success());
    if !ok {
        return Err(TermPdbError::Other(
            "ffmpeg not found in PATH; install ffmpeg to encode MP4 (or use --export-png for a single frame)".into(),
        ));
    }

    let (width, height, ssaa) = fit_render_size(width, height, ssaa);
    let ssaa = ssaa.max(1);
    let sw = width * ssaa;
    let sh = height * ssaa;
    let tmp = std::env::temp_dir().join(format!("termpdb_mp4_{}", std::process::id()));
    std::fs::create_dir_all(&tmp)?;
    let lighting = Lighting::default();
    let mut camera = fit_camera(structure, sw, sh);
    let delta = 2.0 * std::f32::consts::TAU / frames as f32;

    for i in 0..frames {
        let mut fb = Framebuffer::new(sw, sh);
        fb.clear((0, 0, 0));
        render_structure(
            structure, mode, color, &camera, &mut fb, &lighting, visibility, lod,
        );
        let rgba = downsample_rgba(&fb, width, height, ssaa);
        let p = tmp.join(format!("frame_{:06}.png", i));
        write_png(p.to_str().unwrap(), &rgba, width as u32, height as u32)?;
        camera.orbit_angle(delta);
    }

    let status = Command::new("ffmpeg")
        .args([
            "-y",
            "-framerate",
            &fps.to_string(),
            "-i",
            tmp.join("frame_%06d.png").to_str().unwrap(),
            "-pix_fmt",
            "yuv420p",
            out_path,
        ])
        .status()?;
    let _ = std::fs::remove_dir_all(&tmp);
    if !status.success() {
        return Err(TermPdbError::Other(format!(
            "ffmpeg exited with status {:?}",
            status.code()
        )));
    }
    Ok(())
}
