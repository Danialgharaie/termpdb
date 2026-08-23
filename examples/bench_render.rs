use std::time::Instant;
use termpdb::math::Vec3;
use termpdb::model::{Atom, Chain, Residue, Structure, element_by_symbol};
use termpdb::render::{
    Camera, ColorScheme, Framebuffer, Lighting, LodMode, RenderContext, RenderMode,
    RibbonPrimitive, Visibility, build_render_cache, build_ribbon_geometry, render_structure,
    render_structure_ctx,
};

fn build(n_per_chain: u32, chain_ids: &[&str]) -> Structure {
    let mut structure = Structure::new("bench");
    let c = element_by_symbol("C");
    let o = element_by_symbol("O");
    let n = element_by_symbol("N");
    for (ci, cid) in chain_ids.iter().copied().enumerate() {
        let mut chain = Chain::new(cid);
        for i in 0..n_per_chain {
            let mut res = Residue::new((i + 1) as i32, "ALA", cid);
            let x = (i % 80) as f32 * 1.6 + ci as f32 * 160.0;
            let y = ((i / 80) % 80) as f32 * 1.6;
            let z = (i / 6400) as f32 * 1.6;
            let elem = [c, o, n][i as usize % 3];
            let atom = Atom::new(
                0,
                (i + 1) as i32,
                "CA",
                elem,
                Vec3::new(x, y, z),
                20.0,
                "ALA",
                (i + 1) as i32,
                cid,
                false,
            );
            let idx = structure.add_atom(atom);
            res.atom_indices.push(idx);
            chain.residues.push(res);
        }
        structure.add_chain(chain);
    }
    structure.ensure_bonds();
    structure
}

fn fit_camera(structure: &Structure, w: usize, h: usize, zoom: f32) -> Camera {
    let mut camera = Camera::new();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    camera.fit_structure(com, radius);
    camera.aspect = (w as f32) / (h as f32 * 2.0);
    if zoom != 1.0 {
        camera.distance /= zoom;
        camera.near = 0.1;
        camera.far = camera.distance + radius * 3.0;
    }
    camera
}

fn count_drawn(fb: &Framebuffer) -> usize {
    fb.pixels.iter().filter(|&&p| p != (0, 0, 0)).count()
}

/// Per-frame cost rebuilding the cache each call (one-shot / export path).
fn bench_fresh(
    structure: &Structure,
    mode: RenderMode,
    scheme: ColorScheme,
    lod: LodMode,
    w: usize,
    h: usize,
    zoom: f32,
) -> (std::time::Duration, usize) {
    let camera = fit_camera(structure, w, h, zoom);
    let lighting = Lighting::default();
    let mut fb = Framebuffer::new(w, h * 2);
    fb.clear((0, 0, 0));
    render_structure(
        structure,
        mode,
        scheme,
        &camera,
        &mut fb,
        &lighting,
        Visibility::default(),
        lod,
    );
    let iters = 30;
    let start = Instant::now();
    for _ in 0..iters {
        fb.clear((0, 0, 0));
        render_structure(
            structure,
            mode,
            scheme,
            &camera,
            &mut fb,
            &lighting,
            Visibility::default(),
            lod,
        );
    }
    (start.elapsed() / iters as u32, count_drawn(&fb))
}

/// Per-frame cost with the cross-frame cache built ONCE then reused across frames
/// (mirrors the interactive orbit/spin loop, where the camera moves but the
/// structure/color/visibility/LOD do not).
fn bench_cached(
    structure: &Structure,
    mode: RenderMode,
    scheme: ColorScheme,
    lod: LodMode,
    w: usize,
    h: usize,
    zoom: f32,
) -> (std::time::Duration, usize) {
    let camera = fit_camera(structure, w, h, zoom);
    let lighting = Lighting::default();
    let visibility = Visibility::default();
    let (colors, visible, com, radius, max_vdw) =
        build_render_cache(structure, scheme, visibility, lod);
    let level = lod.resolve(structure.atom_count());
    let ribbon_geometry: Vec<RibbonPrimitive> = if mode == RenderMode::Ribbon {
        build_ribbon_geometry(structure, &colors, &visible, visibility, level)
    } else {
        Vec::new()
    };
    let mats = camera.matrices();
    let ctx = RenderContext {
        structure,
        camera: &camera,
        mats,
        lighting: &lighting,
        visibility,
        lod: level,
        colors: &colors,
        visible: &visible,
        com,
        radius,
        max_vdw,
        ribbon_geometry: Some(&ribbon_geometry),
    };
    let mut fb = Framebuffer::new(w, h * 2);
    fb.clear((0, 0, 0));
    render_structure_ctx(&ctx, mode, &mut fb);
    let iters = 30;
    let start = Instant::now();
    for _ in 0..iters {
        fb.clear((0, 0, 0));
        render_structure_ctx(&ctx, mode, &mut fb);
    }
    (start.elapsed() / iters as u32, count_drawn(&fb))
}

fn ms(d: std::time::Duration) -> f32 {
    d.as_secs_f32() * 1000.0
}
fn fps(d: std::time::Duration) -> f32 {
    1.0 / d.as_secs_f32().max(1e-9)
}

fn main() {
    for &n in &[5_000u32, 20_000u32] {
        let chains = if n <= 5_000 {
            vec!["A"]
        } else {
            vec!["A", "B", "C", "D"]
        };
        let structure = build(n, &chains);
        println!(
            "\n=== {} atoms (CPK): fresh vs cached (orbit/spin) ===",
            structure.atom_count()
        );
        for &mode in &[
            RenderMode::Vdw,
            RenderMode::BallAndStick,
            RenderMode::Ribbon,
            RenderMode::Trace,
        ] {
            let (df, drawn) = bench_fresh(
                &structure,
                mode,
                ColorScheme::Cpk,
                LodMode::Full,
                120,
                60,
                1.0,
            );
            let (dc, _) = bench_cached(
                &structure,
                mode,
                ColorScheme::Cpk,
                LodMode::Full,
                120,
                60,
                1.0,
            );
            println!(
                "{:13}  fresh {:7.3} ms ({:5.0} fps) | cached {:7.3} ms ({:5.0} fps)  {:.1}x  drawn={}",
                format!("{:?}", mode),
                ms(df),
                fps(df),
                ms(dc),
                fps(dc),
                ms(df) / ms(dc).max(1e-9),
                drawn
            );
        }
    }

    println!("\n=== BFactor coloring (was O(n^2)) ===");
    for &n in &[5_000u32, 20_000u32] {
        let chains = if n <= 5_000 {
            vec!["A"]
        } else {
            vec!["A", "B", "C", "D"]
        };
        let structure = build(n, &chains);
        let (d, drawn) = bench_cached(
            &structure,
            RenderMode::Vdw,
            ColorScheme::BFactor,
            LodMode::Full,
            120,
            60,
            1.0,
        );
        println!(
            "  {:>6} atoms VDW/BFactor  {:7.3} ms/frame  ({:5.0} fps)  drawn={}",
            structure.atom_count(),
            ms(d),
            fps(d),
            drawn
        );
    }

    println!("\n=== Zoomed-in VDW (heavy overdraw, 1000 atoms, 160x80) ===");
    let structure = build(1_000, &["A"]);
    for &(label, zoom) in &[("1x", 1.0f32), ("4x", 4.0), ("10x", 10.0)] {
        let (d, drawn) = bench_cached(
            &structure,
            RenderMode::Vdw,
            ColorScheme::Cpk,
            LodMode::Full,
            160,
            80,
            zoom,
        );
        println!(
            "  zoom {:3}  {:7.3} ms/frame  ({:5.0} fps)  drawn={}/{}",
            label,
            ms(d),
            fps(d),
            drawn,
            160 * 80 * 2
        );
    }
}
