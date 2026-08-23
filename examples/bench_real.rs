use std::time::Instant;
use termpdb::parser::parse_pdb;
use termpdb::render::{
    Camera, ColorScheme, Framebuffer, Lighting, LodMode, RenderContext, RenderMode, Visibility,
    build_render_cache, build_ribbon_geometry, render_structure_ctx,
};

fn ms(d: std::time::Duration) -> f32 { d.as_secs_f32() * 1000.0 }

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| "/tmp/4egk.pdb".to_string());
    let pdb = std::fs::read_to_string(&path).expect("read pdb");

    let t = Instant::now();
    let mut structure = parse_pdb(&pdb).expect("parse");
    println!("parse:              {:7.2} ms", ms(t.elapsed()));

    let t = Instant::now();
    structure.ensure_bonds();
    println!("ensure_bonds:       {:7.2} ms  (bonds={})", ms(t.elapsed()), structure.bonds().len());

    println!("atoms={} residues={} chains={} models={}",
        structure.atom_count(), structure.residue_count(), structure.chain_count(), structure.model_count());

    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    let mut camera = Camera::new();
    camera.fit_structure(com, radius);
    camera.aspect = 1.0;
    let lighting = Lighting::default();
    let visibility = Visibility::default();

    let t = Instant::now();
    let (colors, visible, com, radius, max_vdw) =
        build_render_cache(&structure, ColorScheme::Cpk, visibility, LodMode::Full);
    println!("build_render_cache: {:7.2} ms", ms(t.elapsed()));

    let level = LodMode::Full.resolve(structure.atom_count());
    let t = Instant::now();
    let rg = build_ribbon_geometry(&structure, &colors, &visible, visibility, level);
    println!("build_ribbon_geom: {:7.2} ms  (prims={})", ms(t.elapsed()), rg.len());

    let mats = camera.matrices();
    for &mode in &[RenderMode::Vdw, RenderMode::BallAndStick, RenderMode::Ribbon, RenderMode::Trace] {
        let rg_ref = if mode == RenderMode::Ribbon { Some(&rg[..]) } else { None };
        let ctx = RenderContext {
            structure: &structure, camera: &camera, mats, lighting: &lighting, visibility,
            lod: level, colors: &colors, visible: &visible, com, radius, max_vdw,
            ribbon_geometry: rg_ref,
        };
        let mut fb = Framebuffer::new(160, 100);
        fb.clear((0, 0, 0));
        render_structure_ctx(&ctx, mode, &mut fb);
        let iters = 200;
        let t = Instant::now();
        for _ in 0..iters {
            fb.clear((0, 0, 0));
            render_structure_ctx(&ctx, mode, &mut fb);
        }
        let d = t.elapsed() / iters;
        println!("{:13} {:7.3} ms/frame  ({:6.1} fps)", format!("{:?}", mode), ms(d), 1.0 / d.as_secs_f32().max(1e-9));
    }

    println!("\n-- with BFactor color scheme --");
    let (colors, visible, com, radius, max_vdw) =
        build_render_cache(&structure, ColorScheme::BFactor, visibility, LodMode::Full);
    let ctx = RenderContext {
        structure: &structure, camera: &camera, mats, lighting: &lighting, visibility,
        lod: level, colors: &colors, visible: &visible, com, radius, max_vdw,
        ribbon_geometry: None,
    };
    let mut fb = Framebuffer::new(160, 100);
    let iters = 200;
    let t = Instant::now();
    for _ in 0..iters { fb.clear((0,0,0)); render_structure_ctx(&ctx, RenderMode::Vdw, &mut fb); }
    let d = t.elapsed()/iters;
    println!("VDW/BFactor       {:7.3} ms/frame  ({:6.1} fps)", ms(d), 1.0/d.as_secs_f32().max(1e-9));
}
