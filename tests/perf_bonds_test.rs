//! Parity and determinism tests for the band-parallel bond/cylinder passes.
//!
//! The bond passes (Ball & Stick bonds, Trace tubes, Wireframe lines) collect
//! drawable primitives in a single-threaded projection pass and then rasterize
//! them either serially ([`draw_band_primitives_serial`]) or split across
//! horizontal row bands ([`draw_band_primitives_parallel`]). These tests pin
//! the invariant that both dispatchers -- and therefore the pre-parallel
//! inline loops they replaced -- produce bit-identical framebuffers.

use termpdb::model::Bond;
use termpdb::model::Structure;
use termpdb::model::bond::BondDetector;
use termpdb::render::buffer::{Framebuffer, PixelColor};
use termpdb::render::camera::Camera;
use termpdb::render::lighting::Lighting;
use termpdb::render::rasterizer::{draw_line_3d, draw_sphere};
use termpdb::render::representations::{
    BandPrimitive, LodLevel, LodMode, RenderContext, RenderMode, Visibility, build_render_cache,
    draw_band_primitives_parallel, draw_band_primitives_serial, project_radius,
    render_structure_ctx,
};
use termpdb::render::{ColorScheme, color_for_atom};

/// Builds a PDB string describing an `nx * ny * nz` cubic lattice of atoms
/// spaced 1.5 A apart with cyclically varying elements.
///
/// Consecutive atoms are ~1.5 A apart, well inside every covalent-bond
/// cutoff, so BondDetector yields roughly `3 * n` bonds; neighboring atoms
/// usually differ in element, so most bonds are drawn as two half-colored
/// cylinders. Every atom gets its own residue (named/guided like a CA) so
/// Trace mode connects them into a long backbone chain.
fn lattice_pdb(nx: usize, ny: usize, nz: usize) -> String {
    const ELEMENTS: [&str; 4] = ["C", "N", "O", "S"];
    let mut pdb = String::from("HEADER    LATTICE\n");
    let mut serial = 1usize;
    let mut resid = 0usize;
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                resid += 1;
                let x = i as f32 * 1.5;
                let y = j as f32 * 1.5;
                let z = k as f32 * 1.5;
                let elem = ELEMENTS[(i + j + k) % ELEMENTS.len()];
                pdb.push_str(&format!(
                    "ATOM  {serial:>5}  CA  LIG A{resid:>4}    {x:8.3}{y:8.3}{z:8.3}  1.00  0.00          {elem:>2}\n"
                ));
                serial += 1;
            }
        }
    }
    pdb.push_str("END\n");
    pdb
}

/// Owns everything a [`RenderContext`] borrows, so tests can build contexts
/// without fighting lifetimes.
struct Scene {
    structure: Structure,
    camera: Camera,
    lighting: Lighting,
    colors: Vec<PixelColor>,
    visible: Vec<bool>,
}

impl Scene {
    /// Parses `pdb`, fits the camera to the structure, and builds the
    /// camera-independent render cache.
    fn new(pdb: &str, width: usize, height: usize) -> Self {
        let structure = termpdb::parser::parse_pdb(pdb).expect("lattice PDB must parse");
        let (colors, visible, com, radius, _max_vdw) = build_render_cache(
            &structure,
            ColorScheme::Cpk,
            Visibility::default(),
            LodMode::Full,
        );

        let mut camera = Camera::new();
        camera.aspect = width as f32 / height as f32;
        camera.fit_structure(com, radius);

        Self {
            structure,
            camera,
            lighting: Lighting::default(),
            colors,
            visible,
        }
    }

    fn ctx(&self) -> RenderContext<'_> {
        let com = self.structure.center_of_mass();
        let radius = self.structure.bounding_sphere_radius();
        RenderContext {
            structure: &self.structure,
            camera: &self.camera,
            mats: self.camera.matrices(),
            lighting: &self.lighting,
            visibility: Visibility::default(),
            lod: LodLevel::Full,
            colors: &self.colors,
            visible: &self.visible,
            com,
            radius,
            max_vdw: 2.0,
            ribbon_geometry: None,
        }
    }
}

fn assert_framebuffers_identical(a: &Framebuffer, b: &Framebuffer) {
    assert_eq!(a.width, b.width);
    assert_eq!(a.height, b.height);
    assert_eq!(a.pixels, b.pixels, "color buffers diverged");
    assert_eq!(a.depth, b.depth, "depth buffers diverged");
}

/// Number of drawable bond/cylinder primitives the Ball & Stick pass would
/// emit for `scene` (projection pass replica used to size the assertions).
fn count_drawn_atoms(scene: &Scene, width: usize, height: usize) -> usize {
    let ctx = scene.ctx();
    ctx.structure
        .atoms()
        .iter()
        .filter(|a| ctx.visible[a.index])
        .filter(|a| {
            ctx.camera
                .project_sphere(&ctx.mats, a.pos, 0.38, width, height)
                .is_some_and(|(_, r)| r >= 0.4)
        })
        .count()
}

#[test]
fn test_lattice_has_many_split_color_bonds() {
    // Sanity guard for the fixtures below: enough atoms, plenty of bonds, and
    // genuinely mixed colors so split-cylinder halves dominate.
    let scene = Scene::new(&lattice_pdb(8, 8, 4), 160, 120);
    assert_eq!(scene.structure.atoms().len(), 256);
    assert_eq!(scene.colors.len(), 256);

    let bonds = BondDetector::detect_bonds(scene.structure.atoms());
    assert!(
        bonds.len() > 300,
        "expected a dense bond graph, got {}",
        bonds.len()
    );

    let mixed = bonds
        .iter()
        .filter(|b| scene.colors[b.atom1_idx] != scene.colors[b.atom2_idx])
        .count();
    assert!(
        mixed > 250,
        "expected mostly split-color bonds, got {mixed}"
    );
}

#[test]
fn test_ball_stick_parallel_deterministic() {
    let scene = Scene::new(&lattice_pdb(8, 8, 4), 160, 120);
    let ctx = scene.ctx();
    assert!(count_drawn_atoms(&scene, 160, 120) > 50);

    let mut run1 = Framebuffer::new(160, 120);
    run1.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::BallAndStick, &mut run1);

    let mut run2 = Framebuffer::new(160, 120);
    run2.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::BallAndStick, &mut run2);

    assert_framebuffers_identical(&run1, &run2);
    let drawn = run1.pixels.iter().filter(|&&p| p != (0, 0, 0)).count();
    assert!(drawn > 200, "ball-and-stick should draw plenty of pixels");
}

#[test]
fn test_trace_parallel_deterministic() {
    let scene = Scene::new(&lattice_pdb(8, 8, 4), 160, 120);
    let ctx = scene.ctx();

    let mut run1 = Framebuffer::new(160, 120);
    run1.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::Trace, &mut run1);

    let mut run2 = Framebuffer::new(160, 120);
    run2.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::Trace, &mut run2);

    assert_framebuffers_identical(&run1, &run2);
    let drawn = run1.pixels.iter().filter(|&&p| p != (0, 0, 0)).count();
    assert!(drawn > 200, "trace should draw plenty of pixels");
}

#[test]
fn test_wireframe_parallel_deterministic() {
    let scene = Scene::new(&lattice_pdb(8, 8, 4), 160, 120);
    let ctx = scene.ctx();

    let mut run1 = Framebuffer::new(160, 120);
    run1.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::Wireframe, &mut run1);

    let mut run2 = Framebuffer::new(160, 120);
    run2.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::Wireframe, &mut run2);

    assert_framebuffers_identical(&run1, &run2);
    let drawn = run1.pixels.iter().filter(|&&p| p != (0, 0, 0)).count();
    assert!(drawn > 200, "wireframe should draw plenty of pixels");
}

/// Fully serial Ball & Stick rendering: the preserved pre-parallel behavior
/// (inline sphere loop, then inline bond loop drawing whole-frame cylinders).
///
/// This is the reference the production parallel path must match pixel for
/// pixel. Keep the rasterization formulas in lockstep with
/// `representations::ball_stick::render_ball_stick`.
fn render_ball_stick_serial_reference(scene: &Scene, buffer: &mut Framebuffer) {
    let ctx = scene.ctx();
    let lighting = &scene.lighting;

    // 1. Atoms as analytical spheres, drawn one by one.
    let mut projected = Vec::new();
    for atom in ctx.structure.atoms() {
        if !ctx.visible[atom.index] {
            continue;
        }
        if let Some((pt, sphere_r)) =
            ctx.camera
                .project_sphere(&ctx.mats, atom.pos, 0.38, buffer.width, buffer.height)
        {
            if sphere_r < 0.4 {
                continue;
            }
            let sphere_r = if ctx.lod == LodLevel::Full {
                sphere_r.max(0.8)
            } else {
                sphere_r
            };
            projected.push((pt, sphere_r, ctx.colors[atom.index]));
        }
    }
    for (pt, r, color) in projected {
        draw_sphere(buffer, pt, r, color, lighting);
    }

    // 2. Bonds as whole-frame cylinders (split-color bonds as two halves).
    let detected;
    let bonds: &[Bond] = if !ctx.structure.bonds().is_empty() {
        ctx.structure.bonds()
    } else {
        detected = BondDetector::detect_bonds(ctx.structure.atoms());
        &detected
    };
    let mut prims: Vec<BandPrimitive> = Vec::with_capacity(bonds.len());
    for bond in bonds {
        let (Some(a1), Some(a2)) = (
            ctx.structure.atoms().get(bond.atom1_idx),
            ctx.structure.atoms().get(bond.atom2_idx),
        ) else {
            continue;
        };
        if !ctx.visible[a1.index] || !ctx.visible[a2.index] {
            continue;
        }
        let Some((p1, p2)) =
            ctx.camera
                .project_segment(&ctx.mats, a1.pos, a2.pos, buffer.width, buffer.height)
        else {
            continue;
        };
        let avg_depth = (p1.2 + p2.2) * 0.5;
        let bond_r = project_radius(0.18, avg_depth, ctx.camera.fov, buffer.height).max(0.5);
        let c1 = ctx.colors[a1.index];
        let c2 = ctx.colors[a2.index];
        if c1 == c2 {
            prims.push(BandPrimitive::Cylinder {
                p1,
                p2,
                radius: bond_r,
                color: c1,
            });
        } else {
            let pmid = ((p1.0 + p2.0) * 0.5, (p1.1 + p2.1) * 0.5, avg_depth);
            prims.push(BandPrimitive::Cylinder {
                p1,
                p2: pmid,
                radius: bond_r,
                color: c1,
            });
            prims.push(BandPrimitive::Cylinder {
                p1: pmid,
                p2,
                radius: bond_r,
                color: c2,
            });
        }
    }
    draw_band_primitives_serial(buffer, &prims, lighting);
}

#[test]
fn test_ball_stick_parallel_matches_serial_reference() {
    // End-to-end parity: production pipeline (parallel sphere pass + parallel
    // bond pass) vs. the fully serial reference on the same scene.
    let scene = Scene::new(&lattice_pdb(8, 8, 4), 160, 120);
    let ctx = scene.ctx();

    let mut par = Framebuffer::new(160, 120);
    par.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::BallAndStick, &mut par);

    let mut ser = Framebuffer::new(160, 120);
    ser.clear((0, 0, 0));
    render_ball_stick_serial_reference(&scene, &mut ser);

    assert_framebuffers_identical(&par, &ser);
}

#[test]
fn test_wireframe_parallel_matches_serial_lines() {
    // Wireframe parity: production path (collected lines, band dispatcher) vs.
    // the old inline loop calling draw_line_3d per bond.
    let scene = Scene::new(&lattice_pdb(8, 8, 4), 160, 120);
    let ctx = scene.ctx();

    let mut par = Framebuffer::new(160, 120);
    par.clear((0, 0, 0));
    render_structure_ctx(&ctx, RenderMode::Wireframe, &mut par);

    let mut ser = Framebuffer::new(160, 120);
    ser.clear((0, 0, 0));
    let atoms = ctx.structure.atoms();
    let detected;
    let bonds: &[Bond] = if !ctx.structure.bonds().is_empty() {
        ctx.structure.bonds()
    } else {
        detected = BondDetector::detect_bonds(atoms);
        &detected
    };
    for bond in bonds {
        let (Some(a1), Some(a2)) = (atoms.get(bond.atom1_idx), atoms.get(bond.atom2_idx)) else {
            continue;
        };
        if !ctx.visible[a1.index] || !ctx.visible[a2.index] {
            continue;
        }
        if let Some((p1, p2)) = ctx
            .camera
            .project_segment(&ctx.mats, a1.pos, a2.pos, ser.width, ser.height)
        {
            draw_line_3d(&mut ser, p1, p2, ctx.colors[a1.index]);
        }
    }

    assert_framebuffers_identical(&par, &ser);
}

#[test]
fn test_band_primitive_dispatcher_parity_edge_cases() {
    // Primitive-layer parity over a torture list: exact 0.5-radius fallback
    // boundary, zero-length capsules, split pairs sharing a midpoint,
    // coincident order-sensitive duplicates, off-screen and far-out-of-bounds
    // segments (screen clipping), band-spanning diagonals, and spheres
    // interleaved with cylinders. Both dispatchers must agree bit-for-bit.
    let red: PixelColor = (220, 40, 40);
    let green: PixelColor = (40, 220, 40);
    let blue: PixelColor = (60, 90, 240);
    let yellow: PixelColor = (230, 210, 40);

    let cyl = |p1: (f32, f32, f32), p2: (f32, f32, f32), radius: f32, color: PixelColor| {
        BandPrimitive::Cylinder {
            p1,
            p2,
            radius,
            color,
        }
    };
    let sph = |center: (f32, f32, f32), radius: f32, color: PixelColor| BandPrimitive::Sphere {
        center,
        radius,
        color,
    };

    let mid = ((10.0 + 70.0) * 0.5, (5.0 + 55.0) * 0.5, (3.0 + 9.0) * 0.5);
    let prims: Vec<BandPrimitive> = vec![
        // Zero-length capsule (degenerate axis, len_sq < 1e-6 branch).
        cyl((40.0, 30.0, 10.0), (40.0, 30.0, 10.0), 6.0, red),
        // Exact fallback threshold and its two sides.
        cyl((8.0, 4.0, 5.0), (72.0, 44.0, 20.0), 0.5, green),
        cyl((8.0, 8.0, 5.0), (72.0, 48.0, 20.0), 0.49, blue),
        cyl((8.0, 12.0, 5.0), (72.0, 52.0, 20.0), 0.51, yellow),
        // Steep diagonal crossing every band, split into two colored halves
        // sharing pmid (exactly what split-color bonds emit).
        cyl((10.0, 5.0, 3.0), mid, 2.5, red),
        cyl(mid, (70.0, 55.0, 9.0), 2.5, green),
        // Coincident duplicates: depth ties must resolve by draw order.
        cyl((20.0, 20.0, 6.0), (60.0, 20.0, 6.0), 4.0, red),
        cyl((20.0, 20.0, 6.0), (60.0, 20.0, 6.0), 4.0, blue),
        // Mostly off-screen segment clipped back into the frame.
        cyl((-500.0, -400.0, 2.0), (600.0, 500.0, 40.0), 3.0, green),
        // Entirely above / below / beside the viewport.
        cyl((10.0, -30.0, 5.0), (60.0, -10.0, 5.0), 2.0, blue),
        cyl((10.0, 90.0, 5.0), (60.0, 130.0, 5.0), 2.0, blue),
        cyl((-80.0, 10.0, 5.0), (-20.0, 50.0, 5.0), 2.0, blue),
        // Near-plane pathology: endpoints millions of pixels out; the thick
        // variant stresses bbox clamping, the thin one the line clip + step cap.
        cyl(
            (-1_000_000.0, -1_000_000.0, 3.0),
            (3_000_000.0, 2_000_000.0, 7.0),
            8.0,
            yellow,
        ),
        cyl(
            (-1_000_000.0, -1_000_000.0, 3.0),
            (3_000_000.0, 2_000_000.0, 7.0),
            0.5,
            yellow,
        ),
        // Spheres interleaved among the cylinders (trace-style ordering).
        sph((35.0, 29.5, 8.0), 5.0, blue),
        sph((45.0, 45.0, 12.0), 2.0, red),
        sph((0.0, 0.0, 4.0), 3.0, green),
        sph((79.0, 59.0, 4.0), 3.0, green),
    ];

    let mut serial = Framebuffer::new(80, 60);
    serial.clear((0, 0, 0));
    let lighting = Lighting::default();
    draw_band_primitives_serial(&mut serial, &prims, &lighting);

    let mut parallel = Framebuffer::new(80, 60);
    parallel.clear((0, 0, 0));
    draw_band_primitives_parallel(&mut parallel, &prims, &lighting);

    assert_framebuffers_identical(&serial, &parallel);
    let drawn = parallel.pixels.iter().filter(|&&p| p != (0, 0, 0)).count();
    assert!(drawn > 500, "torture list should cover the frame densely");
}

#[test]
fn test_band_primitive_parity_various_band_counts() {
    // Parity must hold for any band partition, including a single band and
    // ragged final bands (heights not divisible by the 16-row band size).
    let prims: Vec<BandPrimitive> = (0..24)
        .map(|i| {
            let t = i as f32;
            BandPrimitive::Cylinder {
                p1: (2.0 + t, 2.0 + t * 2.0, 4.0 + t),
                p2: (58.0 + t, 56.0 - t * 2.0, 12.0 + t * 0.5),
                radius: 0.4 + t * 0.25,
                color: (i as u8 * 10, 255 - i as u8 * 9, i as u8 * 7),
            }
        })
        .collect();

    let lighting = Lighting::default();
    for &(w, h) in &[(40u16, 8u16), (64, 17), (80, 31), (100, 47), (33, 64)] {
        let mut serial = Framebuffer::new(w as usize, h as usize);
        serial.clear((10, 10, 10));
        draw_band_primitives_serial(&mut serial, &prims, &lighting);

        let mut parallel = Framebuffer::new(w as usize, h as usize);
        parallel.clear((10, 10, 10));
        draw_band_primitives_parallel(&mut parallel, &prims, &lighting);

        assert_framebuffers_identical(&serial, &parallel);
    }
}

#[test]
fn test_colors_follow_cpk_scheme_for_alternating_elements() {
    // Guards the fixture: the explicit element column drives CPK colors, so
    // neighboring lattice atoms really do render as split-color bonds.
    let scene = Scene::new(&lattice_pdb(2, 1, 1), 80, 60);
    let atoms = scene.structure.atoms();
    let expected_c = color_for_atom(&atoms[0], None, &scene.structure, ColorScheme::Cpk);
    assert_eq!(scene.colors[0], expected_c);
    assert_ne!(scene.colors[0], scene.colors[1]);
}
