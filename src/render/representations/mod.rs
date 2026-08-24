//! 3D Molecular structural representations (Trace, Ball & Stick, Ribbon, VDW).
//!
//! Generates 3D geometric visual primitives from molecular data models:
//! - **Trace**: C-alpha / backbone line/cylinder trace with ligand spheres.
//! - **Ball & Stick**: All-atom analytical spheres and covalent bond cylinders.
//! - **Ribbon**: Secondary structure cartoon (helices as thick tubes, sheets as arrow ribbons, coils as smooth tubes).
//! - **VDW**: Van der Waals space-filling analytical spheres.

pub mod ball_stick;
pub mod nucleic;
pub mod ribbon;
pub mod trace;
pub mod vdw;

use std::collections::HashMap;

pub use ball_stick::render_ball_stick;
pub use ribbon::{RibbonPrimitive, build_ribbon_geometry, render_ribbon};
pub use trace::render_trace;
pub use vdw::render_vdw;

use crate::math::Vec3;
use crate::model::{Atom, Residue, Structure};
use crate::render::buffer::{Framebuffer, PixelColor};
use crate::render::camera::{Camera, CameraMatrices};
use crate::render::color::{ColorScheme, ColorStats, color_for_atom_with_stats};
use crate::render::lighting::Lighting;
use crate::render::rasterizer::{
    ScreenPoint, draw_cylinder, draw_cylinder_band, draw_line_3d, draw_line_3d_band, draw_sphere,
    draw_sphere_band,
};

/// Available molecular representation rendering modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, clap::ValueEnum)]
pub enum RenderMode {
    /// Backbone CA / Nucleic phosphorus trace
    #[value(name = "trace")]
    Trace,
    /// All-atom ball and stick
    #[value(name = "ball-and-stick", alias = "ball_and_stick")]
    BallAndStick,
    /// Secondary structure cartoon ribbon
    #[default]
    #[value(name = "ribbon")]
    Ribbon,
    /// Space-filling Van der Waals spheres
    #[value(name = "vdw")]
    Vdw,
    /// 3D Wireframe / line representation
    #[value(name = "wireframe", alias = "lines")]
    Wireframe,
}

impl RenderMode {
    /// Returns the human-readable display name of the render mode.
    pub fn name(&self) -> &'static str {
        match self {
            RenderMode::Trace => "Trace",
            RenderMode::BallAndStick => "Ball & Stick",
            RenderMode::Ribbon => "Ribbon",
            RenderMode::Vdw => "VDW",
            RenderMode::Wireframe => "Wireframe",
        }
    }

    /// Returns an array of all available render modes in cycle order.
    pub fn all() -> &'static [RenderMode] {
        &[
            RenderMode::Trace,
            RenderMode::BallAndStick,
            RenderMode::Ribbon,
            RenderMode::Vdw,
            RenderMode::Wireframe,
        ]
    }

    /// Cycles to the next render mode.
    pub fn next(&self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&m| m == *self).unwrap_or(0);
        all[(idx + 1) % all.len()]
    }

    /// Cycles to the previous render mode.
    pub fn prev(&self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&m| m == *self).unwrap_or(0);
        all[(idx + all.len() - 1) % all.len()]
    }
}

/// Computes the projected screen pixel radius from a 3D world radius and view depth.
///
/// Delegates to the canonical implementation in [`crate::render::camera`].
pub fn project_radius(world_radius: f32, view_depth: f32, fov: f32, height: usize) -> f32 {
    crate::render::camera::project_radius(world_radius, view_depth, fov, height)
}

/// One drawable screen-space primitive collected by the single-threaded
/// projection pass of a bond/cylinder-dominated representation (Ball & Stick
/// bonds, Trace tubes and ligands, Wireframe lines).
///
/// Endpoints are exactly the values the serial rasterizers would receive; the
/// band-parallel path re-derives everything it needs per band, which is what
/// keeps its output pixel-identical to [`draw_band_primitives_serial`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BandPrimitive {
    /// Thick shaded capsule between two screen points. Radii `<= 0.5` px fall
    /// back to a 1-px line inside the rasterizer, exactly like `draw_cylinder`.
    Cylinder {
        p1: ScreenPoint,
        p2: ScreenPoint,
        radius: f32,
        color: PixelColor,
    },
    /// 1-px depth-tested line (wireframe bonds).
    Line {
        p1: ScreenPoint,
        p2: ScreenPoint,
        color: PixelColor,
    },
    /// Analytical shaded sphere (trace joints / ligand atoms).
    Sphere {
        center: ScreenPoint,
        radius: f32,
        color: PixelColor,
    },
}

impl BandPrimitive {
    /// Conservative screen-space vertical extent `(y_min, y_max)` of the
    /// primitive, used for O(1) per-band rejection in the parallel pass.
    ///
    /// The bounds are deliberately loose (lines pad by 1 px for Bresenham
    /// rounding): a band whose row range cannot contain a written pixel may be
    /// visited needlessly, but a band containing at least one writable pixel
    /// is never skipped. That one-sided guarantee is what makes the per-band
    /// cull safe for bit-exact parity with the serial pass.
    fn y_range(&self) -> (f32, f32) {
        match *self {
            BandPrimitive::Cylinder { p1, p2, radius, .. } => {
                (p1.1.min(p2.1) - radius, p1.1.max(p2.1) + radius)
            }
            BandPrimitive::Line { p1, p2, .. } => (p1.1.min(p2.1) - 1.0, p1.1.max(p2.1) + 1.0),
            BandPrimitive::Sphere { center, radius, .. } => (center.1 - radius, center.1 + radius),
        }
    }
}

/// Draws collected primitives single-threaded in collection order.
///
/// Reference path for the band-parallel renderer (and its fallback for tiny
/// frames / small batches): it replays each primitive through the ordinary
/// whole-frame rasterizers exactly like the pre-parallel inline bond loops did.
pub fn draw_band_primitives_serial(
    buffer: &mut Framebuffer,
    prims: &[BandPrimitive],
    lighting: &Lighting,
) {
    for &prim in prims {
        match prim {
            BandPrimitive::Cylinder {
                p1,
                p2,
                radius,
                color,
            } => {
                draw_cylinder(buffer, p1, p2, radius, color, lighting);
            }
            BandPrimitive::Line { p1, p2, color } => {
                draw_line_3d(buffer, p1, p2, color);
            }
            BandPrimitive::Sphere {
                center,
                radius,
                color,
            } => {
                draw_sphere(buffer, center, radius, color, lighting);
            }
        }
    }
}

/// Band-parallel twin of [`draw_band_primitives_serial`].
///
/// Splits the framebuffer into horizontal bands with `Framebuffer::par_bands_mut`,
/// then scatters every primitive into the bands its y-range intersects (a
/// single-threaded O(primitives) pass that keeps collection order inside each
/// bucket). Each band is drawn independently by rayon, visiting exactly the
/// primitives that can touch its rows, in the same order the serial pass used.
/// Because bands own disjoint rows, every pixel sees the same sequence of
/// depth tests as the serial pass -- the two produce identical framebuffers.
pub fn draw_band_primitives_parallel(
    buffer: &mut Framebuffer,
    prims: &[BandPrimitive],
    lighting: &Lighting,
) {
    if buffer.width == 0 || buffer.height == 0 || prims.is_empty() {
        return;
    }

    const BAND_HEIGHT: usize = 16;
    let full_height = buffer.height;
    let mut bands = buffer.par_bands_mut(BAND_HEIGHT);

    // Scatter pass: bucket primitives per band. Band b covers global rows
    // [b * BAND_HEIGHT, (b + 1) * BAND_HEIGHT); a primitive intersects it iff
    // y_max >= b * BAND_HEIGHT && y_min < (b + 1) * BAND_HEIGHT, which is
    // exactly the clamped band range floor(y_min / BAND_HEIGHT)
    // ..= floor(y_max / BAND_HEIGHT). Pushing in collection order keeps each
    // bucket in serial visitation order.
    let band_count = bands.len();
    let mut buckets: Vec<Vec<&BandPrimitive>> = vec![Vec::new(); band_count];
    for prim in prims {
        let (y_min, y_max) = prim.y_range();
        let first = ((y_min.max(0.0) as usize) / BAND_HEIGHT).min(band_count - 1);
        let last = ((y_max.max(0.0) as usize) / BAND_HEIGHT).min(band_count - 1);
        for bucket in &mut buckets[first..=last] {
            bucket.push(prim);
        }
    }

    use rayon::prelude::*;
    bands.par_iter_mut().enumerate().for_each(|(i, band)| {
        for prim in &buckets[i] {
            match **prim {
                BandPrimitive::Cylinder {
                    p1,
                    p2,
                    radius,
                    color,
                } => {
                    draw_cylinder_band(band, p1, p2, radius, color, lighting, full_height);
                }
                BandPrimitive::Line { p1, p2, color } => {
                    draw_line_3d_band(band, p1, p2, color, full_height);
                }
                BandPrimitive::Sphere {
                    center,
                    radius,
                    color,
                } => {
                    draw_sphere_band(band, center, radius, color, lighting);
                }
            }
        }
    });
}

/// Draws collected primitives, choosing the band-parallel path when the frame
/// spans multiple bands and there are enough primitives to amortize the split
/// (same heuristic as the sphere passes).
pub(crate) fn draw_band_primitives(
    buffer: &mut Framebuffer,
    prims: &[BandPrimitive],
    lighting: &Lighting,
) {
    // Same banding as `Framebuffer::par_bands_mut` with a 16-row band height.
    let band_count = if buffer.width == 0 {
        0
    } else {
        buffer.height.div_ceil(16)
    };
    if band_count > 1 && prims.len() > 50 {
        draw_band_primitives_parallel(buffer, prims, lighting);
    } else {
        draw_band_primitives_serial(buffer, prims, lighting);
    }
}

/// Linearly interpolates between two RGB pixel colors.
pub fn lerp_rgb(c1: PixelColor, c2: PixelColor, t: f32) -> PixelColor {
    let t = t.clamp(0.0, 1.0);
    let r = (c1.0 as f32 + t * (c2.0 as f32 - c1.0 as f32)).round() as u8;
    let g = (c1.1 as f32 + t * (c2.1 as f32 - c1.1 as f32)).round() as u8;
    let b = (c1.2 as f32 + t * (c2.2 as f32 - c1.2 as f32)).round() as u8;
    (r, g, b)
}

/// Which atoms to draw. This is a view filter, not a change to the structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Visibility {
    /// Draw water / solvent residues (HOH, WAT, …). Default: hidden.
    pub show_waters: bool,
    /// Draw hydrogen atoms. Default: shown.
    pub show_hydrogens: bool,
}

impl Default for Visibility {
    fn default() -> Self {
        Self {
            show_waters: false,
            show_hydrogens: true,
        }
    }
}

impl Visibility {
    /// Draw every atom, including waters and hydrogens.
    pub const ALL: Self = Self {
        show_waters: true,
        show_hydrogens: true,
    };

    /// Returns whether `atom` should be rasterized under this filter.
    pub fn atom_visible(&self, atom: &Atom, residue: Option<&Residue>) -> bool {
        if !self.show_hydrogens && atom.is_hydrogen() {
            return false;
        }
        let is_water = residue
            .map(|r| r.is_water())
            .unwrap_or_else(|| Residue::name_is_water(&atom.res_name));
        if !self.show_waters && is_water {
            return false;
        }
        true
    }

    /// Residue-level water skip for backbone traces and ribbons.
    pub fn residue_visible(&self, residue: &Residue) -> bool {
        self.show_waters || !residue.is_water()
    }
}

/// User-facing level-of-detail. `Auto` picks a level from atom count.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, clap::ValueEnum)]
pub enum LodMode {
    /// Full below 25k atoms, backbone below 80k, Cα/P above that.
    #[default]
    #[value(name = "auto")]
    Auto,
    /// Draw every visible atom.
    #[value(name = "full")]
    Full,
    /// Polymer backbone plus non-water heteroatoms.
    #[value(name = "backbone")]
    Backbone,
    /// Cα / nucleic P plus non-water heteroatoms.
    #[value(name = "ca", alias = "calpha")]
    CAlpha,
}

/// Resolved drawing detail (what Auto chose, or a locked mode).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LodLevel {
    Full,
    Backbone,
    CAlpha,
}

/// Atom-count thresholds for [`LodMode::Auto`].
pub const LOD_BACKBONE_ATOMS: usize = 25_000;
pub const LOD_CALPHA_ATOMS: usize = 80_000;

impl LodMode {
    pub fn name(self) -> &'static str {
        match self {
            LodMode::Auto => "Auto",
            LodMode::Full => "Full",
            LodMode::Backbone => "Backbone",
            LodMode::CAlpha => "CA",
        }
    }

    pub fn all() -> &'static [LodMode] {
        &[
            LodMode::Auto,
            LodMode::Full,
            LodMode::Backbone,
            LodMode::CAlpha,
        ]
    }

    pub fn next(self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&m| m == self).unwrap_or(0);
        all[(idx + 1) % all.len()]
    }

    pub fn prev(self) -> Self {
        let all = Self::all();
        let idx = all.iter().position(|&m| m == self).unwrap_or(0);
        all[(idx + all.len() - 1) % all.len()]
    }

    pub fn resolve(self, atom_count: usize) -> LodLevel {
        match self {
            LodMode::Full => LodLevel::Full,
            LodMode::Backbone => LodLevel::Backbone,
            LodMode::CAlpha => LodLevel::CAlpha,
            LodMode::Auto => {
                if atom_count >= LOD_CALPHA_ATOMS {
                    LodLevel::CAlpha
                } else if atom_count >= LOD_BACKBONE_ATOMS {
                    LodLevel::Backbone
                } else {
                    LodLevel::Full
                }
            }
        }
    }

    /// Footer/header label, e.g. `Auto/BB` when Auto selected backbone.
    pub fn hud_label(self, atom_count: usize) -> String {
        match self {
            LodMode::Auto => format!("Auto/{}", self.resolve(atom_count).short()),
            other => other.name().to_string(),
        }
    }
}

impl LodLevel {
    pub fn short(self) -> &'static str {
        match self {
            LodLevel::Full => "Full",
            LodLevel::Backbone => "BB",
            LodLevel::CAlpha => "CA",
        }
    }

    /// Spline samples per residue for cartoon ribbons.
    pub fn ribbon_samples(self, residue_count: usize) -> usize {
        let base = match self {
            LodLevel::Full => 6,
            LodLevel::Backbone => 4,
            LodLevel::CAlpha => 3,
        };
        if residue_count > 3000 {
            (base / 2).max(2)
        } else {
            base
        }
    }
}

fn is_non_water_hetatm(atom: &Atom, residue: Option<&Residue>) -> bool {
    if !atom.is_hetatm {
        return false;
    }
    !residue
        .map(|r| r.is_water())
        .unwrap_or_else(|| Residue::name_is_water(&atom.res_name))
}

fn is_nucleic_guide(atom: &Atom) -> bool {
    !atom.is_hetatm && atom.name.trim().eq_ignore_ascii_case("P")
}

/// Whether an already-visible atom should be drawn at this LOD.
pub fn atom_passes_lod(atom: &Atom, residue: Option<&Residue>, level: LodLevel) -> bool {
    match level {
        LodLevel::Full => true,
        LodLevel::Backbone => atom.is_backbone() || is_non_water_hetatm(atom, residue),
        LodLevel::CAlpha => {
            atom.is_c_alpha() || is_nucleic_guide(atom) || is_non_water_hetatm(atom, residue)
        }
    }
}

pub fn atom_drawn(
    atom: &Atom,
    residue: Option<&Residue>,
    visibility: Visibility,
    level: LodLevel,
) -> bool {
    visibility.atom_visible(atom, residue) && atom_passes_lod(atom, residue, level)
}

/// Builds a fast lookup map from atom index to its parent residue.
pub fn build_atom_residue_map(structure: &Structure) -> HashMap<usize, &Residue> {
    let mut map = HashMap::with_capacity(structure.atoms().len());
    for chain in structure.chains() {
        for res in &chain.residues {
            for &atom_idx in &res.atom_indices {
                map.insert(atom_idx, res);
            }
        }
    }
    map
}

/// Precompute the draw color of every atom for the active scheme in a single pass.
///
/// The returned vector is indexed by Atom::index, which is always equal to the
/// atom's position in Structure::atoms (see how Atom::index is assigned in
/// Structure::add_atom and assembly expansion). Building the color table once
/// per frame removes repeated per-atom work -- notably the per-atom chain lookup
/// in ColorScheme::Chain and the per-bond recoloring of both endpoints.
pub fn precompute_atom_colors(structure: &Structure, scheme: ColorScheme) -> Vec<PixelColor> {
    let residue_map = build_atom_residue_map(structure);
    let stats = ColorStats::for_structure(structure);
    let mut out = Vec::with_capacity(structure.atoms().len());
    for atom in structure.atoms() {
        let res = residue_map.get(&atom.index).copied();
        out.push(color_for_atom_with_stats(atom, res, scheme, &stats));
    }
    out
}

/// Like precompute_atom_colors, but reuses an already-built residue map so the
/// per-frame render path builds the map only once instead of twice.
pub fn precompute_atom_colors_with_map(
    structure: &Structure,
    scheme: ColorScheme,
    residue_map: &HashMap<usize, &Residue>,
) -> Vec<PixelColor> {
    let stats = ColorStats::for_structure(structure);
    let mut out = Vec::with_capacity(structure.atoms().len());
    for atom in structure.atoms() {
        let res = residue_map.get(&atom.index).copied();
        out.push(color_for_atom_with_stats(atom, res, scheme, &stats));
    }
    out
}

/// Camera-independent, per-atom view of the scene that the renderers consume.
///
/// colors and visible are indexed by Atom::index (== position in
/// Structure::atoms) and are rebuilt only when the structure, color scheme,
/// visibility, or LOD changes -- NOT when the camera moves. Letting orbit/spin
/// reuse them across frames avoids rebuilding O(n) tables every frame.
pub struct RenderContext<'a> {
    pub structure: &'a Structure,
    pub camera: &'a Camera,
    pub mats: CameraMatrices,
    pub lighting: &'a Lighting,
    pub visibility: Visibility,
    pub lod: LodLevel,
    pub colors: &'a [PixelColor],
    pub visible: &'a [bool],
    pub com: Vec3,
    pub radius: f32,
    /// Largest Van der Waals radius among visible atoms (floored at the unknown
    /// fallback of 1.5 A); used by the VDW whole-scene sub-pixel cull.
    pub max_vdw: f32,
    /// Cached, camera-independent ribbon geometry (None on the one-shot export
    /// path, where render_ribbon builds it fresh).
    pub ribbon_geometry: Option<&'a [ribbon::RibbonPrimitive]>,
}

/// Builds the camera-independent per-atom cache (colors, visibility flags, and
/// the bounding sphere) in a single O(n) pass. Callers should cache the result
/// and rebuild only when the structure, color scheme, visibility, or LOD changes.
pub fn build_render_cache(
    structure: &Structure,
    color_scheme: ColorScheme,
    visibility: Visibility,
    lod: LodMode,
) -> (Vec<PixelColor>, Vec<bool>, Vec3, f32, f32) {
    let level = lod.resolve(structure.atom_count());
    let residue_map = build_atom_residue_map(structure);
    let colors = precompute_atom_colors_with_map(structure, color_scheme, &residue_map);
    let visible: Vec<bool> = structure
        .atoms()
        .iter()
        .map(|a| {
            let res = residue_map.get(&a.index).copied();
            atom_drawn(a, res, visibility, level)
        })
        .collect();
    let com = structure.center_of_mass();
    let radius = structure.bounding_sphere_radius();
    let max_vdw = structure
        .atoms()
        .iter()
        .map(|a| a.vdw_radius())
        .fold(0.0f32, f32::max)
        .max(1.5);
    (colors, visible, com, radius, max_vdw)
}

/// Renders the scene described by ctx into buffer using the given mode.
///
/// Hot path for interactive rendering: callers supply a RenderContext whose
/// colors/visible/com/radius are cached across frames, so only the
/// camera-dependent work (projection + rasterization) runs each frame.
pub fn render_structure_ctx(ctx: &RenderContext, mode: RenderMode, buffer: &mut Framebuffer) {
    if ctx.structure.atoms().is_empty() {
        return;
    }
    match mode {
        RenderMode::Trace => trace::render_trace(ctx, buffer),
        RenderMode::BallAndStick => ball_stick::render_ball_stick(ctx, buffer),
        RenderMode::Ribbon => ribbon::render_ribbon(ctx, buffer),
        RenderMode::Vdw => vdw::render_vdw(ctx, buffer),
        RenderMode::Wireframe => {
            let atoms = ctx.structure.atoms();
            // Projection pass: collect one line primitive per visible bond
            // (in bond order), then hand the list to the band-parallel
            // rasterizer.
            let mut prims: Vec<BandPrimitive> = Vec::new();
            for bond in ctx.structure.bonds() {
                if bond.atom1_idx < atoms.len() && bond.atom2_idx < atoms.len() {
                    let a1 = &atoms[bond.atom1_idx];
                    let a2 = &atoms[bond.atom2_idx];
                    if !ctx.visible[a1.index] || !ctx.visible[a2.index] {
                        continue;
                    }
                    // Segment clip keeps bonds visible when one endpoint is
                    // nearer than the near plane.
                    if let Some((p1, p2)) = ctx.camera.project_segment(
                        &ctx.mats,
                        a1.pos,
                        a2.pos,
                        buffer.width,
                        buffer.height,
                    ) {
                        let c1 = ctx.colors[a1.index];
                        prims.push(BandPrimitive::Line { p1, p2, color: c1 });
                    }
                }
            }
            draw_band_primitives(buffer, &prims, ctx.lighting);
        }
    }
}

/// Orchestrates 3D rendering of a molecular structure into the given framebuffer.
///
/// Builds a fresh RenderContext (colors, visibility, bounding sphere) each call.
/// For interactive use, prefer render_structure_ctx with a cached context so the
/// O(n) per-atom tables are not rebuilt every frame.
#[allow(clippy::too_many_arguments)]
pub fn render_structure(
    structure: &Structure,
    mode: RenderMode,
    color_scheme: ColorScheme,
    camera: &Camera,
    buffer: &mut Framebuffer,
    lighting: &Lighting,
    visibility: Visibility,
    lod: LodMode,
) {
    if structure.atoms().is_empty() {
        return;
    }
    let (colors, visible, com, radius, max_vdw) =
        build_render_cache(structure, color_scheme, visibility, lod);
    let level = lod.resolve(structure.atom_count());
    let mats = camera.matrices();
    let ctx = RenderContext {
        structure,
        camera,
        mats,
        lighting,
        visibility,
        lod: level,
        colors: &colors,
        visible: &visible,
        com,
        radius,
        max_vdw,
        ribbon_geometry: None,
    };
    render_structure_ctx(&ctx, mode, buffer);
}
