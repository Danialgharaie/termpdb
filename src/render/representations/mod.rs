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
pub fn project_radius(world_radius: f32, view_depth: f32, fov: f32, height: usize) -> f32 {
    let tan_half = (fov * 0.5).tan();
    if view_depth <= 1e-4 || tan_half <= 1e-4 {
        0.0
    } else {
        (world_radius / (view_depth * tan_half)) * (height as f32 * 0.5)
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
            for bond in ctx.structure.bonds() {
                if bond.atom1_idx < atoms.len() && bond.atom2_idx < atoms.len() {
                    let a1 = &atoms[bond.atom1_idx];
                    let a2 = &atoms[bond.atom2_idx];
                    if !ctx.visible[a1.index] || !ctx.visible[a2.index] {
                        continue;
                    }
                    if let (Some(p1), Some(p2)) = (
                        ctx.camera
                            .project(&ctx.mats, a1.pos, buffer.width, buffer.height),
                        ctx.camera
                            .project(&ctx.mats, a2.pos, buffer.width, buffer.height),
                    ) {
                        let c1 = ctx.colors[a1.index];
                        crate::render::rasterizer::draw_line_3d(buffer, p1, p2, c1);
                    }
                }
            }
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
