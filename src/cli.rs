//! Command-line argument parsing for TermPDB.

use std::path::PathBuf;

use crate::render::{ColorScheme, LodMode, RenderMode};
use clap::Parser;

/// 3D Molecular Structure Viewer in your Terminal
#[derive(Parser, Debug, Clone, PartialEq)]
#[command(
    name = "termpdb",
    author,
    version,
    about = "3D Molecular Structure Viewer in your Terminal",
    long_about = "A high-performance terminal PDB/mmCIF molecular viewer featuring 3D software rasterization, ribbon diagrams, ball-and-stick, trace, and VDW space-filling representations with truecolor ANSI rendering."
)]
pub struct Cli {
    /// Path to .pdb / .cif / .gz file(s) or 4-letter RCSB PDB ID (e.g. 1crn, 1ubq)
    #[arg(value_name = "FILES")]
    pub files: Vec<String>,

    /// Enable high-resolution Kitty Graphics Protocol rendering
    #[arg(long, default_value_t = false)]
    pub kitty: bool,

    /// Resolution scale multiplier for Kitty graphics rendering (e.g. 1.0, 1.5, 2.0)
    #[arg(long = "scale", default_value_t = 1.0)]
    pub scale: f32,

    /// Export rendered frame as Kitty Graphics Protocol escape sequence to file or stdout (-)
    #[arg(long = "export-kitty", value_name = "PATH")]
    pub export_kitty: Option<PathBuf>,

    /// Superimpose input structures using Kabsch RMSD alignment
    #[arg(long = "align", default_value_t = false)]
    pub align: bool,

    /// Force DSSP secondary structure recalculation
    #[arg(long = "dssp", default_value_t = false)]
    pub dssp: bool,

    /// Print planar bond angle in degrees between 3 atoms (e.g. `A:1:CA,A:2:CA,A:3:CA`) and exit
    #[arg(long = "angle", value_name = "SPEC,SPEC,SPEC")]
    pub angle: Option<String>,

    /// Print dihedral angle in degrees between 4 atoms (e.g. `A:1:N,A:1:CA,A:1:C,A:2:N`) and exit
    #[arg(long = "dihedral", value_name = "SPEC,SPEC,SPEC,SPEC")]
    pub dihedral: Option<String>,

    /// Representation rendering mode
    #[arg(short = 'm', long = "mode", value_enum, default_value = "ribbon")]
    pub mode: RenderMode,

    /// Color scheme
    #[arg(short = 'c', long = "color", value_enum, default_value = "rainbow")]
    pub color: ColorScheme,

    /// Start with automatic turntable spin enabled
    #[arg(short = 's', long = "spin", default_value_t = false)]
    pub spin: bool,

    /// Auto-spin rotation speed factor
    #[arg(long = "spin-speed", default_value_t = 1.0)]
    pub spin_speed: f32,

    /// Model serial to display (file numbering). Default: lowest serial in the file.
    #[arg(long = "model", value_name = "N")]
    pub model: Option<i32>,

    /// Biological assembly id (`1`, `2`, …). Default: deposited asymmetric unit. `asu`/`0` for ASU.
    #[arg(long = "assembly", value_name = "ID")]
    pub assembly: Option<String>,

    /// Show water / solvent molecules (hidden by default)
    #[arg(long = "show-waters", default_value_t = false)]
    pub show_waters: bool,

    /// Hide hydrogen atoms
    #[arg(long = "hide-hydrogens", default_value_t = false)]
    pub hide_hydrogens: bool,

    /// Level of detail for large structures (`auto`, `full`, `backbone`, `ca`)
    #[arg(long = "lod", value_enum, default_value = "auto")]
    pub lod: LodMode,

    /// Disable silhouette depth outlines (enabled by default)
    #[arg(long = "no-outline", default_value_t = false)]
    pub no_outline: bool,

    /// Disable screen-space ambient occlusion (SSAO) (enabled by default)
    #[arg(long = "no-ssao", default_value_t = false)]
    pub no_ssao: bool,

    /// Render non-covalent interactions (H-bonds & disulfide bridges)
    #[arg(long = "interactions", default_value_t = false)]
    pub interactions: bool,

    /// Depth-of-Field focal distance
    #[arg(long = "dof", value_name = "DISTANCE")]
    pub dof: Option<f32>,

    /// Print distance in Å between two atoms (`A:12:CA,A:40:N`) and exit unless also exporting
    #[arg(long = "dist", value_name = "SPEC,SPEC")]
    pub dist: Option<String>,

    /// Export rendered frame as ANSI text to file (or "-" for stdout) and exit
    #[arg(long = "export-ansi", value_name = "PATH")]
    pub export_ansi: Option<String>,

    /// Export a supersampled PNG image to PATH and exit (use with --width/--height in pixels)
    #[arg(long = "export-png", value_name = "PATH")]
    pub export_png: Option<String>,

    /// Export a vector SVG image to PATH and exit
    #[arg(long = "export-svg", value_name = "PATH")]
    pub export_svg: Option<String>,

    /// Export a turntable-spin animated GIF to PATH in pure Rust (no ffmpeg needed). Use with --frames/--fps/--width/--height.
    #[arg(long = "export-gif", value_name = "PATH")]
    pub export_gif: Option<String>,

    /// Export a turntable-spin MP4 to PATH (requires ffmpeg). Use --frames/--fps.
    #[arg(long = "export-mp4", value_name = "PATH")]
    pub export_mp4: Option<String>,

    /// Supersampling factor for PNG/MP4 (per axis; 2 = 4x samples/pixel). Higher = smoother edges.
    #[arg(long = "ssaa", default_value_t = 2, value_name = "N")]
    pub ssaa: u8,

    /// Number of frames for MP4 turntable spin (one full rotation)
    #[arg(long = "frames", default_value_t = 60, value_name = "N")]
    pub frames: u32,

    /// Frames per second for the exported MP4
    #[arg(long = "fps", default_value_t = 30, value_name = "FPS")]
    pub fps: u32,

    /// Width (pixels for PNG/MP4, columns for ANSI) for headless export
    #[arg(long = "width", default_value_t = 80, value_name = "W")]
    pub width: u16,

    /// Height (pixels for PNG/MP4, rows for ANSI) for headless export
    #[arg(long = "height", default_value_t = 40, value_name = "H")]
    pub height: u16,
}

impl Cli {
    /// Returns the primary structure source if provided.
    pub fn source(&self) -> Option<&str> {
        self.files.first().map(String::as_str)
    }
}
