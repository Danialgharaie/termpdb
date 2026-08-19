//! Command-line argument parsing for TermPDB.

use crate::render::{ColorScheme, RenderMode};
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
    /// Path to .pdb / .cif / .gz file or 4-letter RCSB PDB ID (e.g. 1crn, 1ubq)
    #[arg(value_name = "SOURCE")]
    pub source: Option<String>,

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

    /// Export rendered frame as ANSI text to file (or "-" for stdout) and exit
    #[arg(long = "export-ansi", value_name = "PATH")]
    pub export_ansi: Option<String>,

    /// Terminal width (columns) for headless export
    #[arg(long = "width", default_value_t = 80, value_name = "COLS")]
    pub width: u16,

    /// Terminal height (rows) for headless export
    #[arg(long = "height", default_value_t = 40, value_name = "ROWS")]
    pub height: u16,
}
