//! TermPDB: Terminal PDB/mmCIF 3D Molecular Structure Viewer.
//!
//! Features 60 FPS 3D software rendering, cartoon ribbon diagrams, ball-and-stick,
//! trace, and VDW space-filling models with truecolor half-block ANSI rendering.

pub mod cli;
pub mod error;
pub mod math;
pub mod model;
pub mod parser;
pub mod render;
pub mod select;
pub mod tui;

// Public re-exports for library consumers
pub use cli::Cli;
pub use error::{Result, TermPdbError};
pub use model::Structure;
pub use parser::load_structure;
pub use render::{ColorScheme, Framebuffer, LodMode, RenderMode, Visibility, export_ansi};
pub use select::{
    AtomSpec, Selection, atom_distance, atom_label, distance_report, parse_atom_spec, resolve_atom,
};
