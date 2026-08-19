//! TermPDB binary entry point.

use std::io::Write;
use clap::Parser;
use termpdb::cli::Cli;
use termpdb::parser::load_structure;
use termpdb::render::export_ansi;
use termpdb::tui;

fn main() {
    let cli = Cli::parse();

    let source = match &cli.source {
        Some(s) if !s.trim().is_empty() => s.trim(),
        _ => {
            eprintln!("Error: No structure source provided.");
            eprintln!("Usage: termpdb <SOURCE> [-m <MODE>] [-c <COLOR>] [-s] [--export-ansi <PATH>]");
            eprintln!("Try `termpdb --help` for more information.");
            std::process::exit(1);
        }
    };

    let structure = match load_structure(source) {
        Ok(s) => s,
        Err(err) => {
            eprintln!("Error loading structure '{}': {}", source, err);
            std::process::exit(1);
        }
    };

    if let Some(export_path) = &cli.export_ansi {
        let ansi = export_ansi(&structure, cli.mode, cli.color, cli.width, cli.height);
        if export_path == "-" {
            let mut stdout = std::io::stdout().lock();
            if let Err(e) = stdout.write_all(ansi.as_bytes()) {
                eprintln!("Error writing ANSI output to stdout: {}", e);
                std::process::exit(1);
            }
        } else if let Err(e) = std::fs::write(export_path, ansi) {
            eprintln!("Error writing ANSI output to '{}': {}", export_path, e);
            std::process::exit(1);
        }
        return;
    }

    if let Err(err) = tui::run(structure, cli.mode, cli.color, cli.spin, cli.spin_speed) {
        eprintln!("Error running termpdb: {}", err);
        std::process::exit(1);
    }
}
