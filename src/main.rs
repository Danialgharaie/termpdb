//! TermPDB binary entry point.

use clap::Parser;
use std::io::Write;
use termpdb::cli::Cli;
use termpdb::parser::load_structure;
use termpdb::render::{Visibility, export_ansi_with_visibility};
use termpdb::select::distance_report;
use termpdb::tui;

fn main() {
    let cli = Cli::parse();

    let source = match cli.source() {
        Some(s) if !s.trim().is_empty() => s.trim(),
        _ => {
            eprintln!("Error: No structure source provided.");
            eprintln!(
                "Usage: termpdb <FILES> [-m <MODE>] [-c <COLOR>] [-s] [--align] [--export-ansi <PATH>]"
            );
            eprintln!("Try `termpdb --help` for more information.");
            std::process::exit(1);
        }
    };

    let mut structure = match load_structure(source) {
        Ok(s) => s,
        Err(err) => {
            eprintln!("Error loading structure '{}': {}", source, err);
            std::process::exit(1);
        }
    };

    // If multiple files provided with --align, load target and superimpose
    if cli.files.len() > 1 && cli.align {
        for target_file in &cli.files[1..] {
            match load_structure(target_file) {
                Ok(target_struct) => {
                    if let Some(res) = termpdb::model::align::superimpose_structures(
                        &mut structure,
                        &target_struct,
                    ) {
                        println!(
                            "Superposition ({} onto {}): {} aligned pairs, RMSD = {:.3} Å",
                            source, target_file, res.aligned_pairs, res.kabsch.rmsd
                        );
                    } else {
                        eprintln!(
                            "Warning: could not superimpose {} onto {}",
                            source, target_file
                        );
                    }
                }
                Err(err) => eprintln!("Error loading target structure '{}': {}", target_file, err),
            }
        }
    }

    if cli.dssp {
        let count = termpdb::model::dssp::assign_dssp(&mut structure);
        println!("DSSP: assigned secondary structure to {} residues", count);
    }

    if let Some(serial) = cli.model {
        structure.set_active_model(serial).unwrap_or_else(|err| {
            eprintln!("Error: {}", err);
            std::process::exit(1);
        });
    }

    if let Some(assembly) = &cli.assembly {
        structure
            .set_assembly(Some(assembly.as_str()))
            .unwrap_or_else(|err| {
                eprintln!("Error: {}", err);
                std::process::exit(1);
            });
    }

    // Detect covalent bonds once (if the file lacks CONECT records) so the
    // ball-and-stick / trace / ribbon renderers never re-run detection per frame.
    structure.ensure_bonds();

    let visibility = Visibility {
        show_waters: cli.show_waters,
        show_hydrogens: !cli.hide_hydrogens,
    };

    if let Some(dist) = &cli.dist {
        match distance_report(&structure, dist) {
            Ok(line) => println!("{line}"),
            Err(err) => {
                eprintln!("Error: {err}");
                std::process::exit(1);
            }
        }
        if cli.export_ansi.is_none()
            && cli.export_kitty.is_none()
            && cli.export_png.is_none()
            && cli.export_svg.is_none()
            && cli.export_mp4.is_none()
        {
            return;
        }
    }

    if let Some(angle) = &cli.angle {
        match termpdb::select::angle_report(&structure, angle) {
            Ok(line) => println!("{line}"),
            Err(err) => {
                eprintln!("Error: {err}");
                std::process::exit(1);
            }
        }
        if cli.export_ansi.is_none()
            && cli.export_kitty.is_none()
            && cli.export_png.is_none()
            && cli.export_svg.is_none()
            && cli.export_mp4.is_none()
        {
            return;
        }
    }

    if let Some(dihedral) = &cli.dihedral {
        match termpdb::select::dihedral_report(&structure, dihedral) {
            Ok(line) => println!("{line}"),
            Err(err) => {
                eprintln!("Error: {err}");
                std::process::exit(1);
            }
        }
        if cli.export_ansi.is_none()
            && cli.export_kitty.is_none()
            && cli.export_png.is_none()
            && cli.export_svg.is_none()
            && cli.export_mp4.is_none()
        {
            return;
        }
    }

    if let Some(export_path) = &cli.export_ansi {
        let ansi = export_ansi_with_visibility(
            &structure, cli.mode, cli.color, cli.width, cli.height, visibility, cli.lod,
        );
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

    if let Some(path) = &cli.export_kitty {
        let export_config = termpdb::render::ExportConfig {
            mode: cli.mode,
            color: cli.color,
            visibility,
            lod: cli.lod,
        };
        let (cols, rows) = crossterm::terminal::size().unwrap_or((cli.width, cli.height));
        let (cols, rows) = if cli.width != 80 || cli.height != 40 {
            (cli.width, cli.height)
        } else {
            (cols, rows)
        };
        match termpdb::render::export_kitty_frame(&structure, &export_config, cols, rows) {
            Ok(kitty_str) => {
                if path.as_os_str() == "-" {
                    let mut stdout = std::io::stdout().lock();
                    if let Err(e) = stdout.write_all(kitty_str.as_bytes()) {
                        eprintln!("Error writing Kitty output to stdout: {}", e);
                        std::process::exit(1);
                    }
                } else if let Err(e) = std::fs::write(path, kitty_str) {
                    eprintln!("Error exporting Kitty image to '{}': {}", path.display(), e);
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("Error exporting Kitty image: {}", e);
                std::process::exit(1);
            }
        }
        return;
    }

    if let Some(png_path) = &cli.export_png {
        let rgba = termpdb::render::render_supersampled(
            &structure,
            cli.mode,
            cli.color,
            cli.width as usize,
            cli.height as usize,
            cli.ssaa as usize,
            visibility,
            cli.lod,
        );
        if let Err(e) =
            termpdb::render::write_png(png_path, &rgba, cli.width as u32, cli.height as u32)
        {
            eprintln!("Error exporting PNG image to '{}': {}", png_path, e);
            std::process::exit(1);
        }
        return;
    }

    if let Some(svg_path) = &cli.export_svg {
        let svg = termpdb::render::render_svg(
            &structure,
            cli.mode,
            cli.color,
            cli.width as usize,
            cli.height as usize,
            visibility,
            cli.lod,
        );
        if let Err(e) = std::fs::write(svg_path, svg) {
            eprintln!("Error exporting SVG image to '{}': {}", svg_path, e);
            std::process::exit(1);
        }
        return;
    }

    if let Some(mp4_path) = &cli.export_mp4 {
        if let Err(e) = termpdb::render::export_mp4(
            &structure,
            cli.mode,
            cli.color,
            cli.width as usize,
            cli.height as usize,
            cli.ssaa as usize,
            cli.frames,
            cli.fps,
            visibility,
            cli.lod,
            mp4_path,
        ) {
            eprintln!("Error exporting MP4 video to '{}': {}", mp4_path, e);
            std::process::exit(1);
        }
        return;
    }

    let postprocess = termpdb::render::PostProcessConfig {
        outline: !cli.no_outline,
        ssao: !cli.no_ssao,
        outline_threshold: 0.12,
        ssao_radius: 2,
    };

    let graphics_backend = if cli.kitty {
        termpdb::render::GraphicsBackend::Kitty
    } else {
        termpdb::render::GraphicsBackend::HalfBlock
    };

    if let Err(err) = tui::run(
        structure,
        cli.mode,
        cli.color,
        cli.spin,
        cli.spin_speed,
        visibility,
        cli.lod,
        postprocess,
        cli.interactions,
        cli.dof,
        graphics_backend,
    ) {
        eprintln!("Error running termpdb: {}", err);
        std::process::exit(1);
    }
}
