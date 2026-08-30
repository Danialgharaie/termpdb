use clap::Parser;
use std::path::PathBuf;
use termpdb::cli::Cli;
use termpdb::model::Structure;
use termpdb::render::{
    ColorScheme, ExportConfig, LodMode, RenderMode, Visibility, export_kitty_frame,
    render_structure_to_framebuffer,
};

#[test]
fn test_export_kitty_frame_generates_valid_sequence() {
    let structure = Structure::default();
    let config = ExportConfig::default();
    let output = export_kitty_frame(&structure, &config, 80, 40).expect("kitty export failed");

    assert!(output.starts_with("\x1b_G"));
    assert!(output.contains("a=T"));
    assert!(output.contains("f=32"));
    assert!(output.ends_with("\x1b\\"));
}

const SAMPLE_PDB_1CRN: &str = r#"HEADER    HYDROLASE                               02-AUG-81   1CRN
TITLE     WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION.
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O
END
"#;

#[test]
fn test_export_kitty_frame_with_atoms() {
    let mut structure = termpdb::parser::parse_pdb(SAMPLE_PDB_1CRN).unwrap();
    structure.ensure_bonds();

    for mode in [
        RenderMode::Ribbon,
        RenderMode::BallAndStick,
        RenderMode::Trace,
        RenderMode::Vdw,
        RenderMode::Wireframe,
    ] {
        let config = ExportConfig {
            mode,
            color: ColorScheme::Rainbow,
            visibility: Visibility::default(),
            lod: LodMode::Auto,
        };
        let output =
            export_kitty_frame(&structure, &config, 40, 20).expect("kitty export with atoms failed");
        assert!(output.starts_with("\x1b_G"));
        assert!(output.contains("a=T"));
        assert!(output.contains("f=32"));
        assert!(output.ends_with("\x1b\\"));
    }
}

#[test]
fn test_render_structure_to_framebuffer() {
    let structure = Structure::default();
    let config = ExportConfig::default();
    let fb = render_structure_to_framebuffer(&structure, &config, 100, 50);
    assert_eq!(fb.width, 100);
    assert_eq!(fb.height, 50);
    let rgba = fb.to_rgba_bytes();
    assert_eq!(rgba.len(), 100 * 50 * 4);
}

#[test]
fn test_cli_kitty_flags_parsing() {
    let cli = Cli::try_parse_from(["termpdb", "1crn", "--kitty"]).unwrap();
    assert!(cli.kitty);
    assert_eq!(cli.export_kitty, None);

    let cli2 = Cli::try_parse_from(["termpdb", "1crn", "--export-kitty", "out.kitty"]).unwrap();
    assert!(!cli2.kitty);
    assert_eq!(cli2.export_kitty, Some(PathBuf::from("out.kitty")));

    let cli3 = Cli::try_parse_from(["termpdb", "1crn", "--export-kitty", "-"]).unwrap();
    assert_eq!(cli3.export_kitty, Some(PathBuf::from("-")));
}
