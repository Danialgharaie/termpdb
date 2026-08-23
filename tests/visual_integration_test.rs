use clap::Parser;
use termpdb::cli::Cli;
use termpdb::render::{
    export_ansi_with_visibility, ColorScheme, LodMode, RenderMode, Visibility,
};

const SAMPLE_1CRN: &str = r#"HEADER    PLANT SEED PROTEIN                     30-APR-81   1CRN
COMPND    MOL_ID: 1; MOLECULE: CRAMBIN; CHAIN: A; ENGINEERED: NO
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O
ATOM      5  CB  THR A   1      18.170  12.703   5.337  1.00 13.02           C
ATOM      6  OG1 THR A   1      19.334  12.829   4.463  1.00 15.06           O
ATOM      7  CG2 THR A   1      18.150  11.546   6.304  1.00 19.25           C
ATOM      8  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N
ATOM      9  CA  THR A   2      13.856  11.469   6.066  1.00  7.14           C
ATOM     10  C   THR A   2      14.164  10.785   7.379  1.00  5.84           C
ATOM     11  O   THR A   2      14.993   9.862   7.448  1.00  6.94           O
ATOM     12  CB  THR A   2      12.732  10.711   5.261  1.00  8.85           C
ATOM     13  OG1 THR A   2      12.443  11.442   4.070  1.00 10.70           O
ATOM     14  CG2 THR A   2      11.450  10.589   6.072  1.00 11.46           C
ATOM     15  N   CYS A   3      13.488  11.241   8.417  1.00  5.24           N
ATOM     16  CA  CYS A   3      13.660  10.666   9.727  1.00  4.70           C
ATOM     17  C   CYS A   3      12.600   9.619   9.992  1.00  4.56           C
ATOM     18  O   CYS A   3      11.531   9.655   9.386  1.00  5.21           O
ATOM     19  CB  CYS A   3      13.635  11.758  10.793  1.00  5.37           C
ATOM     20  SG  CYS A   3      15.006  12.923  10.709  1.00  7.07           S
ATOM     21  N   CYS A   4      12.915   8.687  10.887  1.00  4.47           N
ATOM     22  CA  CYS A   4      12.015   7.580  11.196  1.00  4.34           C
ATOM     23  C   CYS A   4      10.826   7.585  10.239  1.00  4.11           C
ATOM     24  O   CYS A   4       9.805   6.963  10.518  1.00  4.61           O
ATOM     25  CB  CYS A   4      12.793   6.262  11.173  1.00  5.04           C
ATOM     26  SG  CYS A   4      14.150   6.257  12.378  1.00  6.25           S
TER      27      CYS A   4
END
"#;

#[test]
fn test_all_12_color_schemes_cycle_and_render() {
    let structure = termpdb::parser::parse_pdb(SAMPLE_1CRN).unwrap();
    let schemes = ColorScheme::all();
    assert_eq!(schemes.len(), 12);

    for &scheme in schemes {
        let ansi = export_ansi_with_visibility(
            &structure,
            RenderMode::Ribbon,
            scheme,
            60,
            30,
            Visibility::default(),
            LodMode::Auto,
        );
        assert!(!ansi.is_empty(), "Color scheme {:?} should export ANSI", scheme);
        assert!(ansi.contains("\x1b[38;2;"));
    }
}

#[test]
fn test_wireframe_mode_ansi_export() {
    let structure = termpdb::parser::parse_pdb(SAMPLE_1CRN).unwrap();
    let ansi = export_ansi_with_visibility(
        &structure,
        RenderMode::Wireframe,
        ColorScheme::Rainbow,
        60,
        30,
        Visibility::default(),
        LodMode::Auto,
    );
    assert!(!ansi.is_empty());
}

#[test]
fn test_cli_parsing_new_visual_flags() {
    let args = vec![
        "termpdb",
        "1crn.pdb",
        "--mode",
        "wireframe",
        "--color",
        "plddt",
        "--no-outline",
        "--no-ssao",
        "--interactions",
        "--dof",
        "15.5",
    ];

    let cli = Cli::try_parse_from(args).unwrap();
    assert_eq!(cli.mode, RenderMode::Wireframe);
    assert_eq!(cli.color, ColorScheme::Plddt);
    assert!(cli.no_outline);
    assert!(cli.no_ssao);
    assert!(cli.interactions);
    assert_eq!(cli.dof, Some(15.5));
}
