use clap::Parser;
use termpdb::cli::Cli;
use termpdb::parser::{parse_cif, parse_pdb};
use termpdb::render::{ColorScheme, RenderMode, export_ansi};

const SAMPLE_PDB_1CRN: &str = r#"HEADER    HYDROLASE                               02-AUG-81   1CRN
TITLE     WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION.
HELIX    1  H1 THR A    7  GLY A   17  1                                  11
SHEET    1  S1 2 CYS A   1  CYS A    4  0
SHEET    2  S1 2 CYS A  32  ILE A   35 -1
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O
ATOM      5  CB  THR A   1      18.170  12.703   5.337  1.00 13.02           C
ATOM      6  OG1 THR A   1      19.334  12.829   4.463  1.00 15.06           O
ATOM      7  CG2 THR A   1      18.150  11.454   6.224  1.00 13.06           C
ATOM      8  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N
ATOM      9  CA  THR A   2      13.856  11.469   6.066  1.00  8.28           C
ATOM     10  CA  GLY A   7      10.000  10.000  10.000  1.00 15.00           C
ATOM     11  CA  ILE A  33      20.000  20.000  20.000  1.00 12.00           C
HETATM   12  O   HOH A  47      15.021  -1.423   2.492  1.00 25.00           O
CONECT   12    1
END
"#;

const SAMPLE_CIF_1CRN: &str = r#"data_1CRN
_entry.id   1CRN
_struct.title 'WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION.'
loop_
_struct_conf.conf_type_id
_struct_conf.id
_struct_conf.beg_auth_asym_id
_struct_conf.beg_auth_seq_id
_struct_conf.end_auth_asym_id
_struct_conf.end_auth_seq_id
HELX_P HELX_P1 A 7 A 17
loop_
_struct_sheet_range.sheet_id
_struct_sheet_range.id
_struct_sheet_range.beg_auth_asym_id
_struct_sheet_range.beg_auth_seq_id
_struct_sheet_range.end_auth_asym_id
_struct_sheet_range.end_auth_seq_id
S1 1 A 1 A 4
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . THR A 1 17.047 14.099 3.625 1.00 13.79 1 THR A N 1
ATOM 2 C CA . THR A 1 16.967 12.784 4.338 1.00 10.80 1 THR A CA 1
ATOM 3 C C . THR A 1 15.685 12.755 5.133 1.00 9.19 1 THR A C 1
ATOM 4 O O . THR A 1 15.268 13.825 5.594 1.00 9.85 1 THR A O 1
ATOM 5 C CB . THR A 1 18.170 12.703 5.337 1.00 13.02 1 THR A CB 1
ATOM 6 O OG1 . THR A 1 19.334 12.829 4.463 1.00 15.06 1 THR A OG1 1
ATOM 7 C CG2 . THR A 1 18.150 11.454 6.224 1.00 13.06 1 THR A CG2 1
ATOM 8 N N . THR A 2 15.115 11.555 5.265 1.00 7.81 2 THR A N 1
ATOM 9 C CA . THR A 2 13.856 11.469 6.066 1.00 8.28 2 THR A CA 1
ATOM 10 C CA . GLY A 7 10.000 10.000 10.000 1.00 15.00 7 GLY A CA 1
ATOM 11 C CA . ILE A 33 20.000 20.000 20.000 1.00 12.00 33 ILE A CA 1
HETATM 12 O O . HOH A . 15.021 -1.423 2.492 1.00 25.00 47 HOH A O 1
#
"#;

#[test]
fn test_pdb_parsing_to_ansi_export_all_modes() {
    let mut structure = parse_pdb(SAMPLE_PDB_1CRN).expect("Failed to parse PDB");
    structure.build_bonds();

    for &mode in RenderMode::all() {
        let ansi = export_ansi(&structure, mode, ColorScheme::Rainbow, 80, 40);
        assert!(
            !ansi.is_empty(),
            "Exported ANSI for mode {:?} should not be empty",
            mode
        );
        assert!(
            ansi.contains('▀'),
            "ANSI output should contain upper half-block character ▀"
        );
        assert!(
            ansi.contains("\x1b[38;2;"),
            "ANSI output should contain truecolor foreground escapes"
        );
        assert!(
            ansi.contains("\x1b[0m"),
            "ANSI output should contain reset escapes"
        );

        let lines: Vec<&str> = ansi.lines().collect();
        assert_eq!(
            lines.len(),
            40,
            "ANSI output line count should match requested height 40"
        );
    }
}

#[test]
fn test_cif_parsing_to_ribbon_ansi_export() {
    let mut structure = parse_cif(SAMPLE_CIF_1CRN).expect("Failed to parse CIF");
    structure.build_bonds();

    let ansi = export_ansi(
        &structure,
        RenderMode::Ribbon,
        ColorScheme::SecondaryStructure,
        60,
        30,
    );
    assert!(!ansi.is_empty());
    assert!(ansi.contains('▀'));
    assert!(ansi.contains("\x1b[38;2;"));
    assert!(ansi.contains("\x1b[0m"));

    let lines: Vec<&str> = ansi.lines().collect();
    assert_eq!(
        lines.len(),
        30,
        "ANSI output line count should match requested height 30"
    );
}

#[test]
fn test_all_color_schemes_ansi_export() {
    let mut structure = parse_pdb(SAMPLE_PDB_1CRN).expect("Failed to parse PDB");
    structure.build_bonds();

    for &scheme in ColorScheme::all() {
        let ansi = export_ansi(&structure, RenderMode::Ribbon, scheme, 80, 40);
        assert!(
            !ansi.is_empty(),
            "Exported ANSI for scheme {:?} should not be empty",
            scheme
        );
        assert!(ansi.contains('▀'));
        assert!(ansi.contains("\x1b[38;2;"));
    }
}

#[test]
fn test_cli_argument_parsing_defaults() {
    let cli = Cli::try_parse_from(["termpdb", "1crn"]).expect("CLI parsing failed");
    assert_eq!(cli.source.as_deref(), Some("1crn"));
    assert_eq!(cli.mode, RenderMode::Ribbon);
    assert_eq!(cli.color, ColorScheme::Rainbow);
    assert!(!cli.spin);
    assert_eq!(cli.spin_speed, 1.0);
    assert_eq!(cli.export_ansi, None);
    assert_eq!(cli.width, 80);
    assert_eq!(cli.height, 40);
}

#[test]
fn test_cli_argument_parsing_flags_and_aliases() {
    let cli = Cli::try_parse_from([
        "termpdb",
        "sample.pdb",
        "-m",
        "trace",
        "-c",
        "ss",
        "-s",
        "--spin-speed",
        "2.5",
        "--export-ansi",
        "-",
        "--width",
        "120",
        "--height",
        "60",
    ])
    .expect("CLI parsing with flags failed");

    assert_eq!(cli.source.as_deref(), Some("sample.pdb"));
    assert_eq!(cli.mode, RenderMode::Trace);
    assert_eq!(cli.color, ColorScheme::SecondaryStructure);
    assert!(cli.spin);
    assert_eq!(cli.spin_speed, 2.5);
    assert_eq!(cli.export_ansi.as_deref(), Some("-"));
    assert_eq!(cli.width, 120);
    assert_eq!(cli.height, 60);

    // Test BallAndStick alias
    let cli_bas = Cli::try_parse_from(["termpdb", "1crn", "-m", "ball-and-stick"]).unwrap();
    assert_eq!(cli_bas.mode, RenderMode::BallAndStick);

    let cli_bas2 = Cli::try_parse_from(["termpdb", "1crn", "-m", "ball_and_stick"]).unwrap();
    assert_eq!(cli_bas2.mode, RenderMode::BallAndStick);

    // Test VDW mode
    let cli_vdw = Cli::try_parse_from(["termpdb", "1crn", "-m", "vdw"]).unwrap();
    assert_eq!(cli_vdw.mode, RenderMode::Vdw);

    // Test SecondaryStructure alias "secondary-structure"
    let cli_ss = Cli::try_parse_from(["termpdb", "1crn", "-c", "secondary-structure"]).unwrap();
    assert_eq!(cli_ss.color, ColorScheme::SecondaryStructure);

    // Test BFactor and b-factor alias
    let cli_bf = Cli::try_parse_from(["termpdb", "1crn", "-c", "bfactor"]).unwrap();
    assert_eq!(cli_bf.color, ColorScheme::BFactor);

    let cli_bf2 = Cli::try_parse_from(["termpdb", "1crn", "-c", "b-factor"]).unwrap();
    assert_eq!(cli_bf2.color, ColorScheme::BFactor);

    // Test Hydrophobicity
    let cli_hyd = Cli::try_parse_from(["termpdb", "1crn", "-c", "hydrophobicity"]).unwrap();
    assert_eq!(cli_hyd.color, ColorScheme::Hydrophobicity);
}

#[test]
fn test_cli_no_arguments() {
    let cli = Cli::try_parse_from(["termpdb"])
        .expect("CLI parsing without arguments should succeed with None source");
    assert_eq!(cli.source, None);
}

#[test]
fn test_export_ansi_zero_dimensions() {
    let structure = parse_pdb(SAMPLE_PDB_1CRN).expect("Failed to parse PDB");
    assert_eq!(
        export_ansi(&structure, RenderMode::Ribbon, ColorScheme::Rainbow, 0, 40),
        ""
    );
    assert_eq!(
        export_ansi(&structure, RenderMode::Ribbon, ColorScheme::Rainbow, 80, 0),
        ""
    );
}

#[test]
fn test_export_ansi_file_roundtrip() {
    let structure = parse_pdb(SAMPLE_PDB_1CRN).expect("Failed to parse PDB");
    let ansi = export_ansi(&structure, RenderMode::Ribbon, ColorScheme::Rainbow, 40, 20);

    let pid = std::process::id();
    let timestamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    let tmp_path =
        std::env::temp_dir().join(format!("termpdb_test_export_{}_{}.ansi", pid, timestamp));
    std::fs::write(&tmp_path, &ansi).expect("Failed to write temporary file");
    let read_back = std::fs::read_to_string(&tmp_path).expect("Failed to read temporary file");
    let _ = std::fs::remove_file(&tmp_path);

    assert_eq!(ansi, read_back);
    assert_eq!(read_back.lines().count(), 20);
}
