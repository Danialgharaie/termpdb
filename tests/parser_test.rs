use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::Write;
use termpdb::model::SecondaryStructure;
use termpdb::parser::rcsb::fetch_pdb;
use termpdb::parser::{load_structure, parse_cif, parse_pdb};

const SAMPLE_PDB_1CRN: &str = r#"HEADER    HYDROLASE                               02-AUG-81   1CRN
TITLE     WATER STRUCTURE OF A HYDROPHOBIC PROTEIN AT ATOMIC RESOLUTION.
TITLE    2 PENTANE-FORMAMIDE CRYSTALS
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

const SAMPLE_MULTI_CHAIN_PDB: &str = r#"HEADER    TEST MULTI-CHAIN                        01-JAN-24   9XYZ
TITLE     MULTI CHAIN COMPLEX
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 10.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 10.00           C
ATOM      4  O   ALA A   1       1.246   2.389   0.000  1.00 10.00           O
ATOM      5  N   GLY B   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      6  CA  GLY B   1      11.458  10.000  10.000  1.00 20.00           C
ATOM      7  C   GLY B   1      12.009  11.420  10.000  1.00 20.00           C
ATOM      8  O   GLY B   1      11.246  12.389  10.000  1.00 20.00           O
END
"#;

#[test]
fn test_parse_pdb_basic_1crn() {
    let structure = parse_pdb(SAMPLE_PDB_1CRN).expect("Failed to parse 1CRN PDB");

    assert_eq!(structure.id_code.as_deref(), Some("1CRN"));
    assert!(structure.title.contains("WATER STRUCTURE OF A HYDROPHOBIC"));
    assert_eq!(structure.atom_count(), 12);
    assert_eq!(structure.chain_count(), 1);

    let chain_a = structure.get_chain("A").expect("Chain A missing");
    assert!(chain_a.residue_count() >= 4);

    // Verify first atom
    let a1 = &structure.atoms()[0];
    assert_eq!(a1.name, "N");
    assert_eq!(a1.res_name, "THR");
    assert_eq!(a1.res_seq, 1);
    assert_eq!(a1.chain_id, "A");
    assert!(!a1.is_hetatm);
    assert_eq!(a1.element.symbol, "N");
    assert!((a1.pos.x - 17.047).abs() < 1e-4);
    assert!((a1.pos.y - 14.099).abs() < 1e-4);
    assert!((a1.pos.z - 3.625).abs() < 1e-4);
    assert!((a1.occupancy - 1.0).abs() < 1e-4);
    assert!((a1.b_factor - 13.79).abs() < 1e-4);

    // Verify HETATM
    let a12 = &structure.atoms()[11];
    assert_eq!(a12.name, "O");
    assert_eq!(a12.res_name, "HOH");
    assert_eq!(a12.res_seq, 47);
    assert!(a12.is_hetatm);
    assert_eq!(a12.element.symbol, "O");

    // Verify Secondary Structure
    // Residue 1 should be Sheet (1..4)
    let res1 = chain_a.get_residue(1).expect("Residue 1 missing");
    assert_eq!(res1.secondary_structure, SecondaryStructure::Sheet);

    // Residue 7 should be Helix (7..17)
    let res7 = chain_a.get_residue(7).expect("Residue 7 missing");
    assert_eq!(res7.secondary_structure, SecondaryStructure::Helix);

    // Verify CA atoms
    let ca_atoms = structure.ca_atoms();
    assert_eq!(ca_atoms.len(), 4); // residues 1, 2, 7, 33
}

#[test]
fn test_parse_pdb_multi_chain() {
    let structure = parse_pdb(SAMPLE_MULTI_CHAIN_PDB).expect("Failed to parse multi-chain PDB");
    assert_eq!(structure.chain_count(), 2);
    assert!(structure.get_chain("A").is_some());
    assert!(structure.get_chain("B").is_some());
    assert_eq!(structure.atom_count(), 8);

    let chain_a = structure.get_chain("A").unwrap();
    let chain_b = structure.get_chain("B").unwrap();
    assert_eq!(chain_a.residues.len(), 1);
    assert_eq!(chain_b.residues.len(), 1);
    assert_eq!(chain_a.residues[0].name, "ALA");
    assert_eq!(chain_b.residues[0].name, "GLY");
}

#[test]
fn test_parse_pdb_element_deduction_when_empty() {
    let pdb_no_elem = r#"ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00
ATOM      2  FE  HEM A   2       1.000   1.000   1.000  1.00 10.00
ATOM      3 1HD1 ILE A   3       2.000   2.000   2.000  1.00 10.00
ATOM      4  ZN   ZN A   4       3.000   3.000   3.000  1.00 10.00
END
"#;
    let structure = parse_pdb(pdb_no_elem).expect("Failed to parse PDB without element columns");
    assert_eq!(structure.atoms()[0].element.symbol, "C");
    assert_eq!(structure.atoms()[1].element.symbol, "Fe");
    assert_eq!(structure.atoms()[2].element.symbol, "H");
    assert_eq!(structure.atoms()[3].element.symbol, "Zn");
}

#[test]
fn test_parse_pdb_alt_loc_and_charge_and_conect() {
    let pdb_text = r#"ATOM      1  N  AALA A   1       0.000   0.000   0.000  0.50 10.00           N1+
ATOM      2  N  BALA A   1       0.100   0.100   0.100  0.50 10.00           N1+
HETATM    3  ZN   ZN A   2      10.000  10.000  10.000  1.00 20.00          ZN2+
CONECT    1    3
END
"#;
    let structure = parse_pdb(pdb_text).expect("Failed to parse PDB with alt loc and charge");
    assert_eq!(structure.atoms()[0].alt_loc, Some('A'));
    assert_eq!(structure.atoms()[0].charge, Some(1));
    assert_eq!(structure.atoms()[1].alt_loc, Some('B'));
    assert_eq!(structure.atoms()[2].charge, Some(2));
    assert_eq!(structure.atoms()[2].element.symbol, "Zn");

    // Verify CONECT bond was added between atom 0 (serial 1) and atom 2 (serial 3)
    assert!(
        structure
            .bonds()
            .iter()
            .any(|b| (b.atom1_idx == 0 && b.atom2_idx == 2)
                || (b.atom1_idx == 2 && b.atom2_idx == 0))
    );
}

#[test]
fn test_parse_cif_basic_1crn() {
    let structure = parse_cif(SAMPLE_CIF_1CRN).expect("Failed to parse 1CRN CIF");

    assert_eq!(structure.id_code.as_deref(), Some("1CRN"));
    assert!(structure.title.contains("WATER STRUCTURE OF A HYDROPHOBIC"));
    assert_eq!(structure.atom_count(), 12);
    assert_eq!(structure.chain_count(), 1);

    let chain_a = structure.get_chain("A").expect("Chain A missing");
    let res1 = chain_a.get_residue(1).expect("Residue 1 missing");
    assert_eq!(res1.secondary_structure, SecondaryStructure::Sheet);

    let res7 = chain_a.get_residue(7).expect("Residue 7 missing");
    assert_eq!(res7.secondary_structure, SecondaryStructure::Helix);

    let a1 = &structure.atoms()[0];
    assert_eq!(a1.name, "N");
    assert_eq!(a1.res_name, "THR");
    assert_eq!(a1.res_seq, 1);
    assert_eq!(a1.chain_id, "A");
    assert_eq!(a1.element.symbol, "N");
    assert!((a1.pos.x - 17.047).abs() < 1e-4);
    assert!((a1.pos.y - 14.099).abs() < 1e-4);
    assert!((a1.pos.z - 3.625).abs() < 1e-4);
}

#[test]
fn test_parse_cif_multiline_and_quotes() {
    let cif_text = r#"data_TEST
_entry.id TEST
_struct.title
;A multi-line
title structure
for testing.
;
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C "C A" ALA A 1 0.0 0.0 0.0
ATOM 2 O 'O B' ALA A 1 1.0 1.0 1.0
#
"#;
    let structure = parse_cif(cif_text).expect("Failed to parse CIF with quotes");
    assert_eq!(structure.id_code.as_deref(), Some("TEST"));
    assert!(structure.title.contains("multi-line"));
    assert_eq!(structure.atom_count(), 2);
    assert_eq!(structure.atoms()[0].name, "C A");
    assert_eq!(structure.atoms()[1].name, "O B");
}

#[test]
fn test_load_structure_from_file_and_gz() {
    let tmp_dir = std::env::temp_dir();
    let pdb_file = tmp_dir.join("test_termpdb_sample.pdb");
    let gz_file = tmp_dir.join("test_termpdb_sample.pdb.gz");
    let cif_file = tmp_dir.join("test_termpdb_sample.cif");

    std::fs::write(&pdb_file, SAMPLE_PDB_1CRN).unwrap();
    std::fs::write(&cif_file, SAMPLE_CIF_1CRN).unwrap();

    // Compress PDB to .gz
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(SAMPLE_PDB_1CRN.as_bytes()).unwrap();
    let gz_bytes = encoder.finish().unwrap();
    std::fs::write(&gz_file, gz_bytes).unwrap();

    // Test loading regular PDB file
    let s_pdb = load_structure(pdb_file.to_str().unwrap()).expect("Failed to load PDB file");
    assert_eq!(s_pdb.atom_count(), 12);

    // Test loading CIF file
    let s_cif = load_structure(cif_file.to_str().unwrap()).expect("Failed to load CIF file");
    assert_eq!(s_cif.atom_count(), 12);

    // Test loading .gz file
    let s_gz = load_structure(gz_file.to_str().unwrap()).expect("Failed to load .pdb.gz file");
    assert_eq!(s_gz.atom_count(), 12);

    // Cleanup
    let _ = std::fs::remove_file(pdb_file);
    let _ = std::fs::remove_file(gz_file);
    let _ = std::fs::remove_file(cif_file);
}

#[test]
fn test_load_structure_invalid_source() {
    let result = load_structure("non_existent_file_path_123456789.pdb");
    assert!(result.is_err());
}

#[test]
fn test_fetch_pdb_invalid_id() {
    let result = fetch_pdb("!INVALID@#");
    assert!(result.is_err());
}

const SAMPLE_MULTI_MODEL_PDB: &str = r#"HEADER    ENSEMBLE                                01-JAN-24   1ENS
TITLE     THREE MODELS WITH GAPPED SERIALS
MODEL        1
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      2  C   ALA A   1       1.500   0.000   0.000  1.00 10.00           C
ENDMDL
MODEL        2
ATOM      1  CA  ALA A   1       5.000   0.000   0.000  1.00 10.00           C
ATOM      2  C   ALA A   1       6.500   0.000   0.000  1.00 10.00           C
ENDMDL
MODEL        5
ATOM      1  CA  ALA A   1      10.000   0.000   0.000  1.00 10.00           C
ENDMDL
END
"#;

const SAMPLE_MULTI_MODEL_CIF: &str = r#"data_1ENS
_entry.id 1ENS
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA ALA A 1 0.0 0.0 0.0 1
ATOM 2 C C ALA A 1 1.5 0.0 0.0 1
ATOM 1 C CA ALA A 1 5.0 0.0 0.0 2
ATOM 2 C C ALA A 1 6.5 0.0 0.0 2
ATOM 1 C CA ALA A 1 10.0 0.0 0.0 5
#
"#;

#[test]
fn test_parse_pdb_multi_model_defaults_to_lowest_serial() {
    let mut structure = parse_pdb(SAMPLE_MULTI_MODEL_PDB).expect("parse multi-model PDB");

    assert_eq!(structure.model_serials(), vec![1, 2, 5]);
    assert_eq!(structure.active_model_serial(), 1);
    assert_eq!(structure.atom_count(), 2);
    assert!((structure.atoms()[0].pos.x - 0.0).abs() < 1e-4);
    assert_eq!(structure.bonds().len(), 1);

    structure.set_active_model(5).expect("model 5 exists");
    assert_eq!(structure.atom_count(), 1);
    assert!((structure.atoms()[0].pos.x - 10.0).abs() < 1e-4);
    assert!(structure.bonds().is_empty());

    let err = structure.set_active_model(3).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("3"), "{msg}");
    assert!(msg.contains("1"), "{msg}");
}

#[test]
fn test_parse_pdb_multi_model_bonds_do_not_cross_models() {
    let structure = parse_pdb(SAMPLE_MULTI_MODEL_PDB).unwrap();
    // Each model has its own atom indices starting at 0; model 1 has a 1.5 Å C-C bond.
    assert_eq!(structure.bonds().len(), 1);
    assert_eq!(structure.bonds()[0].atom1_idx, 0);
    assert_eq!(structure.bonds()[0].atom2_idx, 1);
}

#[test]
fn test_parse_cif_multi_model_by_pdbx_model_num() {
    let mut structure = parse_cif(SAMPLE_MULTI_MODEL_CIF).expect("parse multi-model CIF");
    assert_eq!(structure.model_serials(), vec![1, 2, 5]);
    assert_eq!(structure.active_model_serial(), 1);
    assert!((structure.atoms()[0].pos.x - 0.0).abs() < 1e-4);

    structure.next_model();
    assert_eq!(structure.active_model_serial(), 2);
    assert!((structure.atoms()[0].pos.x - 5.0).abs() < 1e-4);

    structure.next_model();
    assert_eq!(structure.active_model_serial(), 5);
    structure.next_model();
    assert_eq!(structure.active_model_serial(), 1);

    structure.prev_model();
    assert_eq!(structure.active_model_serial(), 5);
}

#[test]
fn test_single_model_file_has_serial_one() {
    let structure = parse_pdb(SAMPLE_PDB_1CRN).unwrap();
    assert_eq!(structure.model_count(), 1);
    assert_eq!(structure.active_model_serial(), 1);
    assert!(!structure.has_multiple_models());
}
