use std::collections::HashMap;
use termpdb::math::{Mat4, Vec3};
use termpdb::model::assembly::{assembly_transforms, parse_oper_expression};
use termpdb::parser::{parse_cif, parse_pdb};

const DIMER_PDB: &str = r#"HEADER    TEST DIMER                              01-JAN-24   1DIM
TITLE     SYNTHETIC DIMER VIA BIOMT TRANSLATION
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2  1.000000  0.000000  0.000000       10.00000
REMARK 350   BIOMT2   2  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00 10.00           C
END
"#;

const DIMER_CIF: &str = r#"data_1DIM
_entry.id 1DIM
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
1 dimer
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1,2 A
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 1 0 0 0 1 0 0 0 1 0 0 0
2 1 0 0 0 1 0 0 0 1 10 0 0
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
_atom_site.auth_asym_id
ATOM 1 N N ALA A 1 0.0 0.0 0.0 A
ATOM 2 C CA ALA A 1 1.5 0.0 0.0 A
#
"#;

#[test]
fn test_parse_oper_expression_groups_and_ranges() {
    assert_eq!(
        parse_oper_expression("1,2").unwrap(),
        vec![vec!["1".to_string(), "2".to_string()]]
    );
    assert_eq!(
        parse_oper_expression("(1-3)").unwrap(),
        vec![vec!["1".to_string(), "2".to_string(), "3".to_string()]]
    );
    assert_eq!(
        parse_oper_expression("(1,2)(3)").unwrap(),
        vec![
            vec!["1".to_string(), "2".to_string()],
            vec!["3".to_string()]
        ]
    );
}

#[test]
fn test_operator_product_composes_translation() {
    let mut ops = HashMap::new();
    ops.insert("1".into(), Mat4::identity());
    ops.insert("2".into(), Mat4::from_translation(Vec3::new(5.0, 0.0, 0.0)));
    let xf = assembly_transforms("(1)(2)", &ops).unwrap();
    assert_eq!(xf.len(), 1);
    let p = xf[0].1.transform_point(Vec3::new(1.0, 0.0, 0.0));
    assert!((p.x - 6.0).abs() < 1e-4);
}

#[test]
fn test_pdb_assembly_defaults_to_asu_and_expands_dimer() {
    let mut s = parse_pdb(DIMER_PDB).expect("parse dimer pdb");
    assert!(s.has_assemblies());
    assert_eq!(s.assembly_ids(), vec!["1"]);
    assert_eq!(s.active_assembly_id(), None);
    assert_eq!(s.atom_count(), 2);
    assert!((s.atoms()[1].pos.x - 1.5).abs() < 1e-4);

    s.set_assembly(Some("1")).unwrap();
    assert_eq!(s.active_assembly_id(), Some("1"));
    assert_eq!(s.atom_count(), 4);
    assert_eq!(s.chain_count(), 2);
    let xs: Vec<f32> = s
        .atoms()
        .iter()
        .filter(|a| a.name == "CA")
        .map(|a| a.pos.x)
        .collect();
    assert!(xs.iter().any(|x| (x - 1.5).abs() < 1e-3), "{xs:?}");
    assert!(xs.iter().any(|x| (x - 11.5).abs() < 1e-3), "{xs:?}");

    s.set_assembly(Some("asu")).unwrap();
    assert_eq!(s.atom_count(), 2);
    assert!(s.set_assembly(Some("9")).is_err());
}

#[test]
fn test_cif_assembly_expands_the_same_way() {
    let mut s = parse_cif(DIMER_CIF).expect("parse dimer cif");
    assert_eq!(s.assembly_ids(), vec!["1"]);
    assert_eq!(s.atom_count(), 2);
    s.set_assembly(Some("1")).unwrap();
    assert_eq!(s.atom_count(), 4);
    let xs: Vec<f32> = s
        .atoms()
        .iter()
        .filter(|a| a.name == "CA")
        .map(|a| a.pos.x)
        .collect();
    assert!(xs.iter().any(|x| (x - 11.5).abs() < 1e-3), "{xs:?}");
}

#[test]
fn test_next_assembly_wraps_through_asu() {
    let mut s = parse_pdb(DIMER_PDB).unwrap();
    s.next_assembly();
    assert_eq!(s.active_assembly_id(), Some("1"));
    assert_eq!(s.atom_count(), 4);
    s.next_assembly();
    assert_eq!(s.active_assembly_id(), None);
    assert_eq!(s.atom_count(), 2);
    s.prev_assembly();
    assert_eq!(s.active_assembly_id(), Some("1"));
}
