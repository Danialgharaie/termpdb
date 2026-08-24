//! Malformed-input edge-case tests for the PDB and mmCIF parsers.
//!
//! These files exist to guarantee the parsers never panic on hostile or
//! corrupted input: multi-byte UTF-8 characters straddling fixed-column
//! boundaries, truncated records, non-numeric coordinates, and ragged CIF
//! loops must all produce either a valid structure or a `TermPdbError`
//! (ragged CIF loops specifically must be rejected, never silently truncated).

use std::sync::atomic::{AtomicUsize, Ordering};

use termpdb::error::TermPdbError;
use termpdb::parser::{load_structure, parse_cif, parse_pdb};

static TEMP_FILE_COUNTER: AtomicUsize = AtomicUsize::new(0);

/// Writes bytes to a uniquely named temp file (pid + counter suffix so
/// concurrent test runs cannot collide).
fn write_temp_file(extension: &str, contents: &[u8]) -> std::path::PathBuf {
    let n = TEMP_FILE_COUNTER.fetch_add(1, Ordering::SeqCst);
    let path = std::env::temp_dir().join(format!(
        "termpdb_edge_case_{}_{n}.{extension}",
        std::process::id()
    ));
    std::fs::write(&path, contents).expect("write temp file");
    path
}

/// A minimal valid 1CRN-style ATOM record used as the base for mutations.
const ATOM_LINE: &str =
    "ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N";

#[test]
fn test_multibyte_atom_name_within_field_width_parses() {
    // "é" is two bytes; padded with spaces the atom-name field still spans
    // exactly four bytes, so every downstream fixed-column slice stays
    // aligned. With no element column present the parser must fall back to
    // name-based inference and report an unknown element instead of panicking
    // on a non-char-boundary slice.
    let base_line = &ATOM_LINE[..66]; // drop the element/charge columns
    let pdb = base_line.replacen(" N  ", "é  ", 1);
    let structure = parse_pdb(&pdb).expect("multibyte atom name should not fail parsing");
    assert_eq!(structure.atom_count(), 1);
    assert_eq!(structure.atoms()[0].element.atomic_number, 0);
}

#[test]
fn test_multibyte_char_straddling_record_type_skips_line() {
    // "Ö" at byte offset 2 makes safe_slice(line, 0, 6) land inside a
    // multi-byte character. The record type becomes unrecognizable and the
    // line is skipped; surrounding records still parse.
    let bad_line = ATOM_LINE.replacen("ATOM", "ATÖM", 1);
    let pdb = format!("{bad_line}\n{ATOM_LINE}\nEND\n");
    let structure = parse_pdb(&pdb).expect("straddled record line must not panic");
    assert_eq!(structure.atom_count(), 1, "only the intact line parses");
}

#[test]
fn test_multibyte_residue_name_returns_error_not_panic() {
    // "É" widens the residue-name field by one byte, shifting every later
    // column; the coordinate slices then land mid-character and inside what
    // is no longer a number. Must yield ParseError, never a panic.
    let bad_line = ATOM_LINE.replacen("THR", "THÉ", 1);
    let result = parse_pdb(&format!("{bad_line}\nEND\n"));
    match result {
        Err(TermPdbError::ParseError(_)) => {}
        Ok(s) => assert_eq!(s.atom_count(), 0, "no atom can be salvaged"),
        Err(other) => panic!("unexpected error kind: {other}"),
    }
}

#[test]
fn test_multibyte_helix_chain_id_does_not_panic() {
    // HELIX chain/seq columns contain a multi-byte char; the fixed-column
    // read fails and the whitespace-token fallback recovers the record.
    let helix = "HELIX    1  HÏ THR A    7  GLY A   17  1";
    let pdb = format!("{}\n{}\nEND\n", helix, ATOM_LINE);
    let result = parse_pdb(&pdb);
    assert!(result.is_ok(), "multibyte HELIX record must not panic");
}

#[test]
fn test_truncated_atom_line_is_parse_error() {
    // Line cut off right after the X coordinate: the Y column is missing
    // entirely, which must yield a ParseError naming the Y field.
    let truncated = &ATOM_LINE[..38];
    let err = parse_pdb(truncated).expect_err("truncated coordinates must be an error");
    assert!(matches!(err, TermPdbError::ParseError(msg) if msg.contains("Y")));
}

#[test]
fn test_non_numeric_coordinates_are_parse_error() {
    // Same byte width as the original X field, so only X is corrupted.
    let garbage_x = ATOM_LINE.replacen("  17.047", "abcdefgh", 1);
    let err = parse_pdb(&garbage_x).expect_err("non-numeric X must be an error");
    assert!(matches!(err, TermPdbError::ParseError(msg) if msg.contains("X")));
}

#[test]
fn test_empty_input_parses_to_empty_structure() {
    // A Structure always carries at least one (possibly empty) model, so an
    // empty file yields a single default model with zero atoms.
    for empty in ["", "\n", "   \n  \n"] {
        let structure = parse_pdb(empty).expect("empty input must parse");
        assert_eq!(structure.atom_count(), 0);
        assert_eq!(structure.chain_count(), 0);
        assert_eq!(structure.model_count(), 1);
    }
}

#[test]
fn test_binary_garbage_text_does_not_panic() {
    // Valid UTF-8 full of control characters, emoji, and RTL marks: no
    // recognizable records, but nothing may panic.
    let garbage = "\u{7f}\u{200f}\u{10ffff}\n\u{1F600}\u{0b}\ndata_??? \u{80}\n";
    let structure = parse_pdb(garbage).expect("garbage lines are simply skipped");
    assert_eq!(structure.atom_count(), 0);
}

#[test]
fn test_cif_multibyte_atom_name_parses_with_unknown_element() {
    // 10 headers; each data row carries exactly 10 values.
    let cif = r#"data_EDGE
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N CA ALA A 1 0.0 0.0 0.0
HETATM 2 ? é HOH A 2 1.0 1.0 1.0
"#;
    let structure = parse_cif(cif).expect("multibyte CIF atom name must parse");
    assert_eq!(structure.atom_count(), 2);
    assert_eq!(structure.atoms()[0].element.symbol, "N");
    // The multi-byte atom name falls through name-based inference to unknown,
    // exercising the char-boundary-safe slicing (previously a panic).
    assert_eq!(structure.atoms()[1].element.atomic_number, 0);
}

#[test]
fn test_ragged_cif_loop_is_rejected() {
    // 10 headers; 23 values (two complete rows plus a 3-value orphan tail).
    // Integer division would silently drop the orphan tokens and shift every
    // downstream row, so the parser must reject the loop outright instead of
    // truncating.
    let cif = r#"data_RAGGED
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N CA ALA A 1 0.0 0.0 0.0
ATOM 2 C CB ALA A 2 1.0 1.0 1.0
ATOM 3 O
"#;
    let err = parse_cif(cif).expect_err("ragged loop must be an error, not a silent truncation");
    match err {
        TermPdbError::ParseError(msg) => {
            assert!(
                msg.contains("_atom_site.group_PDB"),
                "error must name the first header tag: {msg}"
            );
            assert!(
                msg.contains("10"),
                "error must name the column count: {msg}"
            );
            assert!(msg.contains("23"), "error must name the value count: {msg}");
        }
        other => panic!("unexpected error kind: {other:?}"),
    }
}

#[test]
fn test_rectangular_cif_loop_still_parses() {
    // The hardening must not disturb well-formed loops: three complete
    // 10-column rows yield three atoms across two chains.
    let cif = r#"data_OK
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N CA ALA A 1 0.0 0.0 0.0
ATOM 2 C CB ALA A 2 1.0 1.0 1.0
HETATM 3 O O HOH B 1 2.0 2.0 2.0
"#;
    let structure = parse_cif(cif).expect("rectangular loop must parse");
    assert_eq!(structure.atom_count(), 3);
}

#[test]
fn test_stray_data_word_mid_values_preserves_row_count() {
    // `data_...` is a data-block marker only at column 0. Mid-line it is an
    // ordinary value token, so its appearance inside loop values (here as the
    // `_atom_site.id` cell of row 2) must neither terminate value collection
    // early nor drop rows: previously everything from that word onward was
    // discarded, leaving one atom instead of two.
    let cif = r#"data_MAIN
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 N CA ALA A 1 0.0 0.0 0.0
ATOM data_foo C CB ALA B 2 1.0 1.0 1.0
"#;
    let structure = parse_cif(cif).expect("stray mid-line data_ word must not truncate the loop");
    assert_eq!(structure.atom_count(), 2);
}

#[test]
fn test_load_structure_rejects_invalid_utf8_file() {
    let path = write_temp_file("pdb", &[0xFF, 0xFE, b'A', b'T', b'O', b'M']);
    let result = load_structure(path.to_str().unwrap());
    std::fs::remove_file(&path).ok();
    assert!(
        matches!(result, Err(TermPdbError::ParseError(_))),
        "invalid UTF-8 must surface as ParseError"
    );
}
