use termpdb::math::Vec3;
use termpdb::model::{Atom, Chain, Residue, Structure, element_by_symbol};
use termpdb::render::color::{ColorScheme, color_for_atom};

fn make_test_atom(res_name: &str, b_factor: f32, chain_id: &str) -> (Atom, Residue, Structure) {
    let mut structure = Structure::new("test");
    let mut res = Residue::new(1, res_name, chain_id);
    let atom = Atom::new(
        0,
        1,
        "CA",
        element_by_symbol("C"),
        Vec3::new(0.0, 0.0, 0.0),
        b_factor,
        res_name,
        1,
        chain_id,
        false,
    );
    res.atom_indices.push(0);
    let mut chain = Chain::new(chain_id);
    chain.residues.push(res.clone());
    structure.add_chain(chain);
    structure.add_atom(atom.clone());
    (atom, res, structure)
}

#[test]
fn test_plddt_color_scheme_thresholds() {
    // > 90: Deep Blue
    let (a_high, r_high, s_high) = make_test_atom("ALA", 95.0, "A");
    let c_high = color_for_atom(&a_high, Some(&r_high), &s_high, ColorScheme::Plddt);
    assert_eq!(c_high, (0, 83, 214));

    // 70 - 90: Cyan / Light Blue
    let (a_conf, r_conf, s_conf) = make_test_atom("ALA", 80.0, "A");
    let c_conf = color_for_atom(&a_conf, Some(&r_conf), &s_conf, ColorScheme::Plddt);
    assert_eq!(c_conf, (101, 203, 243));

    // 50 - 70: Yellow
    let (a_low, r_low, s_low) = make_test_atom("ALA", 60.0, "A");
    let c_low = color_for_atom(&a_low, Some(&r_low), &s_low, ColorScheme::Plddt);
    assert_eq!(c_low, (255, 219, 19));

    // < 50: Orange
    let (a_vlow, r_vlow, s_vlow) = make_test_atom("ALA", 35.0, "A");
    let c_vlow = color_for_atom(&a_vlow, Some(&r_vlow), &s_vlow, ColorScheme::Plddt);
    assert_eq!(c_vlow, (255, 125, 69));
}

#[test]
fn test_electrostatic_color_scheme() {
    // Acidic / Negative: ASP -> Red
    let (a_asp, r_asp, s_asp) = make_test_atom("ASP", 20.0, "A");
    let c_asp = color_for_atom(&a_asp, Some(&r_asp), &s_asp, ColorScheme::Electrostatic);
    assert!(c_asp.0 > 200 && c_asp.1 < 100 && c_asp.2 < 100);

    // Basic / Positive: LYS -> Blue
    let (a_lys, r_lys, s_lys) = make_test_atom("LYS", 20.0, "A");
    let c_lys = color_for_atom(&a_lys, Some(&r_lys), &s_lys, ColorScheme::Electrostatic);
    assert!(c_lys.2 > 200 && c_lys.0 < 100);
}

#[test]
fn test_theme_color_schemes() {
    let (a, r, s) = make_test_atom("ALA", 20.0, "A");
    let catppuccin = color_for_atom(&a, Some(&r), &s, ColorScheme::Catppuccin);
    let nord = color_for_atom(&a, Some(&r), &s, ColorScheme::Nord);
    let tokyo = color_for_atom(&a, Some(&r), &s, ColorScheme::TokyoNight);
    let gruvbox = color_for_atom(&a, Some(&r), &s, ColorScheme::Gruvbox);

    assert_ne!(catppuccin, (0, 0, 0));
    assert_ne!(nord, (0, 0, 0));
    assert_ne!(tokyo, (0, 0, 0));
    assert_ne!(gruvbox, (0, 0, 0));
}
