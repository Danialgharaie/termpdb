//! Periodic table elements and chemical metadata.
//!
//! Provides element lookup by atomic symbol or atomic number, with standard
//! covalent radii, Van der Waals radii, and CPK RGB color definitions.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Element {
    /// Atomic symbol (e.g. "C", "N", "O", "Fe")
    pub symbol: &'static str,
    /// Full element name (e.g. "Carbon")
    pub name: &'static str,
    /// Atomic number (Z)
    pub atomic_number: u8,
    /// Covalent radius in Angstroms (Å)
    pub covalent_radius: f32,
    /// Van der Waals radius in Angstroms (Å)
    pub vdw_radius: f32,
    /// Default CPK color as (R, G, B) tuple
    pub cpk_color: (u8, u8, u8),
}

impl Element {
    /// Fallback element for unknown or unassigned atom types.
    pub const fn unknown() -> Self {
        Self {
            symbol: "X",
            name: "Unknown",
            atomic_number: 0,
            covalent_radius: 1.50,
            vdw_radius: 1.70,
            cpk_color: (200, 200, 200),
        }
    }
}

impl Default for Element {
    fn default() -> Self {
        Self::unknown()
    }
}

/// Static table of common chemical and biomolecular elements.
pub static ELEMENTS: &[Element] = &[
    Element {
        symbol: "H",
        name: "Hydrogen",
        atomic_number: 1,
        covalent_radius: 0.31,
        vdw_radius: 1.20,
        cpk_color: (255, 255, 255),
    },
    Element {
        symbol: "He",
        name: "Helium",
        atomic_number: 2,
        covalent_radius: 0.28,
        vdw_radius: 1.40,
        cpk_color: (217, 255, 255),
    },
    Element {
        symbol: "Li",
        name: "Lithium",
        atomic_number: 3,
        covalent_radius: 1.28,
        vdw_radius: 1.82,
        cpk_color: (204, 128, 255),
    },
    Element {
        symbol: "Be",
        name: "Beryllium",
        atomic_number: 4,
        covalent_radius: 0.96,
        vdw_radius: 1.53,
        cpk_color: (194, 255, 0),
    },
    Element {
        symbol: "B",
        name: "Boron",
        atomic_number: 5,
        covalent_radius: 0.84,
        vdw_radius: 1.92,
        cpk_color: (255, 181, 181),
    },
    Element {
        symbol: "C",
        name: "Carbon",
        atomic_number: 6,
        covalent_radius: 0.76,
        vdw_radius: 1.70,
        cpk_color: (144, 144, 144),
    },
    Element {
        symbol: "N",
        name: "Nitrogen",
        atomic_number: 7,
        covalent_radius: 0.71,
        vdw_radius: 1.55,
        cpk_color: (48, 80, 248),
    },
    Element {
        symbol: "O",
        name: "Oxygen",
        atomic_number: 8,
        covalent_radius: 0.66,
        vdw_radius: 1.52,
        cpk_color: (255, 13, 13),
    },
    Element {
        symbol: "F",
        name: "Fluorine",
        atomic_number: 9,
        covalent_radius: 0.57,
        vdw_radius: 1.47,
        cpk_color: (144, 224, 80),
    },
    Element {
        symbol: "Ne",
        name: "Neon",
        atomic_number: 10,
        covalent_radius: 0.58,
        vdw_radius: 1.54,
        cpk_color: (179, 227, 245),
    },
    Element {
        symbol: "Na",
        name: "Sodium",
        atomic_number: 11,
        covalent_radius: 1.66,
        vdw_radius: 2.27,
        cpk_color: (171, 92, 242),
    },
    Element {
        symbol: "Mg",
        name: "Magnesium",
        atomic_number: 12,
        covalent_radius: 1.41,
        vdw_radius: 1.73,
        cpk_color: (138, 255, 0),
    },
    Element {
        symbol: "Al",
        name: "Aluminium",
        atomic_number: 13,
        covalent_radius: 1.21,
        vdw_radius: 1.84,
        cpk_color: (191, 166, 166),
    },
    Element {
        symbol: "Si",
        name: "Silicon",
        atomic_number: 14,
        covalent_radius: 1.11,
        vdw_radius: 2.10,
        cpk_color: (240, 200, 160),
    },
    Element {
        symbol: "P",
        name: "Phosphorus",
        atomic_number: 15,
        covalent_radius: 1.07,
        vdw_radius: 1.80,
        cpk_color: (255, 128, 0),
    },
    Element {
        symbol: "S",
        name: "Sulfur",
        atomic_number: 16,
        covalent_radius: 1.05,
        vdw_radius: 1.80,
        cpk_color: (255, 255, 48),
    },
    Element {
        symbol: "Cl",
        name: "Chlorine",
        atomic_number: 17,
        covalent_radius: 1.02,
        vdw_radius: 1.75,
        cpk_color: (31, 240, 31),
    },
    Element {
        symbol: "Ar",
        name: "Argon",
        atomic_number: 18,
        covalent_radius: 1.06,
        vdw_radius: 1.88,
        cpk_color: (128, 209, 227),
    },
    Element {
        symbol: "K",
        name: "Potassium",
        atomic_number: 19,
        covalent_radius: 2.03,
        vdw_radius: 2.75,
        cpk_color: (143, 64, 212),
    },
    Element {
        symbol: "Ca",
        name: "Calcium",
        atomic_number: 20,
        covalent_radius: 1.76,
        vdw_radius: 2.00,
        cpk_color: (61, 255, 0),
    },
    Element {
        symbol: "Sc",
        name: "Scandium",
        atomic_number: 21,
        covalent_radius: 1.70,
        vdw_radius: 2.11,
        cpk_color: (230, 230, 230),
    },
    Element {
        symbol: "Ti",
        name: "Titanium",
        atomic_number: 22,
        covalent_radius: 1.60,
        vdw_radius: 2.00,
        cpk_color: (191, 194, 199),
    },
    Element {
        symbol: "V",
        name: "Vanadium",
        atomic_number: 23,
        covalent_radius: 1.53,
        vdw_radius: 2.00,
        cpk_color: (166, 166, 171),
    },
    Element {
        symbol: "Cr",
        name: "Chromium",
        atomic_number: 24,
        covalent_radius: 1.39,
        vdw_radius: 2.00,
        cpk_color: (138, 153, 199),
    },
    Element {
        symbol: "Mn",
        name: "Manganese",
        atomic_number: 25,
        covalent_radius: 1.39,
        vdw_radius: 1.73,
        cpk_color: (156, 122, 199),
    },
    Element {
        symbol: "Fe",
        name: "Iron",
        atomic_number: 26,
        covalent_radius: 1.32,
        vdw_radius: 1.80,
        cpk_color: (224, 102, 51),
    },
    Element {
        symbol: "Co",
        name: "Cobalt",
        atomic_number: 27,
        covalent_radius: 1.26,
        vdw_radius: 1.70,
        cpk_color: (240, 144, 160),
    },
    Element {
        symbol: "Ni",
        name: "Nickel",
        atomic_number: 28,
        covalent_radius: 1.24,
        vdw_radius: 1.63,
        cpk_color: (80, 208, 80),
    },
    Element {
        symbol: "Cu",
        name: "Copper",
        atomic_number: 29,
        covalent_radius: 1.32,
        vdw_radius: 1.40,
        cpk_color: (200, 115, 51),
    },
    Element {
        symbol: "Zn",
        name: "Zinc",
        atomic_number: 30,
        covalent_radius: 1.22,
        vdw_radius: 1.39,
        cpk_color: (125, 128, 168),
    },
    Element {
        symbol: "Ga",
        name: "Gallium",
        atomic_number: 31,
        covalent_radius: 1.22,
        vdw_radius: 1.87,
        cpk_color: (194, 143, 143),
    },
    Element {
        symbol: "Ge",
        name: "Germanium",
        atomic_number: 32,
        covalent_radius: 1.20,
        vdw_radius: 2.11,
        cpk_color: (102, 143, 143),
    },
    Element {
        symbol: "As",
        name: "Arsenic",
        atomic_number: 33,
        covalent_radius: 1.19,
        vdw_radius: 1.85,
        cpk_color: (189, 128, 227),
    },
    Element {
        symbol: "Se",
        name: "Selenium",
        atomic_number: 34,
        covalent_radius: 1.20,
        vdw_radius: 1.90,
        cpk_color: (255, 161, 0),
    },
    Element {
        symbol: "Br",
        name: "Bromine",
        atomic_number: 35,
        covalent_radius: 1.20,
        vdw_radius: 1.85,
        cpk_color: (166, 41, 41),
    },
    Element {
        symbol: "Kr",
        name: "Krypton",
        atomic_number: 36,
        covalent_radius: 1.16,
        vdw_radius: 2.02,
        cpk_color: (92, 184, 209),
    },
    Element {
        symbol: "Rb",
        name: "Rubidium",
        atomic_number: 37,
        covalent_radius: 2.20,
        vdw_radius: 3.03,
        cpk_color: (112, 46, 176),
    },
    Element {
        symbol: "Sr",
        name: "Strontium",
        atomic_number: 38,
        covalent_radius: 1.95,
        vdw_radius: 2.49,
        cpk_color: (0, 255, 0),
    },
    Element {
        symbol: "Y",
        name: "Yttrium",
        atomic_number: 39,
        covalent_radius: 1.90,
        vdw_radius: 2.40,
        cpk_color: (148, 255, 255),
    },
    Element {
        symbol: "Zr",
        name: "Zirconium",
        atomic_number: 40,
        covalent_radius: 1.75,
        vdw_radius: 2.30,
        cpk_color: (224, 224, 224),
    },
    Element {
        symbol: "Mo",
        name: "Molybdenum",
        atomic_number: 42,
        covalent_radius: 1.54,
        vdw_radius: 2.10,
        cpk_color: (84, 181, 181),
    },
    Element {
        symbol: "Ru",
        name: "Ruthenium",
        atomic_number: 44,
        covalent_radius: 1.46,
        vdw_radius: 2.05,
        cpk_color: (36, 158, 143),
    },
    Element {
        symbol: "Rh",
        name: "Rhodium",
        atomic_number: 45,
        covalent_radius: 1.42,
        vdw_radius: 2.00,
        cpk_color: (10, 125, 140),
    },
    Element {
        symbol: "Pd",
        name: "Palladium",
        atomic_number: 46,
        covalent_radius: 1.39,
        vdw_radius: 1.63,
        cpk_color: (0, 105, 133),
    },
    Element {
        symbol: "Ag",
        name: "Silver",
        atomic_number: 47,
        covalent_radius: 1.45,
        vdw_radius: 1.72,
        cpk_color: (192, 192, 192),
    },
    Element {
        symbol: "Cd",
        name: "Cadmium",
        atomic_number: 48,
        covalent_radius: 1.44,
        vdw_radius: 1.58,
        cpk_color: (255, 217, 143),
    },
    Element {
        symbol: "In",
        name: "Indium",
        atomic_number: 49,
        covalent_radius: 1.42,
        vdw_radius: 1.93,
        cpk_color: (166, 117, 115),
    },
    Element {
        symbol: "Sn",
        name: "Tin",
        atomic_number: 50,
        covalent_radius: 1.39,
        vdw_radius: 2.17,
        cpk_color: (102, 128, 128),
    },
    Element {
        symbol: "I",
        name: "Iodine",
        atomic_number: 53,
        covalent_radius: 1.39,
        vdw_radius: 1.98,
        cpk_color: (148, 0, 148),
    },
    Element {
        symbol: "Cs",
        name: "Caesium",
        atomic_number: 55,
        covalent_radius: 2.44,
        vdw_radius: 3.43,
        cpk_color: (87, 23, 143),
    },
    Element {
        symbol: "Ba",
        name: "Barium",
        atomic_number: 56,
        covalent_radius: 2.15,
        vdw_radius: 2.68,
        cpk_color: (0, 201, 0),
    },
    Element {
        symbol: "La",
        name: "Lanthanum",
        atomic_number: 57,
        covalent_radius: 2.07,
        vdw_radius: 2.50,
        cpk_color: (112, 212, 255),
    },
    Element {
        symbol: "Ce",
        name: "Cerium",
        atomic_number: 58,
        covalent_radius: 2.04,
        vdw_radius: 2.48,
        cpk_color: (255, 255, 199),
    },
    Element {
        symbol: "W",
        name: "Tungsten",
        atomic_number: 74,
        covalent_radius: 1.62,
        vdw_radius: 2.10,
        cpk_color: (33, 148, 214),
    },
    Element {
        symbol: "Os",
        name: "Osmium",
        atomic_number: 76,
        covalent_radius: 1.44,
        vdw_radius: 2.00,
        cpk_color: (38, 102, 150),
    },
    Element {
        symbol: "Ir",
        name: "Iridium",
        atomic_number: 77,
        covalent_radius: 1.41,
        vdw_radius: 2.00,
        cpk_color: (23, 84, 135),
    },
    Element {
        symbol: "Pt",
        name: "Platinum",
        atomic_number: 78,
        covalent_radius: 1.36,
        vdw_radius: 1.75,
        cpk_color: (208, 208, 224),
    },
    Element {
        symbol: "Au",
        name: "Gold",
        atomic_number: 79,
        covalent_radius: 1.36,
        vdw_radius: 1.66,
        cpk_color: (255, 209, 35),
    },
    Element {
        symbol: "Hg",
        name: "Mercury",
        atomic_number: 80,
        covalent_radius: 1.32,
        vdw_radius: 1.55,
        cpk_color: (184, 184, 208),
    },
    Element {
        symbol: "Pb",
        name: "Lead",
        atomic_number: 82,
        covalent_radius: 1.46,
        vdw_radius: 2.02,
        cpk_color: (87, 89, 97),
    },
    Element {
        symbol: "U",
        name: "Uranium",
        atomic_number: 92,
        covalent_radius: 1.96,
        vdw_radius: 1.86,
        cpk_color: (0, 143, 255),
    },
];

/// Looks up an element by its symbol (case-insensitive, trims whitespace).
/// Returns [`Element::unknown()`] if not found.
pub fn element_by_symbol(symbol: &str) -> Element {
    let trimmed = symbol.trim();
    if trimmed.is_empty() {
        return Element::unknown();
    }

    for elem in ELEMENTS {
        if elem.symbol.eq_ignore_ascii_case(trimmed) {
            return *elem;
        }
    }

    Element::unknown()
}

/// Looks up an element by atomic number (Z).
/// Returns [`Element::unknown()`] if not found.
pub fn element_by_atomic_number(z: u8) -> Element {
    if z == 0 {
        return Element::unknown();
    }

    for elem in ELEMENTS {
        if elem.atomic_number == z {
            return *elem;
        }
    }

    Element::unknown()
}
