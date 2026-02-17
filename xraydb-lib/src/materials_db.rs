/// Embedded materials database (from XrayDB materials.dat).
///
/// Each entry: (name, density_g_per_cm3, formula)
pub(crate) const MATERIALS: &[(&str, f64, &str)] = &[
    // Gases
    ("hydrogen", 0.0000899, "H"),
    ("helium", 0.0001786, "He"),
    ("nitrogen", 0.00125, "N"),
    ("oxygen", 0.001429, "O"),
    ("neon", 0.0009002, "Ne"),
    ("argon", 0.001784, "Ar"),
    ("krypton", 0.003749, "Kr"),
    ("xenon", 0.005894, "Xe"),
    (
        "air",
        0.001225,
        "(N2)0.7808(O2)0.2095Ar9.34e-3(CO2)4.1e-4Ne1.82e-5He5.24e-6(CH4)1.8e-6Kr1.0e-6(H2)0.5e-6Xe9.e-8",
    ),
    ("methane", 0.000657, "CH4"),
    ("carbon dioxide", 0.001562, "CO2"),
    // Solvents
    ("water", 1.0, "H2O"),
    ("ethanol", 0.789, "C2H5OH"),
    ("acetone", 0.785, "C3H6O"),
    ("methanol", 0.791, "CH3OH"),
    ("isopropanol", 0.803, "C3H8O"),
    ("toluene", 0.867, "C7H8"),
    ("xylene", 0.844, "C6H4(CH3)2"),
    ("benzene", 0.877, "C6H6"),
    ("butanol", 0.810, "C4H10O"),
    ("chlorobenzene", 1.106, "C6H5Cl"),
    ("cyclohexane", 0.774, "C6H12"),
    ("dimethyl sulfoxide", 1.09, "C2H6OS"),
    ("ethylene glycol", 1.115, "C2H6O2"),
    ("glycerin", 1.261, "C3H8O3"),
    ("heptane", 0.684, "C7H16"),
    ("hexane", 0.659, "C6H14"),
    // Polymers
    ("kapton", 1.42, "C22H10N2O5"),
    ("polyimide", 1.42, "C22H10N2O5"),
    ("polypropylene", 0.86, "C3H6"),
    ("pmma", 1.18, "C5H8O2"),
    ("polycarbonate", 1.2, "C16H14O3"),
    ("kimol", 1.2, "C16H14O3"),
    ("mylar", 1.4, "C10H8O4"),
    ("teflon", 2.2, "C2F4"),
    ("parylene-c", 1.29, "C8H7Cl"),
    ("parylene-n", 1.11, "C8H8"),
    ("peek", 1.32, "C19H14O3"),
    // Ceramics & minerals
    ("boron nitride", 2.1, "BN"),
    ("cubic boron nitride", 3.45, "BN"),
    ("silicon nitride", 3.17, "Si3N4"),
    ("yag", 4.56, "Y3Al5O12"),
    ("sapphire", 4.0, "Al2O3"),
    ("ule glass", 2.205, "Si0.925Ti0.075O2"),
    (
        "zerodur",
        2.53,
        "Si0.56Al0.5P0.16Li0.04Ti0.02Zr0.02Zn0.03O2.46",
    ),
    ("fluorite", 3.18, "CaF2"),
    ("mica", 2.83, "KAl3Si3O12H2"),
    ("fayalite", 4.392, "Fe2SiO4"),
    ("forsterite", 3.27, "Mg2SiO4"),
    ("wustite", 5.7, "FeO"),
    ("salt", 2.165, "NaCl"),
    ("silica", 2.2, "SiO2"),
    ("quartz", 2.65, "SiO2"),
    ("cristobalite", 2.27, "SiO2"),
    ("rutile", 4.23, "TiO2"),
    ("magnesium dioxide", 3.6, "MgO"),
    ("galena", 7.60, "PbS"),
    // Semiconductors
    ("cadmium telluride", 5.85, "CdTe"),
    ("gallium arsenide", 5.318, "GaAs"),
    // Metals & elements
    ("beryllium copper", 8.4, "Cu0.98Be0.02"),
    ("diamond carbon", 3.52, "C"),
    ("graphite carbon", 2.23, "C"),
    ("beryllium", 1.85, "Be"),
    ("aluminum", 2.70, "Al"),
    ("silicon", 2.329, "Si"),
    ("titanium", 4.506, "Ti"),
    ("chromium", 7.15, "Cr"),
    ("iron", 7.88, "Fe"),
    ("cobalt", 8.90, "Co"),
    ("nickel", 8.908, "Ni"),
    ("copper", 8.96, "Cu"),
    ("zinc", 7.14, "Zn"),
    ("gallium", 5.91, "Ga"),
    ("germanium", 5.323, "Ge"),
    ("molybdenum", 10.28, "Mo"),
    ("ruthenium", 12.45, "Ru"),
    ("rhodium", 12.41, "Rh"),
    ("palladium", 12.02, "Pd"),
    ("silver", 10.49, "Ag"),
    ("indium", 7.31, "In"),
    ("tin", 7.265, "Sn"),
    ("tantalum", 16.69, "Ta"),
    ("tungsten", 19.25, "W"),
    ("rhenium", 21.02, "Re"),
    ("osmium", 22.59, "Os"),
    ("iridium", 22.56, "Ir"),
    ("platinum", 21.45, "Pt"),
    ("gold", 19.3, "Au"),
    ("mercury", 13.534, "Hg"),
    ("lead", 11.34, "Pb"),
    ("bismuth", 9.78, "Bi"),
    ("uranium", 19.1, "U"),
    ("zirconium", 6.5, "Zr"),
];

/// Find a material by name (case-insensitive) or formula.
/// Returns (formula, density).
pub(crate) fn find_material(name: &str) -> Option<(&'static str, f64)> {
    let lower = name.to_lowercase();
    // Try by name first
    for &(mat_name, density, formula) in MATERIALS {
        if mat_name == lower {
            return Some((formula, density));
        }
    }
    // Try by formula
    for &(_, density, formula) in MATERIALS {
        if formula.eq_ignore_ascii_case(name) {
            return Some((formula, density));
        }
    }
    None
}
