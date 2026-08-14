//! WASM bindings for [`xraydb`].
//!
//! Build with:
//!
//! ```sh
//! wasm-pack build --target web --release xraydb-wasm --out-dir ../web/pkg
//! ```
//!
//! Call [`init`] once after loading the module to install a panic hook that turns any
//! Rust panic into a readable JavaScript error.

#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(clippy::unwrap_used, clippy::expect_used)]

use wasm_bindgen::prelude::*;

use xraydb::{ChantlerKind, CrossSectionKind, XrayDb};

/// Install a panic hook so Rust panics surface as readable JS errors.
///
/// Safe to call more than once.
#[wasm_bindgen]
pub fn init() {
    console_error_panic_hook::set_once();
}

fn db() -> Result<XrayDb, JsError> {
    XrayDb::try_new().map_err(to_js)
}

fn to_js(e: xraydb::XrayDbError) -> JsError {
    JsError::new(&e.to_string())
}

fn parse_kind(kind: &str) -> Result<CrossSectionKind, JsError> {
    match kind.to_lowercase().as_str() {
        "total" => Ok(CrossSectionKind::Total),
        "photo" => Ok(CrossSectionKind::Photo),
        "coherent" | "coh" => Ok(CrossSectionKind::Coherent),
        "incoherent" | "incoh" => Ok(CrossSectionKind::Incoherent),
        _ => Err(JsError::new(&format!(
            "unknown cross-section kind '{kind}'; expected total, photo, coherent, or incoherent"
        ))),
    }
}

fn parse_chantler_kind(kind: &str) -> Result<ChantlerKind, JsError> {
    match kind.to_lowercase().as_str() {
        "total" => Ok(ChantlerKind::Total),
        "photo" => Ok(ChantlerKind::Photo),
        "incoherent" | "incoh" => Ok(ChantlerKind::Incoherent),
        _ => Err(JsError::new(&format!(
            "unknown Chantler kind '{kind}'; expected total, photo, or incoherent"
        ))),
    }
}

// ── Exported value types ──────────────────────────────────────────────────

/// Basic facts about an element.
#[wasm_bindgen(getter_with_clone)]
pub struct Element {
    /// Atomic number.
    pub z: u16,
    /// Element symbol, e.g. `"Fe"`.
    pub symbol: String,
    /// Element name, e.g. `"iron"`.
    pub name: String,
    /// Molar mass in g/mol.
    pub molar_mass: f64,
    /// Density in g/cm³.
    pub density: f64,
}

/// One X-ray absorption edge.
#[wasm_bindgen(getter_with_clone)]
pub struct Edge {
    /// IUPAC label, e.g. `"K"`.
    pub label: String,
    /// Edge energy in eV.
    pub energy: f64,
    /// Fluorescence yield.
    pub fluorescence_yield: f64,
    /// Absorption jump ratio.
    pub jump_ratio: f64,
}

/// One X-ray emission line.
#[derive(Clone)]
#[wasm_bindgen(getter_with_clone)]
pub struct EmissionLine {
    /// Siegbahn label, e.g. `"Ka1"`.
    pub label: String,
    /// Emission energy in eV.
    pub energy: f64,
    /// Relative intensity.
    pub intensity: f64,
    /// Initial (core-hole) level.
    pub initial_level: String,
    /// Final level.
    pub final_level: String,
}

/// Emission lines sharing one initial level.
#[wasm_bindgen(getter_with_clone)]
pub struct LineGroup {
    /// IUPAC label of the initial level, e.g. `"K"`.
    pub level: String,
    /// Lines from that level, strongest first.
    pub lines: Vec<EmissionLine>,
}

/// Refractive index decrements: `n = 1 - delta - i*beta`.
#[wasm_bindgen]
pub struct DeltaBeta {
    /// Real decrement δ.
    pub delta: f64,
    /// Imaginary part β.
    pub beta: f64,
    /// 1/e attenuation length in cm.
    pub attenuation_length_cm: f64,
}

/// Compton scattering energies in eV.
#[wasm_bindgen]
pub struct Compton {
    /// Energy of a photon scattered through 90°.
    pub xray_90deg: f64,
    /// Mean scattered-photon energy.
    pub xray_mean: f64,
    /// Mean recoil-electron energy.
    pub electron_mean: f64,
}

/// A material from the built-in materials table.
#[wasm_bindgen(getter_with_clone)]
pub struct MaterialInfo {
    /// Common name.
    pub name: String,
    /// Chemical formula.
    pub formula: String,
    /// Density in g/cm³.
    pub density: f64,
}

/// A candidate edge match from [`guess_edge`].
#[wasm_bindgen(getter_with_clone)]
pub struct EdgeMatch {
    /// Element symbol.
    pub element: String,
    /// IUPAC edge label.
    pub edge: String,
    /// Tabulated edge energy in eV.
    pub energy: f64,
    /// Signed difference from the queried energy, in eV.
    pub difference: f64,
}

// ── Metadata ──────────────────────────────────────────────────────────────

/// Number of elements in the database.
#[wasm_bindgen]
pub fn element_count() -> Result<usize, JsError> {
    Ok(db()?.element_count())
}

/// Upstream XrayDB version tag the embedded data was built from.
#[wasm_bindgen]
pub fn data_version() -> Result<Option<String>, JsError> {
    Ok(db()?.data_version().map(str::to_string))
}

// ── Element lookups ───────────────────────────────────────────────────────

/// Atomic number for a symbol, name, or Z string.
#[wasm_bindgen]
pub fn atomic_number(element: &str) -> Result<u16, JsError> {
    db()?.atomic_number(element).map_err(to_js)
}

/// Canonical element symbol.
#[wasm_bindgen]
pub fn symbol(element: &str) -> Result<String, JsError> {
    db()?.symbol(element).map(str::to_string).map_err(to_js)
}

/// Full element name.
#[wasm_bindgen]
pub fn atomic_name(element: &str) -> Result<String, JsError> {
    db()?
        .atomic_name(element)
        .map(str::to_string)
        .map_err(to_js)
}

/// Molar mass in g/mol.
#[wasm_bindgen]
pub fn molar_mass(element: &str) -> Result<f64, JsError> {
    db()?.molar_mass(element).map_err(to_js)
}

/// Elemental density in g/cm³.
#[wasm_bindgen]
pub fn element_density(element: &str) -> Result<f64, JsError> {
    db()?.density(element).map_err(to_js)
}

/// All facts about an element in one call.
#[wasm_bindgen]
pub fn element(element: &str) -> Result<Element, JsError> {
    let db = db()?;
    Ok(Element {
        z: db.atomic_number(element).map_err(to_js)?,
        symbol: db.symbol(element).map_err(to_js)?.to_string(),
        name: db.atomic_name(element).map_err(to_js)?.to_string(),
        molar_mass: db.molar_mass(element).map_err(to_js)?,
        density: db.density(element).map_err(to_js)?,
    })
}

// ── Cross-sections ────────────────────────────────────────────────────────

/// Elam mass attenuation coefficients in cm²/g.
///
/// `kind` is one of `total`, `photo`, `coherent`, `incoherent`.
#[wasm_bindgen]
pub fn mu_elam(element: &str, energies: &[f64], kind: &str) -> Result<Vec<f64>, JsError> {
    let k = parse_kind(kind)?;
    db()?.mu_elam(element, energies, k).map_err(to_js)
}

/// f1 — real part of the anomalous scattering factor (Chantler).
#[wasm_bindgen]
pub fn f1_chantler(element: &str, energies: &[f64]) -> Result<Vec<f64>, JsError> {
    db()?.f1_chantler(element, energies).map_err(to_js)
}

/// f2 — imaginary part of the anomalous scattering factor (Chantler).
#[wasm_bindgen]
pub fn f2_chantler(element: &str, energies: &[f64]) -> Result<Vec<f64>, JsError> {
    db()?.f2_chantler(element, energies).map_err(to_js)
}

/// Chantler mass attenuation coefficient in cm²/g.
///
/// `kind` is one of `total`, `photo`, `incoherent`.
#[wasm_bindgen]
pub fn mu_chantler(element: &str, energies: &[f64], kind: &str) -> Result<Vec<f64>, JsError> {
    let k = parse_chantler_kind(kind)?;
    db()?.mu_chantler(element, energies, k).map_err(to_js)
}

/// Valid Chantler energy range for an element, as `[min, max]` in eV.
#[wasm_bindgen]
pub fn chantler_energy_range(element: &str) -> Result<Vec<f64>, JsError> {
    let (min, max) = db()?.chantler_energy_range(element).map_err(to_js)?;
    Ok(vec![min, max])
}

/// Elam energy range, as `[min, max]` in eV.
#[wasm_bindgen]
pub fn elam_energy_range() -> Vec<f64> {
    let (min, max) = XrayDb::elam_energy_range();
    vec![min, max]
}

/// f0 elastic scattering factor at the given q values in Å⁻¹.
#[wasm_bindgen]
pub fn f0(ion: &str, q: &[f64]) -> Result<Vec<f64>, JsError> {
    db()?.f0(ion, q).map_err(to_js)
}

/// Ion names with tabulated f0 coefficients, optionally filtered to one element.
#[wasm_bindgen]
pub fn f0_ions(element: Option<String>) -> Result<Vec<String>, JsError> {
    let db = db()?;
    let ions = db.f0_ions(element.as_deref()).map_err(to_js)?;
    Ok(ions.into_iter().map(str::to_string).collect())
}

// ── Edges and lines ───────────────────────────────────────────────────────

/// All absorption edges for an element, sorted by descending energy.
#[wasm_bindgen]
pub fn xray_edges(element: &str) -> Result<Vec<Edge>, JsError> {
    let db = db()?;
    let mut edges: Vec<Edge> = db
        .xray_edges(element)
        .map_err(to_js)?
        .into_iter()
        .map(|(label, e)| Edge {
            label,
            energy: e.energy,
            fluorescence_yield: e.fluorescence_yield,
            jump_ratio: e.jump_ratio,
        })
        .collect();
    edges.sort_by(|a, b| b.energy.total_cmp(&a.energy));
    Ok(edges)
}

/// Energy in eV of one absorption edge.
#[wasm_bindgen]
pub fn xray_edge_energy(element: &str, edge: &str) -> Result<f64, JsError> {
    db()?
        .xray_edge(element, edge)
        .map(|e| e.energy)
        .map_err(to_js)
}

/// Fluorescence yield for one edge.
#[wasm_bindgen]
pub fn fluorescence_yield(element: &str, edge: &str) -> Result<f64, JsError> {
    db()?
        .xray_edge(element, edge)
        .map(|e| e.fluorescence_yield)
        .map_err(to_js)
}

/// Absorption jump ratio for one edge.
#[wasm_bindgen]
pub fn jump_ratio(element: &str, edge: &str) -> Result<f64, JsError> {
    db()?
        .xray_edge(element, edge)
        .map(|e| e.jump_ratio)
        .map_err(to_js)
}

/// Emission lines for an element, sorted by descending intensity.
///
/// Pass `level` to restrict to one initial level, and `excitation_energy` to drop
/// lines that cannot be excited at that energy.
#[wasm_bindgen]
pub fn xray_lines(
    element: &str,
    level: Option<String>,
    excitation_energy: Option<f64>,
) -> Result<Vec<EmissionLine>, JsError> {
    let db = db()?;
    let lines = db
        .xray_lines_by_intensity(element, level.as_deref(), excitation_energy)
        .map_err(to_js)?;
    Ok(lines
        .into_iter()
        .map(|(label, l)| EmissionLine {
            label,
            energy: l.energy,
            intensity: l.intensity,
            initial_level: l.initial_level,
            final_level: l.final_level,
        })
        .collect())
}

/// Emission lines grouped by initial level (K, L1, L2, ...), each group sorted by
/// descending intensity.
///
/// Elam intensities are relative *within* a level, so this grouping is the correct
/// way to rank them; a global sort would put weak L lines above Ka1.
#[wasm_bindgen]
pub fn xray_lines_by_level(
    element: &str,
    excitation_energy: Option<f64>,
) -> Result<Vec<LineGroup>, JsError> {
    let db = db()?;
    let groups = db
        .xray_lines_by_level(element, excitation_energy)
        .map_err(to_js)?;
    Ok(groups
        .into_iter()
        .map(|(level, lines)| LineGroup {
            level,
            lines: lines
                .into_iter()
                .map(|(label, l)| EmissionLine {
                    label,
                    energy: l.energy,
                    intensity: l.intensity,
                    initial_level: l.initial_level,
                    final_level: l.final_level,
                })
                .collect(),
        })
        .collect())
}

/// Nearest absorption edge to an energy, or `undefined` if none is within tolerance.
#[wasm_bindgen]
pub fn guess_edge(energy: f64, max_difference: Option<f64>) -> Result<Option<EdgeMatch>, JsError> {
    Ok(db()?
        .guess_edge(energy, None, max_difference)
        .map(|g| EdgeMatch {
            element: g.element,
            edge: g.edge,
            energy: g.energy,
            difference: g.difference,
        }))
}

/// Core-hole width in eV.
#[wasm_bindgen]
pub fn corehole_width(element: &str, edge: &str) -> Result<f64, JsError> {
    db()?.core_width_at(element, edge).map_err(to_js)
}

// ── Materials ─────────────────────────────────────────────────────────────

/// Linear attenuation coefficient in 1/cm for a compound.
#[wasm_bindgen]
pub fn material_mu(
    formula: &str,
    density: f64,
    energies: &[f64],
    kind: &str,
) -> Result<Vec<f64>, JsError> {
    let k = parse_kind(kind)?;
    db()?
        .material_mu(formula, density, energies, k)
        .map_err(to_js)
}

/// Linear attenuation coefficient in 1/cm from elemental mass fractions.
#[wasm_bindgen]
pub fn material_mu_from_mass_fractions(
    symbols: Vec<String>,
    fractions: &[f64],
    density: f64,
    energies: &[f64],
    kind: &str,
) -> Result<Vec<f64>, JsError> {
    if symbols.len() != fractions.len() {
        return Err(JsError::new(&format!(
            "symbols/fractions length mismatch: {} != {}",
            symbols.len(),
            fractions.len()
        )));
    }
    let k = parse_kind(kind)?;
    let composition: Vec<(&str, f64)> = symbols
        .iter()
        .map(String::as_str)
        .zip(fractions.iter().copied())
        .collect();
    db()?
        .material_mu_from_mass_fractions(&composition, density, energies, k)
        .map_err(to_js)
}

/// Refractive index decrements for a compound.
#[wasm_bindgen]
pub fn xray_delta_beta(formula: &str, density: f64, energy: f64) -> Result<DeltaBeta, JsError> {
    let n = db()?
        .xray_delta_beta(formula, density, energy)
        .map_err(to_js)?;
    Ok(DeltaBeta {
        delta: n.delta,
        beta: n.beta,
        attenuation_length_cm: n.attenuation_length_cm,
    })
}

/// Look up a material by name or formula.
#[wasm_bindgen]
pub fn find_material(name: &str) -> Result<Option<MaterialInfo>, JsError> {
    Ok(db()?.find_material(name).map(|m| MaterialInfo {
        name: m.name.to_string(),
        formula: m.formula.to_string(),
        density: m.density,
    }))
}

/// Every material in the built-in table.
#[wasm_bindgen]
pub fn materials() -> Result<Vec<MaterialInfo>, JsError> {
    Ok(db()?
        .materials()
        .into_iter()
        .map(|m| MaterialInfo {
            name: m.name.to_string(),
            formula: m.formula.to_string(),
            density: m.density,
        })
        .collect())
}

/// True if a chemical formula parses.
#[wasm_bindgen]
pub fn validate_formula(formula: &str) -> Result<bool, JsError> {
    Ok(db()?.validate_formula(formula))
}

// ── Misc ──────────────────────────────────────────────────────────────────

/// Compton scattering energies for an incident energy in eV.
#[wasm_bindgen]
pub fn compton_energies(incident_energy: f64) -> Result<Compton, JsError> {
    let c = db()?.compton_energies(incident_energy);
    Ok(Compton {
        xray_90deg: c.xray_90deg,
        xray_mean: c.xray_mean,
        electron_mean: c.electron_mean,
    })
}

/// Ionization potential in eV per ion pair for a counting gas.
#[wasm_bindgen]
pub fn ionization_potential(gas: &str) -> Result<f64, JsError> {
    db()?.ionization_potential(gas).map_err(to_js)
}
