//! WASM bindings for XrayDB.
//!
//! Build with:
//! ```sh
//! wasm-pack build -p xraydb-wasm
//! ```

use wasm_bindgen::prelude::*;

use xraydb::{ChantlerKind, CrossSectionKind, XrayDb};

fn db() -> XrayDb {
    XrayDb::new()
}

fn parse_kind(kind: &str) -> Result<CrossSectionKind, JsError> {
    match kind.to_lowercase().as_str() {
        "total" => Ok(CrossSectionKind::Total),
        "photo" => Ok(CrossSectionKind::Photo),
        "coherent" | "coh" => Ok(CrossSectionKind::Coherent),
        "incoherent" | "incoh" => Ok(CrossSectionKind::Incoherent),
        _ => Err(JsError::new(&format!("unknown cross-section kind: {kind}"))),
    }
}

fn to_js(e: xraydb::XrayDbError) -> JsError {
    JsError::new(&e.to_string())
}

// ── Element lookups ──

#[wasm_bindgen]
pub fn atomic_number(element: &str) -> Result<u16, JsError> {
    db().atomic_number(element).map_err(to_js)
}

#[wasm_bindgen]
pub fn symbol(element: &str) -> Result<String, JsError> {
    db().symbol(element).map(|s| s.to_string()).map_err(to_js)
}

#[wasm_bindgen]
pub fn atomic_name(element: &str) -> Result<String, JsError> {
    db().atomic_name(element).map(|s| s.to_string()).map_err(to_js)
}

#[wasm_bindgen]
pub fn molar_mass(element: &str) -> Result<f64, JsError> {
    db().molar_mass(element).map_err(to_js)
}

#[wasm_bindgen]
pub fn element_density(element: &str) -> Result<f64, JsError> {
    db().density(element).map_err(to_js)
}

// ── Elam cross-sections ──

/// Returns mass attenuation coefficients (cm²/g) from Elam tables.
///
/// `kind` is one of: "total", "photo", "coherent", "incoherent".
#[wasm_bindgen]
pub fn mu_elam(element: &str, energies: &[f64], kind: &str) -> Result<Vec<f64>, JsError> {
    let k = parse_kind(kind)?;
    db().mu_elam(element, energies, k).map_err(to_js)
}

// ── Chantler data ──

/// Returns f1 (anomalous scattering factor, real part) from Chantler tables.
#[wasm_bindgen]
pub fn f1_chantler(element: &str, energies: &[f64]) -> Result<Vec<f64>, JsError> {
    db().f1_chantler(element, energies).map_err(to_js)
}

/// Returns f2 (anomalous scattering factor, imaginary part) from Chantler tables.
#[wasm_bindgen]
pub fn f2_chantler(element: &str, energies: &[f64]) -> Result<Vec<f64>, JsError> {
    db().f2_chantler(element, energies).map_err(to_js)
}

/// Returns Chantler mass attenuation coefficient (cm²/g).
///
/// `kind` is one of: "total", "photo", "incoherent".
#[wasm_bindgen]
pub fn mu_chantler(element: &str, energies: &[f64], kind: &str) -> Result<Vec<f64>, JsError> {
    let k = match kind.to_lowercase().as_str() {
        "total" => ChantlerKind::Total,
        "photo" => ChantlerKind::Photo,
        "incoherent" | "incoh" => ChantlerKind::Incoherent,
        _ => return Err(JsError::new(&format!("unknown Chantler kind: {kind}"))),
    };
    db().mu_chantler(element, energies, k).map_err(to_js)
}

// ── Waasmaier-Kirfel f0 ──

/// Returns f0 elastic scattering factor at given q values (Å⁻¹).
#[wasm_bindgen]
pub fn f0(ion: &str, q: &[f64]) -> Result<Vec<f64>, JsError> {
    db().f0(ion, q).map_err(to_js)
}

// ── X-ray edges and lines ──

/// Returns X-ray edge energy (eV) for an element and edge label.
#[wasm_bindgen]
pub fn xray_edge_energy(element: &str, edge: &str) -> Result<f64, JsError> {
    db().xray_edge(element, edge)
        .map(|e| e.energy)
        .map_err(to_js)
}

/// Returns fluorescence yield for an element and edge label.
#[wasm_bindgen]
pub fn fluorescence_yield(element: &str, edge: &str) -> Result<f64, JsError> {
    db().xray_edge(element, edge)
        .map(|e| e.fluorescence_yield)
        .map_err(to_js)
}

/// Returns jump ratio for an element and edge label.
#[wasm_bindgen]
pub fn jump_ratio(element: &str, edge: &str) -> Result<f64, JsError> {
    db().xray_edge(element, edge)
        .map(|e| e.jump_ratio)
        .map_err(to_js)
}

// ── Materials ──

/// Returns material linear attenuation coefficient (1/cm).
///
/// `kind` is one of: "total", "photo", "coherent", "incoherent".
#[wasm_bindgen]
pub fn material_mu(
    formula: &str,
    density: f64,
    energies: &[f64],
    kind: &str,
) -> Result<Vec<f64>, JsError> {
    let k = parse_kind(kind)?;
    db().material_mu(formula, density, energies, k)
        .map_err(to_js)
}

/// Returns [delta, beta, attenuation_length_cm] for a material.
///
/// The complex refractive index is n = 1 - delta - i*beta.
#[wasm_bindgen]
pub fn xray_delta_beta(formula: &str, density: f64, energy: f64) -> Result<Vec<f64>, JsError> {
    let (delta, beta, atlen) = db()
        .xray_delta_beta(formula, density, energy)
        .map_err(to_js)?;
    Ok(vec![delta, beta, atlen])
}

// ── Compton energies ──

/// Returns [xray_90deg, xray_mean, electron_mean] for a given incident energy.
#[wasm_bindgen]
pub fn compton_energies(incident_energy: f64) -> Vec<f64> {
    let c = db().compton_energies(incident_energy);
    vec![c.xray_90deg, c.xray_mean, c.electron_mean]
}

// ── Core widths ──

/// Returns core-hole width (eV) for an element and edge.
#[wasm_bindgen]
pub fn corehole_width(element: &str, edge: &str) -> Result<f64, JsError> {
    let widths = db().core_width(element, Some(edge)).map_err(to_js)?;
    widths
        .get(edge)
        .copied()
        .ok_or_else(|| JsError::new(&format!("no width for edge '{edge}'")))
}

/// Returns ionization potential (eV per ion pair) for a gas.
#[wasm_bindgen]
pub fn ionization_potential(gas: &str) -> Result<f64, JsError> {
    db().ionization_potential(gas).map_err(to_js)
}
