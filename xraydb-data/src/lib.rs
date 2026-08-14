//! Serde data model shared by the `xraydb` crates.
//!
//! This crate is `#![no_std]` and holds only plain data structures; all lookup and
//! calculation lives in `xraydb`. It exists so the generator binary and the library
//! agree on the on-disk format without the library depending on the generator.

#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(clippy::unwrap_used, clippy::expect_used)]

extern crate alloc;

use alloc::string::String;
use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

/// The complete XrayDB database, deserialized from the embedded blob.
#[derive(Debug, Serialize, Deserialize)]
pub struct XrayDatabase {
    /// Upstream XrayDB version records.
    pub version: Vec<VersionRecord>,
    /// One record per element, ordered by atomic number.
    pub elements: Vec<ElementRecord>,
    /// Absorption edges (one row per element and IUPAC level).
    pub xray_levels: Vec<XrayLevelRecord>,
    /// Emission lines (one row per element and Siegbahn label).
    pub xray_transitions: Vec<XrayTransitionRecord>,
    /// Coster-Kronig transition probabilities.
    pub coster_kronig: Vec<CosterKronigRecord>,
    /// Elam photoabsorption cross-section splines.
    pub photoabsorption: Vec<PhotoabsorptionRecord>,
    /// Elam coherent and incoherent scattering splines.
    pub scattering: Vec<ScatteringRecord>,
    /// Chantler anomalous-scattering and attenuation tables.
    pub chantler: Vec<ChantlerRecord>,
    /// Waasmaier-Kirfel f0 coefficients, one row per ion.
    pub waasmaier: Vec<WaasmaierRecord>,
    /// Tabulated Compton scattering energies.
    pub compton_energies: ComptonEnergiesRecord,
    /// Core widths from Keski-Rahkonen & Krause.
    pub keski_rahkonen_krause: Vec<CoreWidthRecord>,
    /// Core widths from Krause & Oliver (1979).
    pub krause_oliver: Vec<CoreWidthRecord>,
    /// Merged core widths: KK updated with KO for K, L1, L2, L3.
    pub corelevel_widths: Vec<CoreWidthRecord>,
    /// Ion-chamber gas ionization potentials, eV per ion pair.
    pub ionization_potentials: Vec<IonizationPotentialRecord>,
}

/// One upstream XrayDB version entry.
#[derive(Debug, Serialize, Deserialize)]
pub struct VersionRecord {
    /// Version tag, e.g. `"1.0"`.
    pub tag: String,
    /// Release date as written upstream.
    pub date: String,
    /// Free-text release notes.
    pub notes: String,
}

/// Basic facts about one element.
#[derive(Debug, Serialize, Deserialize)]
pub struct ElementRecord {
    /// Atomic number Z.
    pub atomic_number: u16,
    /// Element symbol, e.g. `"Fe"`.
    pub symbol: String,
    /// Element name, e.g. `"iron"`.
    pub name: String,
    /// Molar mass in g/mol.
    pub molar_mass: f64,
    /// Density in g/cm³ at standard conditions.
    pub density: f64,
}

/// One absorption edge of one element.
#[derive(Debug, Serialize, Deserialize)]
pub struct XrayLevelRecord {
    /// Element symbol this row belongs to.
    pub element: String,
    /// IUPAC level label, e.g. `"K"`, `"L3"`.
    pub iupac_symbol: String,
    /// Absorption edge energy in eV.
    pub absorption_edge: f64,
    /// Fluorescence yield, 0-1.
    pub fluorescence_yield: f64,
    /// Absorption jump ratio across the edge.
    pub jump_ratio: f64,
}

/// One emission line of one element.
#[derive(Debug, Serialize, Deserialize)]
pub struct XrayTransitionRecord {
    /// Element symbol this row belongs to.
    pub element: String,
    /// IUPAC level label, e.g. `"K"`, `"L3"`.
    pub iupac_symbol: String,
    /// Siegbahn line label, e.g. `"Ka1"`.
    pub siegbahn_symbol: String,
    /// IUPAC label of the initial (core-hole) level.
    pub initial_level: String,
    /// IUPAC label of the final level.
    pub final_level: String,
    /// Emission energy in eV.
    pub emission_energy: f64,
    /// Relative intensity, normalized within the initial level.
    pub intensity: f64,
}

/// One Coster-Kronig transition between levels of the same shell.
#[derive(Debug, Serialize, Deserialize)]
pub struct CosterKronigRecord {
    /// Element symbol this row belongs to.
    pub element: String,
    /// IUPAC label of the initial (core-hole) level.
    pub initial_level: String,
    /// IUPAC label of the final level.
    pub final_level: String,
    /// Direct transition probability.
    pub transition_probability: f64,
    /// Total probability including intermediate states.
    pub total_transition_probability: f64,
}

/// Elam photoabsorption spline for one element, in log-log space.
#[derive(Debug, Serialize, Deserialize)]
pub struct PhotoabsorptionRecord {
    /// Element symbol this row belongs to.
    pub element: String,
    /// Natural log of the tabulated energies in eV.
    pub log_energy: Vec<f64>,
    /// Natural log of the photoabsorption cross-section in cm²/g.
    pub log_photoabsorption: Vec<f64>,
    /// Second derivatives for `log_photoabsorption`.
    pub log_photoabsorption_spline: Vec<f64>,
}

/// Elam coherent and incoherent scattering splines for one element.
#[derive(Debug, Serialize, Deserialize)]
pub struct ScatteringRecord {
    /// Element symbol this row belongs to.
    pub element: String,
    /// Natural log of the tabulated energies in eV.
    pub log_energy: Vec<f64>,
    /// Natural log of the coherent scattering cross-section in cm²/g.
    pub log_coherent_scatter: Vec<f64>,
    /// Second derivatives for `log_coherent_scatter`.
    pub log_coherent_scatter_spline: Vec<f64>,
    /// Natural log of the incoherent scattering cross-section in cm²/g.
    pub log_incoherent_scatter: Vec<f64>,
    /// Second derivatives for `log_incoherent_scatter`.
    pub log_incoherent_scatter_spline: Vec<f64>,
}

/// Chantler tables for one element.
#[derive(Debug, Serialize, Deserialize)]
pub struct ChantlerRecord {
    /// Element symbol this row belongs to.
    pub element: String,
    /// Chantler sigma_mu normalization constant.
    pub sigma_mu: f64,
    /// Chantler mu_e/f2 conversion constant.
    pub mue_f2: f64,
    /// Density in g/cm³ at standard conditions.
    pub density: f64,
    /// Henke-scale correction term.
    pub corr_henke: f64,
    /// Cromer-Liberman 3.5 keV correction term.
    pub corr_cl35: f64,
    /// Nuclear Thomson correction term.
    pub corr_nucl: f64,
    /// Tabulated energies in eV.
    pub energy: Vec<f64>,
    /// Anomalous scattering factor f' (Z already subtracted).
    pub f1: Vec<f64>,
    /// Anomalous scattering factor f".
    pub f2: Vec<f64>,
    /// Photoelectric mass attenuation in cm²/g.
    pub mu_photo: Vec<f64>,
    /// Incoherent mass attenuation in cm²/g.
    pub mu_incoh: Vec<f64>,
    /// Total mass attenuation in cm²/g.
    pub mu_total: Vec<f64>,
}

/// Waasmaier-Kirfel f0 coefficients for one ion.
#[derive(Debug, Serialize, Deserialize)]
pub struct WaasmaierRecord {
    /// Atomic number Z.
    pub atomic_number: u16,
    /// Element symbol this row belongs to.
    pub element: String,
    /// Ion label, e.g. `"Fe"`, `"Fe2+"`.
    pub ion: String,
    /// Constant term c in f0(q) = c + sum(a_i exp(-b_i q²)).
    pub offset: f64,
    /// Amplitudes a_i.
    pub scale: Vec<f64>,
    /// Exponents b_i in Å².
    pub exponents: Vec<f64>,
}

/// Tabulated Compton scattering energies.
#[derive(Debug, Serialize, Deserialize)]
pub struct ComptonEnergiesRecord {
    /// Incident photon energies in eV.
    pub incident: Vec<f64>,
    /// Energy of a photon scattered through 90°, in eV.
    pub xray_90deg: Vec<f64>,
    /// Mean scattered-photon energy in eV.
    pub xray_mean: Vec<f64>,
    /// Mean recoil-electron energy in eV.
    pub electron_mean: Vec<f64>,
}

/// Core-hole width for one element and edge.
#[derive(Debug, Serialize, Deserialize)]
pub struct CoreWidthRecord {
    /// Atomic number Z.
    pub atomic_number: u16,
    /// Element symbol this row belongs to.
    pub element: String,
    /// IUPAC edge label.
    pub edge: String,
    /// Core-hole width in eV.
    pub width: f64,
}

/// Ionization potential for one counting gas.
#[derive(Debug, Serialize, Deserialize)]
pub struct IonizationPotentialRecord {
    /// Gas name.
    pub gas: String,
    /// Ionization potential in eV per ion pair.
    pub potential: f64,
}
