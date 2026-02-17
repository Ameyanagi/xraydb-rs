#![no_std]

extern crate alloc;

use alloc::string::String;
use alloc::vec::Vec;
use serde::{Deserialize, Serialize};

/// The complete XrayDB database, deserialized from the embedded blob.
#[derive(Debug, Serialize, Deserialize)]
pub struct XrayDatabase {
    pub version: Vec<VersionRecord>,
    pub elements: Vec<ElementRecord>,
    pub xray_levels: Vec<XrayLevelRecord>,
    pub xray_transitions: Vec<XrayTransitionRecord>,
    pub coster_kronig: Vec<CosterKronigRecord>,
    pub photoabsorption: Vec<PhotoabsorptionRecord>,
    pub scattering: Vec<ScatteringRecord>,
    pub chantler: Vec<ChantlerRecord>,
    pub waasmaier: Vec<WaasmaierRecord>,
    pub compton_energies: ComptonEnergiesRecord,
    pub keski_rahkonen_krause: Vec<CoreWidthRecord>,
    pub krause_oliver: Vec<CoreWidthRecord>,
    pub corelevel_widths: Vec<CoreWidthRecord>,
    pub ionization_potentials: Vec<IonizationPotentialRecord>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct VersionRecord {
    pub tag: String,
    pub date: String,
    pub notes: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ElementRecord {
    pub atomic_number: u16,
    pub symbol: String,
    pub name: String,
    pub molar_mass: f64,
    pub density: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct XrayLevelRecord {
    pub element: String,
    pub iupac_symbol: String,
    pub absorption_edge: f64,
    pub fluorescence_yield: f64,
    pub jump_ratio: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct XrayTransitionRecord {
    pub element: String,
    pub iupac_symbol: String,
    pub siegbahn_symbol: String,
    pub initial_level: String,
    pub final_level: String,
    pub emission_energy: f64,
    pub intensity: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CosterKronigRecord {
    pub element: String,
    pub initial_level: String,
    pub final_level: String,
    pub transition_probability: f64,
    pub total_transition_probability: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhotoabsorptionRecord {
    pub element: String,
    pub log_energy: Vec<f64>,
    pub log_photoabsorption: Vec<f64>,
    pub log_photoabsorption_spline: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ScatteringRecord {
    pub element: String,
    pub log_energy: Vec<f64>,
    pub log_coherent_scatter: Vec<f64>,
    pub log_coherent_scatter_spline: Vec<f64>,
    pub log_incoherent_scatter: Vec<f64>,
    pub log_incoherent_scatter_spline: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ChantlerRecord {
    pub element: String,
    pub sigma_mu: f64,
    pub mue_f2: f64,
    pub density: f64,
    pub corr_henke: f64,
    pub corr_cl35: f64,
    pub corr_nucl: f64,
    pub energy: Vec<f64>,
    pub f1: Vec<f64>,
    pub f2: Vec<f64>,
    pub mu_photo: Vec<f64>,
    pub mu_incoh: Vec<f64>,
    pub mu_total: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct WaasmaierRecord {
    pub atomic_number: u16,
    pub element: String,
    pub ion: String,
    pub offset: f64,
    pub scale: Vec<f64>,
    pub exponents: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ComptonEnergiesRecord {
    pub incident: Vec<f64>,
    pub xray_90deg: Vec<f64>,
    pub xray_mean: Vec<f64>,
    pub electron_mean: Vec<f64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CoreWidthRecord {
    pub atomic_number: u16,
    pub element: String,
    pub edge: String,
    pub width: f64,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct IonizationPotentialRecord {
    pub gas: String,
    pub potential: f64,
}
