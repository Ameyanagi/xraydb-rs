use crate::db::XrayDb;
use crate::interp::interp_one;

/// Compton scattering energies for a given incident energy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComptonEnergies {
    pub incident: f64,
    pub xray_90deg: f64,
    pub xray_mean: f64,
    pub electron_mean: f64,
}

impl XrayDb {
    /// Returns Compton scattering energies for a given incident X-ray energy (eV).
    pub fn compton_energies(&self, incident_energy: f64) -> ComptonEnergies {
        let data = &self.raw().compton_energies;

        ComptonEnergies {
            incident: incident_energy,
            xray_90deg: interp_one(incident_energy, &data.incident, &data.xray_90deg),
            xray_mean: interp_one(incident_energy, &data.incident, &data.xray_mean),
            electron_mean: interp_one(incident_energy, &data.incident, &data.electron_mean),
        }
    }
}
