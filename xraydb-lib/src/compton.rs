//! Compton scattering energies.

use crate::db::XrayDb;
use crate::interp::interp_one;

/// Compton scattering energies for one incident photon energy, all in eV.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComptonEnergies {
    /// The incident photon energy that was queried.
    pub incident: f64,
    /// Energy of a photon Compton-scattered through 90°.
    pub xray_90deg: f64,
    /// Mean energy of the Compton-scattered photon over all angles.
    pub xray_mean: f64,
    /// Mean energy imparted to the recoil electron.
    pub electron_mean: f64,
}

impl XrayDb {
    /// Compton scattering energies for a given incident X-ray energy in eV.
    ///
    /// Values are interpolated from the tabulated grid and clamped at its endpoints.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let c = db.compton_energies(20_000.0);
    /// assert!(c.xray_90deg < c.incident);
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
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
