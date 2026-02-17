use crate::db::XrayDb;
use crate::elam::CrossSectionKind;
use crate::error::{Result, XrayDbError};
use crate::materials_db::find_material;

/// Ion chamber flux results.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IonChamberFluxes {
    pub incident: f64,
    pub transmitted: f64,
    pub photo: f64,
    pub incoherent: f64,
    pub coherent: f64,
}

impl XrayDb {
    /// Lookup a material name, returning (formula, density).
    ///
    /// Looks up the embedded materials database by name (case-insensitive)
    /// or by chemical formula. Returns `None` if not found.
    pub fn find_material(&self, name: &str) -> Option<(&'static str, f64)> {
        find_material(name)
    }

    /// Material mu by name: looks up formula and density from the materials database.
    ///
    /// If the name is not found in the database, it is treated as a formula
    /// and `density` must be provided.
    pub fn material_mu_named(
        &self,
        name: &str,
        energies: &[f64],
        kind: CrossSectionKind,
        density: Option<f64>,
    ) -> Result<Vec<f64>> {
        let (formula, dens) = if let Some((f, d)) = find_material(name) {
            (f, density.unwrap_or(d))
        } else {
            let d = density.ok_or_else(|| {
                XrayDbError::DataError(format!(
                    "unknown material '{name}', density must be provided"
                ))
            })?;
            (name, d)
        };
        self.material_mu(formula, dens, energies, kind)
    }

    /// Calculate ion chamber fluxes from measured voltage.
    ///
    /// # Arguments
    /// * `gases` - Gas mixture as (name, fraction) pairs. Use `&[("nitrogen", 1.0)]` for pure N₂.
    /// * `volts` - Measured voltage
    /// * `length_cm` - Active length of ion chamber in cm
    /// * `energy` - X-ray energy in eV
    /// * `sensitivity` - Current sensitivity in A/V
    /// * `with_compton` - Include Compton electron energy contribution
    /// * `both_carriers` - Count both electron and ion carriers (true for most chambers)
    #[allow(clippy::too_many_arguments)]
    pub fn ionchamber_fluxes(
        &self,
        gases: &[(&str, f64)],
        volts: f64,
        length_cm: f64,
        energy: f64,
        sensitivity: f64,
        with_compton: bool,
        both_carriers: bool,
    ) -> Result<IonChamberFluxes> {
        let ncarriers: f64 = if both_carriers { 2.0 } else { 1.0 };

        // Normalize fractions
        let gas_total: f64 = gases.iter().map(|(_, f)| f).sum();
        if gas_total <= 0.0 {
            return Err(XrayDbError::DataError(
                "gas fractions must sum to > 0".to_string(),
            ));
        }

        // Compton electron mean energy
        let energy_compton = if with_compton {
            self.compton_energies(energy).electron_mean
        } else {
            0.0
        };

        // Weighted sums of mu values and ionization potential
        let mut mu_photo = 0.0;
        let mut mu_incoh = 0.0;
        let mut mu_total = 0.0;
        let mut mu_coh = 0.0;
        let mut ion_pot = 0.0;

        let e_arr = [energy];

        for &(gas_name, frac) in gases {
            let weight = frac / gas_total;

            // Resolve gas name: "N2" -> "nitrogen", "O2" -> "oxygen"
            let lookup_name = match gas_name {
                "N2" => "nitrogen",
                "O2" => "oxygen",
                other => other,
            };

            // Get ionization potential
            let ip = self
                .ionization_potential(gas_name)
                .or_else(|_| self.ionization_potential(lookup_name))
                .unwrap_or(32.0); // default fallback

            // Compute material_mu for each kind
            let photo =
                self.material_mu_named(lookup_name, &e_arr, CrossSectionKind::Photo, None)?[0];
            let total =
                self.material_mu_named(lookup_name, &e_arr, CrossSectionKind::Total, None)?[0];
            let incoh =
                self.material_mu_named(lookup_name, &e_arr, CrossSectionKind::Incoherent, None)?[0];
            let coh =
                self.material_mu_named(lookup_name, &e_arr, CrossSectionKind::Coherent, None)?[0];

            mu_photo += photo * weight;
            mu_total += total * weight;
            mu_incoh += incoh * weight;
            mu_coh += coh * weight;
            ion_pot += ip * weight;
        }

        let atten_total = 1.0 - (-length_cm * mu_total).exp();
        let atten_photo = if mu_total > 0.0 {
            atten_total * mu_photo / mu_total
        } else {
            0.0
        };
        let atten_incoh = if mu_total > 0.0 {
            atten_total * mu_incoh / mu_total
        } else {
            0.0
        };
        let atten_coh = if mu_total > 0.0 {
            atten_total * mu_coh / mu_total
        } else {
            0.0
        };

        let absorbed_energy = ncarriers * (energy * atten_photo + energy_compton * atten_incoh);

        let qcharge = crate::constants::ELEMENTARY_CHARGE;
        let flux_in = if absorbed_energy > 0.0 {
            volts * sensitivity * ion_pot / (qcharge * absorbed_energy)
        } else {
            0.0
        };

        Ok(IonChamberFluxes {
            incident: flux_in,
            transmitted: flux_in * (1.0 - atten_total),
            photo: flux_in * atten_photo,
            incoherent: flux_in * atten_incoh,
            coherent: flux_in * atten_coh,
        })
    }
}
