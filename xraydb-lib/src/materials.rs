use std::collections::HashMap;

use crate::chemparser::chemparse;
use crate::constants::{AVOGADRO, R_ELECTRON_CM};
use crate::db::XrayDb;
use crate::elam::CrossSectionKind;
use crate::error::{Result, XrayDbError};

impl XrayDb {
    /// Returns mass attenuation coefficient for a material in 1/cm.
    ///
    /// # Arguments
    /// * `formula` - Chemical formula (e.g., "H2O", "SiO2")
    /// * `density` - Density in g/cm³
    /// * `energies` - X-ray energies in eV
    /// * `kind` - Cross-section kind (total, photo, coherent, incoherent)
    pub fn material_mu(
        &self,
        formula: &str,
        density: f64,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        let composition = chemparse(formula)?;

        // Calculate total formula weight
        let total_weight: f64 = composition
            .iter()
            .map(|(sym, &count)| count * self.molar_mass(sym).unwrap_or(0.0))
            .sum();

        if total_weight <= 0.0 {
            return Err(XrayDbError::InvalidFormula(format!(
                "zero weight formula: {formula}"
            )));
        }

        // Calculate mass fractions
        let mass_fracs: HashMap<&str, f64> = composition
            .iter()
            .map(|(sym, &count)| {
                let frac = count * self.molar_mass(sym).unwrap_or(0.0) / total_weight;
                (sym.as_str(), frac)
            })
            .collect();

        // Sum weighted cross-sections
        let mut mu = vec![0.0_f64; energies.len()];
        for (sym, &frac) in &mass_fracs {
            let elem_mu = self.mu_elam(sym, energies, kind)?;
            for (i, &val) in elem_mu.iter().enumerate() {
                mu[i] += frac * val;
            }
        }

        // Convert from cm²/g to 1/cm by multiplying by density
        for val in &mut mu {
            *val *= density;
        }

        Ok(mu)
    }

    /// Returns X-ray refractive index components (delta, beta, attenuation_length_cm).
    ///
    /// The complex index of refraction is: n = 1 - delta - i*beta
    ///
    /// # Arguments
    /// * `formula` - Chemical formula
    /// * `density` - Density in g/cm³
    /// * `energy` - X-ray energy in eV
    pub fn xray_delta_beta(
        &self,
        formula: &str,
        density: f64,
        energy: f64,
    ) -> Result<(f64, f64, f64)> {
        let composition = chemparse(formula)?;

        // wavelength in cm
        let wavelength = 1.0e-7 * crate::constants::PLANCK_HC / energy;

        let total_weight: f64 = composition
            .iter()
            .map(|(sym, &count)| count * self.molar_mass(sym).unwrap_or(0.0))
            .sum();

        if total_weight <= 0.0 {
            return Err(XrayDbError::InvalidFormula(format!(
                "zero weight formula: {formula}"
            )));
        }

        // Sum f1 and f2 contributions weighted by composition
        let mut sum_f1 = 0.0_f64;
        let mut sum_f2 = 0.0_f64;

        for (sym, &count) in &composition {
            let z = self.atomic_number(sym)? as f64;
            let f1_vals = self.f1_chantler(sym, &[energy])?;
            let f2_vals = self.f2_chantler(sym, &[energy])?;
            // f1 from Chantler is f' (anomalous correction), so full f1 = Z + f'
            sum_f1 += count * (z + f1_vals[0]);
            sum_f2 += count * f2_vals[0];
        }

        // Number density factor
        let prefactor =
            R_ELECTRON_CM * wavelength * wavelength * density * AVOGADRO / (2.0 * std::f64::consts::PI * total_weight);

        let delta = prefactor * sum_f1;
        let beta = prefactor * sum_f2;
        let atlen = if beta > 0.0 {
            wavelength / (4.0 * std::f64::consts::PI * beta)
        } else {
            f64::INFINITY
        };

        Ok((delta, beta, atlen))
    }
}
