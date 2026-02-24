use std::collections::HashMap;

use chemical_formula::prelude::parse_formula;

use crate::chemparser::chemparse;
use crate::constants::{AVOGADRO, R_ELECTRON_CM};
use crate::db::XrayDb;
use crate::elam::CrossSectionKind;
use crate::error::{Result, XrayDbError};

fn validate_material_density(density: f64) -> Result<()> {
    if density.is_finite() && density > 0.0 {
        return Ok(());
    }

    Err(XrayDbError::DataError(format!(
        "density must be finite and > 0 g/cm^3, got {density}"
    )))
}

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
        validate_material_density(density)?;

        let molecular_formula = parse_formula(formula)
            .and_then(|parsed| parsed.to_molecular_formula())
            .map_err(|err| XrayDbError::InvalidFormula(format!("{formula}: {err}")))?;

        if molecular_formula.stoichiometry.is_empty() {
            return Err(XrayDbError::InvalidFormula(format!(
                "zero weight formula: {formula}"
            )));
        }

        let mut mass_fractions =
            HashMap::<String, f64>::with_capacity(molecular_formula.stoichiometry.len());
        let mut total_mass = 0.0_f64;
        for (element, count) in molecular_formula.stoichiometry {
            if !count.is_finite() || count < 0.0 {
                return Err(XrayDbError::InvalidFormula(format!(
                    "invalid stoichiometry for {element} in formula '{formula}': {count}"
                )));
            }
            if count == 0.0 {
                continue;
            }

            let symbol = element.to_string();
            if symbol.is_empty() {
                return Err(XrayDbError::InvalidFormula(format!(
                    "invalid element symbol in formula: {formula}"
                )));
            }

            let mass = count * self.molar_mass(&symbol)?;
            if !mass.is_finite() || mass < 0.0 {
                return Err(XrayDbError::InvalidFormula(format!(
                    "invalid mass contribution for element {symbol} in formula '{formula}'"
                )));
            }

            *mass_fractions.entry(symbol).or_insert(0.0) += mass;
            total_mass += mass;
        }

        if !total_mass.is_finite() || total_mass <= 0.0 {
            return Err(XrayDbError::InvalidFormula(format!(
                "zero weight formula: {formula}"
            )));
        }

        let normalized: Vec<(String, f64)> = mass_fractions
            .into_iter()
            .map(|(symbol, mass)| (symbol, mass / total_mass))
            .collect();

        self.material_mu_from_normalized_mass_fractions(&normalized, density, energies, kind)
    }

    /// Returns material linear attenuation coefficient from elemental mass fractions.
    ///
    /// # Arguments
    /// * `composition` - Elemental mass fractions as `(symbol, fraction)` pairs.
    /// * `density` - Density in g/cm³
    /// * `energies` - X-ray energies in eV
    /// * `kind` - Cross-section kind (total, photo, coherent, incoherent)
    pub fn material_mu_from_mass_fractions(
        &self,
        composition: &[(&str, f64)],
        density: f64,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        validate_material_density(density)?;

        if composition.is_empty() {
            return Err(XrayDbError::DataError(
                "composition mass fractions must not be empty".to_string(),
            ));
        }

        let mut by_symbol = HashMap::<String, f64>::with_capacity(composition.len());
        let mut total_fraction = 0.0_f64;
        for &(symbol, fraction) in composition {
            if !fraction.is_finite() || fraction < 0.0 {
                return Err(XrayDbError::DataError(format!(
                    "invalid mass fraction for '{symbol}': {fraction}"
                )));
            }
            if fraction == 0.0 {
                continue;
            }

            let canonical = self.symbol(symbol)?.to_string();
            *by_symbol.entry(canonical).or_insert(0.0) += fraction;
            total_fraction += fraction;
        }

        if !total_fraction.is_finite() || total_fraction <= 0.0 {
            return Err(XrayDbError::DataError(
                "composition mass fractions must include at least one positive value".to_string(),
            ));
        }

        let normalized: Vec<(String, f64)> = by_symbol
            .into_iter()
            .map(|(symbol, fraction)| (symbol, fraction / total_fraction))
            .collect();

        self.material_mu_from_normalized_mass_fractions(&normalized, density, energies, kind)
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
        let mut species = Vec::with_capacity(composition.len());
        for (sym, &count) in &composition {
            let molar_mass = self.molar_mass(sym)?;
            species.push((sym.as_str(), count, molar_mass));
        }

        // wavelength in cm
        let wavelength = 1.0e-7 * crate::constants::PLANCK_HC / energy;

        let total_weight: f64 = species
            .iter()
            .map(|(_, count, molar_mass)| count * molar_mass)
            .sum();

        if total_weight <= 0.0 {
            return Err(XrayDbError::InvalidFormula(format!(
                "zero weight formula: {formula}"
            )));
        }

        // Sum f1 and f2 contributions weighted by composition
        let mut sum_f1 = 0.0_f64;
        let mut sum_f2 = 0.0_f64;

        for &(sym, count, _) in &species {
            let z = self.atomic_number(sym)? as f64;
            let f1_vals = self.f1_chantler(sym, &[energy])?;
            let f2_vals = self.f2_chantler(sym, &[energy])?;
            // f1 from Chantler is f' (anomalous correction), so full f1 = Z + f'
            sum_f1 += count * (z + f1_vals[0]);
            sum_f2 += count * f2_vals[0];
        }

        // Number density factor
        let prefactor = R_ELECTRON_CM * wavelength * wavelength * density * AVOGADRO
            / (2.0 * std::f64::consts::PI * total_weight);

        let delta = prefactor * sum_f1;
        let beta = prefactor * sum_f2;
        let atlen = if beta > 0.0 {
            wavelength / (4.0 * std::f64::consts::PI * beta)
        } else {
            f64::INFINITY
        };

        Ok((delta, beta, atlen))
    }

    fn material_mu_from_normalized_mass_fractions(
        &self,
        normalized_mass_fractions: &[(String, f64)],
        density: f64,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        let mut mu = vec![0.0_f64; energies.len()];
        let mut total_fraction = 0.0_f64;

        for (symbol, fraction) in normalized_mass_fractions {
            if !fraction.is_finite() || *fraction < 0.0 {
                return Err(XrayDbError::DataError(format!(
                    "invalid normalized mass fraction for '{symbol}': {fraction}"
                )));
            }
            if *fraction == 0.0 {
                continue;
            }

            total_fraction += *fraction;
            let elem_mu = self.mu_elam(symbol, energies, kind)?;
            for (acc, &value) in mu.iter_mut().zip(elem_mu.iter()) {
                *acc += fraction * value;
            }
        }

        if !total_fraction.is_finite() || total_fraction <= 0.0 {
            return Err(XrayDbError::DataError(
                "composition mass fractions must include at least one positive value".to_string(),
            ));
        }

        for value in &mut mu {
            *value *= density;
        }

        Ok(mu)
    }
}
