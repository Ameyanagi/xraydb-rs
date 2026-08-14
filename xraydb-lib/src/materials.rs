//! Compound-material attenuation and refractive index.

use std::f64::consts::PI;

use crate::constants::{AVOGADRO, PLANCK_HC, R_ELECTRON_CM};
use crate::db::XrayDb;
use crate::elam::{CrossSectionKind, clamp_log_energy};
use crate::error::{Result, XrayDbError};

/// Complex refractive index decrements for a material at one energy.
///
/// The index of refraction is `n = 1 − delta − i·beta`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefractiveIndex {
    /// Real decrement δ (dimensionless, typically ~1e-5 to 1e-6).
    pub delta: f64,
    /// Imaginary part β (dimensionless).
    pub beta: f64,
    /// 1/e attenuation length in cm; infinite when β is zero.
    pub attenuation_length_cm: f64,
}

fn validate_density(density: f64) -> Result<()> {
    if density.is_finite() && density > 0.0 {
        return Ok(());
    }
    Err(XrayDbError::invalid_input(format!(
        "density must be finite and > 0 g/cm³, got {density}"
    )))
}

impl XrayDb {
    /// Linear attenuation coefficient µ in 1/cm for a compound, at one energy.
    ///
    /// * `formula` — any formula accepted by [`XrayDb::parse_formula`]
    /// * `density` — g/cm³
    /// * `energy` — eV, clamped to the Elam range
    ///
    /// ```
    /// # use xraydb::{XrayDb, CrossSectionKind};
    /// let db = XrayDb::try_new()?;
    /// let mu = db.material_mu_at("H2O", 1.0, 10_000.0, CrossSectionKind::Total)?;
    /// assert!((mu - 5.33).abs() < 0.05, "got {mu}");
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn material_mu_at(
        &self,
        formula: &str,
        density: f64,
        energy: f64,
        kind: CrossSectionKind,
    ) -> Result<f64> {
        validate_density(density)?;
        let fractions = self.mass_fractions(formula)?;
        self.mu_from_fractions_at(&fractions, density, energy, kind)
    }

    /// Linear attenuation coefficient µ in 1/cm for a compound, at several energies.
    ///
    /// See [`XrayDb::material_mu_at`].
    pub fn material_mu(
        &self,
        formula: &str,
        density: f64,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        validate_density(density)?;
        let fractions = self.mass_fractions(formula)?;
        self.mu_from_fractions(&fractions, density, energies, kind)
    }

    /// Linear attenuation coefficient µ in 1/cm from elemental **mass** fractions.
    ///
    /// Fractions need not sum to 1; they are normalized. Duplicate symbols are summed.
    pub fn material_mu_from_mass_fractions(
        &self,
        composition: &[(&str, f64)],
        density: f64,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        validate_density(density)?;
        let fractions = self.normalize_mass_fractions(composition)?;
        self.mu_from_fractions(&fractions, density, energies, kind)
    }

    /// Normalized elemental mass fractions for a formula.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let f = db.mass_fractions("H2O")?;
    /// let o = f.iter().find(|(s, _)| s == "O").map(|&(_, v)| v).unwrap_or(0.0);
    /// assert!((o - 0.888).abs() < 0.001, "got {o}");
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    ///
    /// # Weight-percent mixtures
    ///
    /// A formula may also describe a mixture by weight, e.g. `"Ru1wt%SiO2"` for 1 wt%
    /// ruthenium in silica. Each `wt%` term takes the percentage written immediately
    /// before it; the trailing component is the balance.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let f = db.mass_fractions("Ru1wt%SiO2")?;
    /// let ru = f.iter().find(|(s, _)| s == "Ru").map(|&(_, v)| v).unwrap_or(0.0);
    /// assert!((ru - 0.01).abs() < 1e-9, "got {ru}");
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn mass_fractions(&self, formula: &str) -> Result<Vec<(String, f64)>> {
        if formula.contains("wt%") {
            return self.weight_percent_mass_fractions(formula);
        }
        self.stoichiometric_mass_fractions(formula)
    }

    /// Mass fractions for a `wt%` mixture such as `"Ru1wt%SiO2"`.
    fn weight_percent_mass_fractions(&self, formula: &str) -> Result<Vec<(String, f64)>> {
        let compact: String = formula.chars().filter(|c| !c.is_whitespace()).collect();

        let mut components: Vec<(String, Option<f64>)> = Vec::new();
        let mut rest = compact.as_str();
        while let Some(marker) = rest.find("wt%") {
            let (head, tail) = rest.split_at(marker);
            // The percentage is the numeric literal immediately before "wt%".
            let split = head
                .rfind(|c: char| !(c.is_ascii_digit() || c == '.'))
                .map_or(0, |i| i + 1);
            let (sub_formula, percent_text) = head.split_at(split);
            let percent: f64 = percent_text.parse().map_err(|_| {
                XrayDbError::invalid_formula(
                    formula,
                    format!("expected a number before 'wt%', found '{percent_text}'"),
                )
            })?;
            if !percent.is_finite() || percent < 0.0 {
                return Err(XrayDbError::invalid_formula(
                    formula,
                    format!("weight percent must be finite and >= 0, got {percent}"),
                ));
            }
            components.push((sub_formula.to_string(), Some(percent)));
            rest = &tail["wt%".len()..];
        }
        if !rest.is_empty() {
            components.push((rest.to_string(), None));
        }

        let named_total: f64 = components.iter().filter_map(|&(_, p)| p).sum();
        let balance_count = components.iter().filter(|(_, p)| p.is_none()).count();
        if balance_count > 1 {
            return Err(XrayDbError::invalid_formula(
                formula,
                "a wt% mixture may have at most one balance component",
            ));
        }
        if named_total > 100.0 + 1e-9 {
            return Err(XrayDbError::invalid_formula(
                formula,
                format!("weight percentages sum to {named_total}%, which exceeds 100%"),
            ));
        }

        let balance_weight = if balance_count == 1 {
            100.0 - named_total
        } else {
            0.0
        };

        let mut totals: Vec<(String, f64)> = Vec::new();
        let mut grand_total = 0.0_f64;
        for (sub_formula, percent) in &components {
            let weight = percent.unwrap_or(balance_weight) / 100.0;
            if weight <= 0.0 {
                continue;
            }
            for (symbol, fraction) in self.stoichiometric_mass_fractions(sub_formula)? {
                let scaled = weight * fraction;
                match totals.iter_mut().find(|(s, _)| *s == symbol) {
                    Some((_, acc)) => *acc += scaled,
                    None => totals.push((symbol, scaled)),
                }
                grand_total += scaled;
            }
        }

        if !grand_total.is_finite() || grand_total <= 0.0 {
            return Err(XrayDbError::invalid_formula(
                formula,
                "wt% mixture has zero total weight",
            ));
        }
        totals.sort_by(|a, b| a.0.cmp(&b.0));
        for (_, fraction) in &mut totals {
            *fraction /= grand_total;
        }
        Ok(totals)
    }

    fn stoichiometric_mass_fractions(&self, formula: &str) -> Result<Vec<(String, f64)>> {
        let composition = self.parse_formula(formula)?;

        let mut masses: Vec<(String, f64)> = Vec::with_capacity(composition.len());
        let mut total_mass = 0.0_f64;
        for (symbol, count) in &composition {
            if !count.is_finite() || *count < 0.0 {
                return Err(XrayDbError::invalid_formula(
                    formula,
                    format!("invalid stoichiometry for {symbol}: {count}"),
                ));
            }
            if *count == 0.0 {
                continue;
            }
            let mass = count * self.molar_mass(symbol)?;
            if !mass.is_finite() || mass < 0.0 {
                return Err(XrayDbError::invalid_formula(
                    formula,
                    format!("invalid mass contribution for {symbol}"),
                ));
            }
            masses.push((symbol.clone(), mass));
            total_mass += mass;
        }

        if !total_mass.is_finite() || total_mass <= 0.0 {
            return Err(XrayDbError::invalid_formula(formula, "zero total mass"));
        }
        for (_, mass) in &mut masses {
            *mass /= total_mass;
        }
        Ok(masses)
    }

    /// Refractive index decrements for a compound at one energy.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let n = db.xray_delta_beta("SiO2", 2.2, 10_000.0)?;
    /// assert!(n.delta > 0.0 && n.beta > 0.0);
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn xray_delta_beta(
        &self,
        formula: &str,
        density: f64,
        energy: f64,
    ) -> Result<RefractiveIndex> {
        validate_density(density)?;

        // Go through mass fractions so this accepts exactly the same formulas as
        // `material_mu`, including `wt%` mixtures. Relative molar amounts are
        // `w_i / M_i`; only their ratios matter here.
        let fractions = self.mass_fractions(formula)?;
        let mut moles: Vec<(String, f64)> = Vec::with_capacity(fractions.len());
        for (symbol, weight) in fractions {
            let molar_mass = self.molar_mass(&symbol)?;
            if molar_mass.is_nan() || molar_mass <= 0.0 {
                return Err(XrayDbError::invalid_formula(
                    formula,
                    format!("element {symbol} has non-positive molar mass"),
                ));
            }
            moles.push((symbol, weight / molar_mass));
        }
        self.delta_beta_from_amounts(formula, &moles, density, energy)
    }

    fn delta_beta_from_amounts(
        &self,
        formula: &str,
        composition: &[(String, f64)],
        density: f64,
        energy: f64,
    ) -> Result<RefractiveIndex> {
        // Wavelength in cm.
        let wavelength = 1.0e-7 * PLANCK_HC / energy;

        let mut total_weight = 0.0_f64;
        let mut sum_f1 = 0.0_f64;
        let mut sum_f2 = 0.0_f64;

        // `composition` is sorted, so this sum is reproducible run to run.
        for (symbol, count) in composition {
            total_weight += count * self.molar_mass(symbol)?;

            let z = f64::from(self.atomic_number(symbol)?);
            // Chantler f1 is f′ — the anomalous correction with Z already subtracted.
            sum_f1 += count * (z + self.f1_chantler_at(symbol, energy)?);
            sum_f2 += count * self.f2_chantler_at(symbol, energy)?;
        }

        if !total_weight.is_finite() || total_weight <= 0.0 {
            return Err(XrayDbError::invalid_formula(formula, "zero total mass"));
        }

        let prefactor = R_ELECTRON_CM * wavelength * wavelength * density * AVOGADRO
            / (2.0 * PI * total_weight);

        let delta = prefactor * sum_f1;
        let beta = prefactor * sum_f2;
        let attenuation_length_cm = if beta > 0.0 {
            wavelength / (4.0 * PI * beta)
        } else {
            f64::INFINITY
        };

        Ok(RefractiveIndex {
            delta,
            beta,
            attenuation_length_cm,
        })
    }

    // ── Shared internals ──────────────────────────────────────────────────

    fn normalize_mass_fractions(&self, composition: &[(&str, f64)]) -> Result<Vec<(String, f64)>> {
        if composition.is_empty() {
            return Err(XrayDbError::invalid_input(
                "composition mass fractions must not be empty",
            ));
        }

        let mut by_symbol: Vec<(String, f64)> = Vec::with_capacity(composition.len());
        let mut total = 0.0_f64;
        for &(symbol, fraction) in composition {
            if !fraction.is_finite() || fraction < 0.0 {
                return Err(XrayDbError::invalid_input(format!(
                    "invalid mass fraction for '{symbol}': {fraction}"
                )));
            }
            if fraction == 0.0 {
                continue;
            }
            let canonical = self.symbol(symbol)?.to_string();
            match by_symbol.iter_mut().find(|(s, _)| *s == canonical) {
                Some((_, acc)) => *acc += fraction,
                None => by_symbol.push((canonical, fraction)),
            }
            total += fraction;
        }

        if !total.is_finite() || total <= 0.0 {
            return Err(XrayDbError::invalid_input(
                "composition mass fractions must include at least one positive value",
            ));
        }
        by_symbol.sort_by(|a, b| a.0.cmp(&b.0));
        for (_, fraction) in &mut by_symbol {
            *fraction /= total;
        }
        Ok(by_symbol)
    }

    pub(crate) fn mu_from_fractions_at(
        &self,
        fractions: &[(String, f64)],
        density: f64,
        energy: f64,
        kind: CrossSectionKind,
    ) -> Result<f64> {
        let log_energy = clamp_log_energy(energy);
        let mut mu = 0.0_f64;
        for (symbol, fraction) in fractions {
            if *fraction == 0.0 {
                continue;
            }
            mu += fraction * self.mu_elam_sym(symbol, symbol, log_energy, kind)?;
        }
        Ok(mu * density)
    }

    pub(crate) fn mu_from_fractions(
        &self,
        fractions: &[(String, f64)],
        density: f64,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        let mut mu = vec![0.0_f64; energies.len()];
        for (symbol, fraction) in fractions {
            if *fraction == 0.0 {
                continue;
            }
            let elem_mu = self.mu_elam(symbol, energies, kind)?;
            for (acc, &value) in mu.iter_mut().zip(elem_mu.iter()) {
                *acc += fraction * value;
            }
        }
        for value in &mut mu {
            *value *= density;
        }
        Ok(mu)
    }
}
