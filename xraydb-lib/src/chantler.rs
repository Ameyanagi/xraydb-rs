use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};
use crate::interp::{interp, interp_loglog};

/// Kind of Chantler cross-section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChantlerKind {
    Total,
    Photo,
    Incoherent,
}

impl XrayDb {
    fn chantler_record(&self, element: &str) -> Result<&xraydb_data::ChantlerRecord> {
        let sym = self.symbol(element)?;
        self.chantler_by_symbol(sym)
            .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))
    }

    /// Returns tabulated energy points for Chantler data for an element.
    ///
    /// Optionally filtered to a range [emin, emax] in eV.
    pub fn chantler_energies(
        &self,
        element: &str,
        emin: Option<f64>,
        emax: Option<f64>,
    ) -> Result<Vec<f64>> {
        let row = self.chantler_record(element)?;
        let emin = emin.unwrap_or(0.0);
        let emax = emax.unwrap_or(1e9);
        Ok(row
            .energy
            .iter()
            .copied()
            .filter(|&e| e >= emin && e <= emax)
            .collect())
    }

    /// Returns f1 — real part of anomalous X-ray scattering factor (Chantler).
    ///
    /// Uses linear interpolation (matching Python's UnivariateSpline with s=0).
    pub fn f1_chantler(&self, element: &str, energies: &[f64]) -> Result<Vec<f64>> {
        let row = self.chantler_record(element)?;
        let (emin, emax) = chantler_energy_bounds(row);

        // Clamp energies to valid range
        let clamped = clamp_energies(energies, emin, emax);

        // For f1, use linear interpolation in linear space
        // (Python uses UnivariateSpline; linear interp is a reasonable approximation)
        Ok(interp(&clamped, &row.energy, &row.f1))
    }

    /// Returns f2 — imaginary part of anomalous X-ray scattering factor (Chantler).
    ///
    /// Uses log-log linear interpolation.
    pub fn f2_chantler(&self, element: &str, energies: &[f64]) -> Result<Vec<f64>> {
        let row = self.chantler_record(element)?;
        let (emin, emax) = chantler_energy_bounds(row);

        let clamped = clamp_energies(energies, emin, emax);

        // Clamp values to avoid log(0)
        let f2_safe = safe_for_log(&row.f2);
        Ok(interp_loglog(&clamped, &row.energy, &f2_safe))
    }

    /// Returns X-ray mass attenuation coefficient (mu/rho) in cm²/g (Chantler).
    ///
    /// Uses log-log linear interpolation.
    pub fn mu_chantler(
        &self,
        element: &str,
        energies: &[f64],
        kind: ChantlerKind,
    ) -> Result<Vec<f64>> {
        let row = self.chantler_record(element)?;
        let (emin, emax) = chantler_energy_bounds(row);

        let clamped = clamp_energies(energies, emin, emax);

        let values = match kind {
            ChantlerKind::Total => &row.mu_total,
            ChantlerKind::Photo => &row.mu_photo,
            ChantlerKind::Incoherent => &row.mu_incoh,
        };

        // Clamp values to avoid log(0)
        let safe = safe_for_log(values);
        Ok(interp_loglog(&clamped, &row.energy, &safe))
    }
}

#[inline]
fn chantler_energy_bounds(row: &xraydb_data::ChantlerRecord) -> (f64, f64) {
    let emin = row.energy.first().copied().unwrap_or(0.0);
    let emax = row.energy.last().copied().unwrap_or(emin).min(1e6_f64);
    (emin, emax)
}

#[inline]
fn clamp_energies(energies: &[f64], emin: f64, emax: f64) -> Vec<f64> {
    energies.iter().map(|&e| e.clamp(emin, emax)).collect()
}

#[inline]
fn safe_for_log(values: &[f64]) -> Vec<f64> {
    values
        .iter()
        .map(|&v| if v.abs() < 1e-99 { 1e-99 } else { v })
        .collect()
}
