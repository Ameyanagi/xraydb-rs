//! Chantler tables: anomalous scattering factors and mass attenuation coefficients.

use crate::db::{ChantlerLogs, XrayDb};
use crate::error::{Result, XrayDbError};
use crate::interp::{interp_loglog_pre, interp_one};

/// Which Chantler mass attenuation coefficient to evaluate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChantlerKind {
    /// Total attenuation.
    Total,
    /// Photoelectric absorption only.
    Photo,
    /// Incoherent (Compton) scattering only.
    Incoherent,
}

/// Upper bound applied to Chantler energies, matching upstream XrayDB.
const CHANTLER_EMAX_CAP: f64 = 1e6;

impl XrayDb {
    fn chantler_record(
        &self,
        element: &str,
    ) -> Result<(&'static xraydb_data::ChantlerRecord, &'static ChantlerLogs)> {
        let sym = self.symbol(element)?;
        self.chantler_by_symbol(sym)
            .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))
    }

    /// Tabulated Chantler energy grid for an element, in eV.
    ///
    /// Optionally filtered to `[emin, emax]`.
    pub fn chantler_energies(
        &self,
        element: &str,
        emin: Option<f64>,
        emax: Option<f64>,
    ) -> Result<Vec<f64>> {
        let (row, _) = self.chantler_record(element)?;
        let emin = emin.unwrap_or(0.0);
        let emax = emax.unwrap_or(1e9);
        Ok(row
            .energy
            .iter()
            .copied()
            .filter(|&e| e >= emin && e <= emax)
            .collect())
    }

    /// Valid energy range `(min, max)` in eV for this element's Chantler table.
    ///
    /// Queries outside this range are silently clamped to the endpoints; use this to
    /// detect that before it happens.
    pub fn chantler_energy_range(&self, element: &str) -> Result<(f64, f64)> {
        let (row, _) = self.chantler_record(element)?;
        Ok(energy_bounds(row))
    }

    /// f1 — the real part of the anomalous scattering factor, at one energy.
    ///
    /// # Note on accuracy
    ///
    /// Upstream XrayDB evaluates f1 with a smoothing spline (`UnivariateSpline`, `s=0`);
    /// this port uses linear interpolation on the same grid. Agreement is close between
    /// tabulated points and exact at them, but the two are not bit-identical.
    ///
    /// The tabulated quantity is f′ (the anomalous correction with Z already subtracted),
    /// not the full f1 — add Z for the total.
    ///
    /// Energies outside the tabulated range are clamped to its endpoints; see
    /// [`XrayDb::chantler_energy_range`].
    pub fn f1_chantler_at(&self, element: &str, energy: f64) -> Result<f64> {
        let (row, _) = self.chantler_record(element)?;
        let (emin, emax) = energy_bounds(row);
        Ok(interp_one(energy.clamp(emin, emax), &row.energy, &row.f1))
    }

    /// f1 at each of several energies. See [`XrayDb::f1_chantler_at`].
    pub fn f1_chantler(&self, element: &str, energies: &[f64]) -> Result<Vec<f64>> {
        let (row, _) = self.chantler_record(element)?;
        let (emin, emax) = energy_bounds(row);
        Ok(energies
            .iter()
            .map(|&e| interp_one(e.clamp(emin, emax), &row.energy, &row.f1))
            .collect())
    }

    /// f2 — the imaginary part of the anomalous scattering factor, at one energy.
    ///
    /// Evaluated by log-log linear interpolation. Energies outside the tabulated range
    /// are clamped; see [`XrayDb::chantler_energy_range`].
    pub fn f2_chantler_at(&self, element: &str, energy: f64) -> Result<f64> {
        let (row, logs) = self.chantler_record(element)?;
        let (emin, emax) = energy_bounds(row);
        Ok(interp_loglog_pre(
            energy.clamp(emin, emax),
            &logs.log_energy,
            &logs.log_f2,
        ))
    }

    /// f2 at each of several energies. See [`XrayDb::f2_chantler_at`].
    pub fn f2_chantler(&self, element: &str, energies: &[f64]) -> Result<Vec<f64>> {
        let (row, logs) = self.chantler_record(element)?;
        let (emin, emax) = energy_bounds(row);
        Ok(energies
            .iter()
            .map(|&e| interp_loglog_pre(e.clamp(emin, emax), &logs.log_energy, &logs.log_f2))
            .collect())
    }

    /// Chantler mass attenuation coefficient µ/ρ in cm²/g, at one energy.
    ///
    /// Evaluated by log-log linear interpolation. Energies outside the tabulated range
    /// are clamped; see [`XrayDb::chantler_energy_range`].
    pub fn mu_chantler_at(&self, element: &str, energy: f64, kind: ChantlerKind) -> Result<f64> {
        let (row, logs) = self.chantler_record(element)?;
        let (emin, emax) = energy_bounds(row);
        Ok(interp_loglog_pre(
            energy.clamp(emin, emax),
            &logs.log_energy,
            log_values(logs, kind),
        ))
    }

    /// Chantler µ/ρ at each of several energies. See [`XrayDb::mu_chantler_at`].
    pub fn mu_chantler(
        &self,
        element: &str,
        energies: &[f64],
        kind: ChantlerKind,
    ) -> Result<Vec<f64>> {
        let (row, logs) = self.chantler_record(element)?;
        let (emin, emax) = energy_bounds(row);
        let values = log_values(logs, kind);
        Ok(energies
            .iter()
            .map(|&e| interp_loglog_pre(e.clamp(emin, emax), &logs.log_energy, values))
            .collect())
    }
}

#[inline]
fn log_values(logs: &ChantlerLogs, kind: ChantlerKind) -> &[f64] {
    match kind {
        ChantlerKind::Total => &logs.log_total,
        ChantlerKind::Photo => &logs.log_photo,
        ChantlerKind::Incoherent => &logs.log_incoh,
    }
}

#[inline]
fn energy_bounds(row: &xraydb_data::ChantlerRecord) -> (f64, f64) {
    let emin = row.energy.first().copied().unwrap_or(0.0);
    let emax = row
        .energy
        .last()
        .copied()
        .unwrap_or(emin)
        .min(CHANTLER_EMAX_CAP);
    (emin, emax)
}
