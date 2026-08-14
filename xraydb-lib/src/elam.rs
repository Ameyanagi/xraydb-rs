//! Elam/Ravel/Sieber cross-section tables.

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};
use crate::spline::elam_spline_one;

/// Which cross-section to evaluate from the Elam tables.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossSectionKind {
    /// Photoelectric absorption.
    Photo,
    /// Coherent (Rayleigh) scattering.
    Coherent,
    /// Incoherent (Compton) scattering.
    Incoherent,
    /// Sum of the three above.
    Total,
}

/// Lowest energy covered by the Elam tables, in eV.
pub const ELAM_ENERGY_MIN: f64 = 100.0;
/// Highest energy covered by the Elam tables, in eV.
pub const ELAM_ENERGY_MAX: f64 = 800_000.0;

impl XrayDb {
    /// Energy range `(min, max)` in eV covered by the Elam tables.
    ///
    /// Queries outside this range are silently clamped to the endpoints.
    pub const fn elam_energy_range() -> (f64, f64) {
        (ELAM_ENERGY_MIN, ELAM_ENERGY_MAX)
    }

    /// Elam mass attenuation cross-section µ/ρ in cm²/g at one energy.
    ///
    /// Energies are in eV and are clamped to [`ELAM_ENERGY_MIN`], [`ELAM_ENERGY_MAX`].
    ///
    /// ```
    /// # use xraydb::{XrayDb, CrossSectionKind};
    /// let db = XrayDb::try_new()?;
    /// let mu = db.mu_elam_at("Fe", 10_000.0, CrossSectionKind::Total)?;
    /// assert!((mu - 170.6).abs() < 1.0, "got {mu}");
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn mu_elam_at(&self, element: &str, energy: f64, kind: CrossSectionKind) -> Result<f64> {
        let sym = self.symbol(element)?;
        self.mu_elam_sym(sym, element, clamp_log_energy(energy), kind)
    }

    /// Elam µ/ρ in cm²/g at each of several energies.
    ///
    /// See [`XrayDb::mu_elam_at`]. Resolving the element and looking up its table rows
    /// happens once for the whole batch.
    pub fn mu_elam(
        &self,
        element: &str,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        let sym = self.symbol(element)?;
        match kind {
            CrossSectionKind::Total => {
                let photo = self.photo_row(sym, element)?;
                let scatter = self.scatter_row(sym, element)?;
                Ok(energies
                    .iter()
                    .map(|&e| {
                        let le = clamp_log_energy(e);
                        total_at(photo, scatter, le)
                    })
                    .collect())
            }
            CrossSectionKind::Photo => {
                let row = self.photo_row(sym, element)?;
                Ok(energies
                    .iter()
                    .map(|&e| {
                        elam_spline_one(
                            &row.log_energy,
                            &row.log_photoabsorption,
                            &row.log_photoabsorption_spline,
                            clamp_log_energy(e),
                        )
                        .exp()
                    })
                    .collect())
            }
            CrossSectionKind::Coherent | CrossSectionKind::Incoherent => {
                let row = self.scatter_row(sym, element)?;
                let (y, spl) = if kind == CrossSectionKind::Coherent {
                    (&row.log_coherent_scatter, &row.log_coherent_scatter_spline)
                } else {
                    (
                        &row.log_incoherent_scatter,
                        &row.log_incoherent_scatter_spline,
                    )
                };
                Ok(energies
                    .iter()
                    .map(|&e| elam_spline_one(&row.log_energy, y, spl, clamp_log_energy(e)).exp())
                    .collect())
            }
        }
    }

    fn photo_row(
        &self,
        sym: &str,
        element: &str,
    ) -> Result<&'static xraydb_data::PhotoabsorptionRecord> {
        self.photo_by_symbol(sym)
            .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))
    }

    fn scatter_row(
        &self,
        sym: &str,
        element: &str,
    ) -> Result<&'static xraydb_data::ScatteringRecord> {
        self.scatter_by_symbol(sym)
            .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))
    }

    /// Shared scalar evaluation against an already-resolved symbol and log-energy.
    pub(crate) fn mu_elam_sym(
        &self,
        sym: &str,
        element: &str,
        log_energy: f64,
        kind: CrossSectionKind,
    ) -> Result<f64> {
        match kind {
            CrossSectionKind::Total => {
                let photo = self.photo_row(sym, element)?;
                let scatter = self.scatter_row(sym, element)?;
                Ok(total_at(photo, scatter, log_energy))
            }
            CrossSectionKind::Photo => {
                let row = self.photo_row(sym, element)?;
                Ok(elam_spline_one(
                    &row.log_energy,
                    &row.log_photoabsorption,
                    &row.log_photoabsorption_spline,
                    log_energy,
                )
                .exp())
            }
            CrossSectionKind::Coherent => {
                let row = self.scatter_row(sym, element)?;
                Ok(elam_spline_one(
                    &row.log_energy,
                    &row.log_coherent_scatter,
                    &row.log_coherent_scatter_spline,
                    log_energy,
                )
                .exp())
            }
            CrossSectionKind::Incoherent => {
                let row = self.scatter_row(sym, element)?;
                Ok(elam_spline_one(
                    &row.log_energy,
                    &row.log_incoherent_scatter,
                    &row.log_incoherent_scatter_spline,
                    log_energy,
                )
                .exp())
            }
        }
    }
}

fn total_at(
    photo: &xraydb_data::PhotoabsorptionRecord,
    scatter: &xraydb_data::ScatteringRecord,
    log_energy: f64,
) -> f64 {
    let p = elam_spline_one(
        &photo.log_energy,
        &photo.log_photoabsorption,
        &photo.log_photoabsorption_spline,
        log_energy,
    );
    let c = elam_spline_one(
        &scatter.log_energy,
        &scatter.log_coherent_scatter,
        &scatter.log_coherent_scatter_spline,
        log_energy,
    );
    let i = elam_spline_one(
        &scatter.log_energy,
        &scatter.log_incoherent_scatter,
        &scatter.log_incoherent_scatter_spline,
        log_energy,
    );
    p.exp() + c.exp() + i.exp()
}

#[inline]
pub(crate) fn clamp_log_energy(energy: f64) -> f64 {
    energy.clamp(ELAM_ENERGY_MIN, ELAM_ENERGY_MAX).ln()
}
