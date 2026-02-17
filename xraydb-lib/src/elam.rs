use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};
use crate::spline::elam_spline;

/// Kind of cross-section for Elam calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossSectionKind {
    Photo,
    Coherent,
    Incoherent,
    Total,
}

impl XrayDb {
    /// Returns mass attenuation cross-section in cm²/g using Elam tables.
    ///
    /// Energies are in eV. Values are clamped to [100 eV, 800 keV].
    pub fn mu_elam(
        &self,
        element: &str,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        let sym = self.symbol(element)?;
        let log_en = clamp_log_energies(energies);

        match kind {
            CrossSectionKind::Total => {
                let photo_row = self
                    .photo_by_symbol(sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;
                let scatter_row = self
                    .scatter_by_symbol(sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;

                let photo_log = elam_spline(
                    &photo_row.log_energy,
                    &photo_row.log_photoabsorption,
                    &photo_row.log_photoabsorption_spline,
                    &log_en,
                );
                let coh_log = elam_spline(
                    &scatter_row.log_energy,
                    &scatter_row.log_coherent_scatter,
                    &scatter_row.log_coherent_scatter_spline,
                    &log_en,
                );
                let incoh_log = elam_spline(
                    &scatter_row.log_energy,
                    &scatter_row.log_incoherent_scatter,
                    &scatter_row.log_incoherent_scatter_spline,
                    &log_en,
                );

                Ok(photo_log
                    .iter()
                    .zip(coh_log.iter())
                    .zip(incoh_log.iter())
                    .map(|((&p, &c), &i)| p.exp() + c.exp() + i.exp())
                    .collect())
            }
            other => self.cross_section_elam_with_symbol(sym, element, &log_en, other),
        }
    }

    /// Returns Elam cross-section for a specific kind (photo, coh, or incoh).
    fn cross_section_elam_with_symbol(
        &self,
        sym: &str,
        element: &str,
        log_en: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        match kind {
            CrossSectionKind::Photo => {
                let row = self
                    .photo_by_symbol(sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;

                let result = elam_spline(
                    &row.log_energy,
                    &row.log_photoabsorption,
                    &row.log_photoabsorption_spline,
                    log_en,
                );
                Ok(result.into_iter().map(|v| v.exp()).collect())
            }
            CrossSectionKind::Coherent => {
                let row = self
                    .scatter_by_symbol(sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;

                let result = elam_spline(
                    &row.log_energy,
                    &row.log_coherent_scatter,
                    &row.log_coherent_scatter_spline,
                    log_en,
                );
                Ok(result.into_iter().map(|v| v.exp()).collect())
            }
            CrossSectionKind::Incoherent => {
                let row = self
                    .scatter_by_symbol(sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;

                let result = elam_spline(
                    &row.log_energy,
                    &row.log_incoherent_scatter,
                    &row.log_incoherent_scatter_spline,
                    log_en,
                );
                Ok(result.into_iter().map(|v| v.exp()).collect())
            }
            CrossSectionKind::Total => unreachable!(),
        }
    }
}

#[inline]
fn clamp_log_energies(energies: &[f64]) -> Vec<f64> {
    energies
        .iter()
        .map(|&e| e.clamp(100.0, 800_000.0).ln())
        .collect()
}
