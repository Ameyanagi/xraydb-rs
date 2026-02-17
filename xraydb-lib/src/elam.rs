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
    pub fn mu_elam(&self, element: &str, energies: &[f64], kind: CrossSectionKind) -> Result<Vec<f64>> {
        match kind {
            CrossSectionKind::Total => {
                let photo = self.cross_section_elam(element, energies, CrossSectionKind::Photo)?;
                let coh = self.cross_section_elam(element, energies, CrossSectionKind::Coherent)?;
                let incoh = self.cross_section_elam(element, energies, CrossSectionKind::Incoherent)?;
                Ok(photo
                    .iter()
                    .zip(coh.iter())
                    .zip(incoh.iter())
                    .map(|((&p, &c), &i)| p + c + i)
                    .collect())
            }
            other => self.cross_section_elam(element, energies, other),
        }
    }

    /// Returns Elam cross-section for a specific kind (photo, coh, or incoh).
    fn cross_section_elam(
        &self,
        element: &str,
        energies: &[f64],
        kind: CrossSectionKind,
    ) -> Result<Vec<f64>> {
        let sym = self.symbol(element)?;

        // Clamp energies to valid range [100 eV, 800 keV]
        let clamped: Vec<f64> = energies
            .iter()
            .map(|&e| e.clamp(100.0, 800_000.0))
            .collect();

        let log_en: Vec<f64> = clamped.iter().map(|e| e.ln()).collect();

        match kind {
            CrossSectionKind::Photo => {
                let row = self
                    .raw()
                    .photoabsorption
                    .iter()
                    .find(|r| r.element == sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;

                let result = elam_spline(
                    &row.log_energy,
                    &row.log_photoabsorption,
                    &row.log_photoabsorption_spline,
                    &log_en,
                );
                Ok(result.into_iter().map(|v| v.exp()).collect())
            }
            CrossSectionKind::Coherent => {
                let row = self
                    .raw()
                    .scattering
                    .iter()
                    .find(|r| r.element == sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;

                let result = elam_spline(
                    &row.log_energy,
                    &row.log_coherent_scatter,
                    &row.log_coherent_scatter_spline,
                    &log_en,
                );
                Ok(result.into_iter().map(|v| v.exp()).collect())
            }
            CrossSectionKind::Incoherent => {
                let row = self
                    .raw()
                    .scattering
                    .iter()
                    .find(|r| r.element == sym)
                    .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))?;

                let result = elam_spline(
                    &row.log_energy,
                    &row.log_incoherent_scatter,
                    &row.log_incoherent_scatter_spline,
                    &log_en,
                );
                Ok(result.into_iter().map(|v| v.exp()).collect())
            }
            CrossSectionKind::Total => unreachable!(),
        }
    }
}
