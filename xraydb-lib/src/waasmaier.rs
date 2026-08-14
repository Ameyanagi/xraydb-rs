//! Waasmaier–Kirfel elastic (Thomson) scattering factors, f0.

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

impl XrayDb {
    /// Ion names with tabulated f0 coefficients.
    ///
    /// If `element` is given, only that element's ions are returned (e.g. `Fe`, `Fe2+`,
    /// `Fe3+`); otherwise every tabulated ion is returned.
    pub fn f0_ions(&self, element: Option<&str>) -> Result<Vec<&'static str>> {
        let ions: Vec<&'static str> = match element {
            Some(elem) => {
                let sym = self.symbol(elem)?;
                self.waasmaier_indices_by_symbol(sym)
                    .iter()
                    .filter_map(|&idx| self.raw().waasmaier.get(idx).map(|w| w.ion.as_str()))
                    .collect()
            }
            None => self
                .raw()
                .waasmaier
                .iter()
                .map(|w| w.ion.as_str())
                .collect(),
        };
        Ok(ions)
    }

    /// f0 elastic scattering factor for an ion at one momentum transfer.
    ///
    /// `q` is `sin(θ)/λ` in Å⁻¹, valid for `q` up to about 6 Å⁻¹.
    ///
    /// f0(q) = c + Σᵢ aᵢ·exp(−bᵢ·q²)
    pub fn f0_at(&self, ion: &str, q: f64) -> Result<f64> {
        let record = self
            .waasmaier_by_ion(ion)
            .ok_or_else(|| XrayDbError::UnknownIon(ion.to_string()))?;
        Ok(f0_from(record, q))
    }

    /// f0 at each of several `q` values. See [`XrayDb::f0_at`].
    pub fn f0(&self, ion: &str, q: &[f64]) -> Result<Vec<f64>> {
        let record = self
            .waasmaier_by_ion(ion)
            .ok_or_else(|| XrayDbError::UnknownIon(ion.to_string()))?;
        Ok(q.iter().map(|&qi| f0_from(record, qi)).collect())
    }
}

#[inline]
fn f0_from(record: &xraydb_data::WaasmaierRecord, q: f64) -> f64 {
    let q2 = q * q;
    let mut val = record.offset;
    for (a, b) in record.scale.iter().zip(record.exponents.iter()) {
        val += a * (-b * q2).exp();
    }
    val
}
