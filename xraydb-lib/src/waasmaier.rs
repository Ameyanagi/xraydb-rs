use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

impl XrayDb {
    /// Returns list of supported ion names for f0 calculations.
    ///
    /// If `element` is provided, returns only ions for that element.
    pub fn f0_ions(&self, element: Option<&str>) -> Result<Vec<&str>> {
        let ions: Vec<&str> = match element {
            Some(elem) => {
                let sym = self.symbol(elem)?;
                self.waasmaier_indices_by_symbol(sym)
                    .unwrap_or(&[])
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

    /// Returns f0 elastic X-ray scattering factor for an ion at given q values.
    ///
    /// q = sin(theta) / lambda in Angstroms^-1.
    ///
    /// Formula: f0(q) = c + sum_i(a_i * exp(-b_i * q^2))
    /// where c = offset, a_i = scale, b_i = exponents.
    pub fn f0(&self, ion: &str, q: &[f64]) -> Result<Vec<f64>> {
        let record = self
            .waasmaier_by_ion(ion)
            .ok_or_else(|| XrayDbError::UnknownIon(ion.to_string()))?;

        Ok(q.iter()
            .map(|&qi| {
                let q2 = qi * qi;
                let mut val = record.offset;
                for (a, b) in record.scale.iter().zip(record.exponents.iter()) {
                    val += a * (-b * q2).exp();
                }
                val
            })
            .collect())
    }
}
