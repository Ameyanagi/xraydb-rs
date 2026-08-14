//! Coster–Kronig transition probabilities.

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

impl XrayDb {
    /// Coster–Kronig transition probability between two levels of the same shell.
    ///
    /// With `total = true` this is the total probability including routes via
    /// intermediate states; otherwise it is the direct transition probability.
    pub fn ck_probability(
        &self,
        element: &str,
        initial: &str,
        final_level: &str,
        total: bool,
    ) -> Result<f64> {
        let sym = self.symbol(element)?;
        let record = self
            .coster_kronig_record(sym, initial, final_level)
            .ok_or_else(|| XrayDbError::UnknownEdge {
                element: element.to_string(),
                edge: format!("{initial}->{final_level}"),
                available: Vec::new(),
            })?;

        Ok(if total {
            record.total_transition_probability
        } else {
            record.transition_probability
        })
    }
}
