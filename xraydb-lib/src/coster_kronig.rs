use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

impl XrayDb {
    /// Returns Coster-Kronig transition probability.
    ///
    /// If `total` is true, returns the total transition probability
    /// (including via intermediate states). Otherwise returns the
    /// direct transition probability.
    pub fn ck_probability(
        &self,
        element: &str,
        initial: &str,
        final_level: &str,
        total: bool,
    ) -> Result<f64> {
        let sym = self.symbol(element)?;
        let record = self
            .raw()
            .coster_kronig
            .iter()
            .find(|ck| {
                ck.element == sym && ck.initial_level == initial && ck.final_level == final_level
            })
            .ok_or_else(|| XrayDbError::UnknownEdge {
                element: element.to_string(),
                edge: format!("{initial}->{final_level}"),
            })?;

        if total {
            Ok(record.total_transition_probability)
        } else {
            Ok(record.transition_probability)
        }
    }
}
