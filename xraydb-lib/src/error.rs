//! Error type for all fallible operations in this crate.

use thiserror::Error;

/// Errors returned by [`XrayDb`](crate::XrayDb) queries and calculations.
///
/// This enum is `#[non_exhaustive]`: new variants may be added in a minor release,
/// so match on it with a `_ => …` arm.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum XrayDbError {
    /// The element identifier did not resolve to a symbol, name, or atomic number.
    #[error("unknown element: {0}")]
    UnknownElement(String),

    /// The requested absorption edge does not exist for this element.
    #[error("unknown edge '{edge}' for element '{element}'{}", format_available(.available))]
    UnknownEdge {
        /// The element that was queried.
        element: String,
        /// The edge label that was requested.
        edge: String,
        /// Edge labels that *are* available for this element.
        available: Vec<String>,
    },

    /// The ion has no Waasmaier–Kirfel coefficients.
    #[error("unknown ion: {0}")]
    UnknownIon(String),

    /// The gas has no tabulated ionization potential.
    #[error("unknown gas: {gas}{}", format_available(.available))]
    UnknownGas {
        /// The gas name that was requested.
        gas: String,
        /// Gas names that *are* tabulated.
        available: Vec<String>,
    },

    /// A chemical formula could not be parsed.
    #[error("invalid chemical formula '{formula}': {reason}")]
    InvalidFormula {
        /// The formula as supplied by the caller.
        formula: String,
        /// Why it could not be parsed.
        reason: String,
    },

    /// An argument supplied by the caller was out of range or inconsistent.
    #[error("invalid input: {0}")]
    InvalidInput(String),

    /// The embedded database could not be decoded, or is internally inconsistent.
    #[error("data error: {0}")]
    DataError(String),
}

impl XrayDbError {
    /// Convenience constructor for [`XrayDbError::InvalidInput`].
    pub(crate) fn invalid_input(msg: impl Into<String>) -> Self {
        Self::InvalidInput(msg.into())
    }

    /// Convenience constructor for [`XrayDbError::InvalidFormula`].
    pub(crate) fn invalid_formula(formula: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidFormula {
            formula: formula.into(),
            reason: reason.into(),
        }
    }
}

/// Render an "available options" hint, truncated so an error message stays readable.
fn format_available(available: &[String]) -> String {
    if available.is_empty() {
        return String::new();
    }
    const MAX: usize = 12;
    let shown: Vec<&str> = available.iter().take(MAX).map(String::as_str).collect();
    if available.len() > MAX {
        format!(
            " (available: {}, … and {} more)",
            shown.join(", "),
            available.len() - MAX
        )
    } else {
        format!(" (available: {})", shown.join(", "))
    }
}

/// Result alias used throughout this crate.
pub type Result<T> = core::result::Result<T, XrayDbError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unknown_edge_lists_available_edges() {
        let err = XrayDbError::UnknownEdge {
            element: "Fe".into(),
            edge: "Q9".into(),
            available: vec!["K".into(), "L1".into()],
        };
        let msg = err.to_string();
        assert!(msg.contains("unknown edge 'Q9'"), "{msg}");
        assert!(msg.contains("available: K, L1"), "{msg}");
    }

    #[test]
    fn empty_available_adds_no_hint() {
        let err = XrayDbError::UnknownGas {
            gas: "argonn".into(),
            available: Vec::new(),
        };
        assert_eq!(err.to_string(), "unknown gas: argonn");
    }

    #[test]
    fn available_list_is_truncated() {
        let many: Vec<String> = (0..20).map(|i| format!("e{i}")).collect();
        let msg = XrayDbError::UnknownGas {
            gas: "x".into(),
            available: many,
        }
        .to_string();
        assert!(msg.contains("… and 8 more"), "{msg}");
    }
}
