use std::fmt;

#[derive(Debug)]
pub enum XrayDbError {
    UnknownElement(String),
    UnknownEdge { element: String, edge: String },
    UnknownIon(String),
    UnknownGas(String),
    EnergyOutOfRange { energy: f64, min: f64, max: f64 },
    InvalidFormula(String),
    DataError(String),
}

pub type Result<T> = std::result::Result<T, XrayDbError>;

impl fmt::Display for XrayDbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnknownElement(e) => write!(f, "unknown element: {e}"),
            Self::UnknownEdge { element, edge } => {
                write!(f, "unknown edge '{edge}' for element '{element}'")
            }
            Self::UnknownIon(ion) => write!(f, "unknown ion: {ion}"),
            Self::UnknownGas(gas) => write!(f, "unknown gas: {gas}"),
            Self::EnergyOutOfRange { energy, min, max } => {
                write!(f, "energy {energy} eV out of range [{min}, {max}]")
            }
            Self::InvalidFormula(formula) => write!(f, "invalid chemical formula: {formula}"),
            Self::DataError(msg) => write!(f, "data error: {msg}"),
        }
    }
}

impl std::error::Error for XrayDbError {}
