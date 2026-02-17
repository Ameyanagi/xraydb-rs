pub mod chantler;
pub mod chemparser;
pub mod compton;
pub mod constants;
pub mod core_widths;
pub mod coster_kronig;
pub mod db;
pub mod elam;
pub mod error;
pub mod interp;
pub mod ionchamber;
pub mod materials;
pub(crate) mod materials_db;
#[cfg(feature = "optics")]
pub mod optics;
pub mod spline;
pub mod transitions;
pub mod waasmaier;

pub use chantler::ChantlerKind;
pub use compton::ComptonEnergies;
pub use db::XrayDb;
pub use elam::CrossSectionKind;
pub use error::{Result, XrayDbError};
pub use ionchamber::IonChamberFluxes;
#[cfg(feature = "optics")]
pub use optics::{DarwinWidth, Polarization};
pub use transitions::{XrayEdge, XrayLine};
pub use xraydb_data;
