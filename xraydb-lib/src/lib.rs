//! X-ray reference data for the elements — a pure-Rust port of
//! [XrayDB](https://github.com/xraypy/XrayDB).
//!
//! The complete database (elements, absorption edges, emission lines, Elam and Chantler
//! cross-sections, Waasmaier–Kirfel form factors) is compiled into the binary as a
//! 3 MB zstd blob and decoded lazily on first use. There are no runtime files to ship
//! and no network access.
//!
//! # Getting started
//!
//! ```
//! use xraydb::{CrossSectionKind, XrayDb};
//!
//! let db = XrayDb::try_new()?;
//!
//! // Element facts
//! assert_eq!(db.atomic_number("Fe")?, 26);
//! assert_eq!(db.symbol("iron")?, "Fe");
//!
//! // Absorption edges
//! assert_eq!(db.xray_edge("Fe", "K")?.energy, 7112.0);
//!
//! // Mass attenuation, single energy
//! let mu = db.mu_elam_at("Fe", 10_000.0, CrossSectionKind::Total)?;
//! assert!(mu > 0.0);
//!
//! // Compounds, by formula and density
//! let mu_water = db.material_mu_at("H2O", 1.0, 10_000.0, CrossSectionKind::Total)?;
//! assert!((mu_water - 5.33).abs() < 0.05);
//! # Ok::<(), xraydb::XrayDbError>(())
//! ```
//!
//! # Scalar and batch forms
//!
//! Every array-valued calculation has an `_at` scalar counterpart that allocates
//! nothing. Prefer the batch form when evaluating many energies at once — it resolves
//! the element and its table rows a single time for the whole slice.
//!
//! ```
//! # use xraydb::{CrossSectionKind, XrayDb};
//! # let db = XrayDb::try_new()?;
//! let energies: Vec<f64> = (0..100).map(|i| 5_000.0 + 100.0 * i as f64).collect();
//! let mu = db.mu_elam("Fe", &energies, CrossSectionKind::Total)?;
//! assert_eq!(mu.len(), energies.len());
//! # Ok::<(), xraydb::XrayDbError>(())
//! ```
//!
//! # Energy clamping
//!
//! Following upstream XrayDB, energies outside a table's range are **clamped** to its
//! endpoints rather than rejected. Use [`XrayDb::elam_energy_range`] and
//! [`XrayDb::chantler_energy_range`] to check before querying.
//!
//! # Features
//!
//! * `optics` — crystal Darwin widths and mirror/multilayer reflectivity
//!   ([`optics`]). Adds a dependency on `num-complex`.
//!
//! # Attribution
//!
//! X-ray data from the [XrayDB](https://github.com/xraypy/XrayDB) project by
//! Matt Newville et al., released into the public domain (CC0 1.0).

#![doc = ""]
#![doc = "# README"]
#![doc = ""]
#![doc = include_str!("../README.md")]
#![forbid(unsafe_code)]
#![warn(missing_docs)]
#![deny(clippy::unwrap_used, clippy::expect_used)]
#![cfg_attr(test, allow(clippy::unwrap_used, clippy::expect_used))]

pub mod chantler;
pub mod chemparser;
pub mod compton;
pub mod constants;
pub mod core_widths;
pub mod coster_kronig;
pub(crate) mod cubic_spline;
pub mod db;
pub mod elam;
pub mod error;
pub(crate) mod interp;
pub mod ionchamber;
pub mod materials;
pub(crate) mod materials_db;
#[cfg(feature = "optics")]
pub mod optics;
pub(crate) mod spline;
pub mod transitions;
pub mod waasmaier;

pub use chantler::ChantlerKind;
pub use chemparser::{Composition, chemparse};
pub use compton::ComptonEnergies;
pub use db::XrayDb;
pub use elam::{CrossSectionKind, ELAM_ENERGY_MAX, ELAM_ENERGY_MIN};
pub use error::{Result, XrayDbError};
pub use ionchamber::IonChamberFluxes;
pub use materials::RefractiveIndex;
pub use materials_db::Material;
#[cfg(feature = "optics")]
pub use optics::{
    CoatedMirrorParams, DarwinParams, DarwinWidth, MirrorParams, MultilayerParams, Polarization,
};
pub use transitions::{EdgeGuess, LevelLines, XrayEdge, XrayLine};
pub use xraydb_data;
