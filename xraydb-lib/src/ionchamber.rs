//! Ion-chamber flux calculations and named-material lookup.

use crate::constants::ELEMENTARY_CHARGE;
use crate::db::XrayDb;
use crate::elam::CrossSectionKind;
use crate::error::{Result, XrayDbError};
use crate::materials_db::{Material, all_materials, find_material};

/// Photon fluxes derived from an ion-chamber current, in photons/s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IonChamberFluxes {
    /// Flux incident on the chamber.
    pub incident: f64,
    /// Flux leaving the chamber.
    pub transmitted: f64,
    /// Flux absorbed photoelectrically.
    pub photo: f64,
    /// Flux lost to incoherent (Compton) scattering.
    pub incoherent: f64,
    /// Flux lost to coherent (Rayleigh) scattering.
    pub coherent: f64,
}

/// Fallback ionization potential in eV per ion pair, used when a gas has no
/// tabulated value. Close to the value for air and most common counting gases.
pub const DEFAULT_ION_POTENTIAL_EV: f64 = 32.0;

/// Common aliases for gases whose chamber name differs from the materials-table name.
const GAS_ALIASES: &[(&str, &str)] = &[
    ("n2", "nitrogen"),
    ("o2", "oxygen"),
    ("ar", "argon"),
    ("he", "helium"),
    ("kr", "krypton"),
    ("xe", "xenon"),
    ("ne", "neon"),
    ("co2", "carbon dioxide"),
    ("ch4", "methane"),
];

fn resolve_gas_name(name: &str) -> &str {
    let lower = name.to_ascii_lowercase();
    GAS_ALIASES
        .iter()
        .find(|(alias, _)| *alias == lower)
        .map_or(name, |&(_, canonical)| canonical)
}

impl XrayDb {
    /// Look up a material by name (case-insensitive) or by exact formula.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let kapton = db.find_material("kapton").expect("kapton is in the table");
    /// assert_eq!(kapton.formula, "C22H10N2O5");
    /// assert_eq!(kapton.density, 1.42);
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn find_material(&self, name: &str) -> Option<Material> {
        find_material(name)
    }

    /// Every material in the embedded materials table.
    pub fn materials(&self) -> Vec<Material> {
        all_materials().collect()
    }

    /// Linear attenuation coefficient for a named material or a bare formula.
    ///
    /// If `name` is in the materials table, its tabulated density is used unless
    /// `density` overrides it. Otherwise `name` is treated as a formula and `density`
    /// is required.
    pub fn material_mu_named(
        &self,
        name: &str,
        energies: &[f64],
        kind: CrossSectionKind,
        density: Option<f64>,
    ) -> Result<Vec<f64>> {
        let (formula, dens) = self.resolve_named_material(name, density)?;
        self.material_mu(formula, dens, energies, kind)
    }

    /// Single-energy form of [`XrayDb::material_mu_named`].
    pub fn material_mu_named_at(
        &self,
        name: &str,
        energy: f64,
        kind: CrossSectionKind,
        density: Option<f64>,
    ) -> Result<f64> {
        let (formula, dens) = self.resolve_named_material(name, density)?;
        self.material_mu_at(formula, dens, energy, kind)
    }

    fn resolve_named_material<'a>(
        &self,
        name: &'a str,
        density: Option<f64>,
    ) -> Result<(&'a str, f64)> {
        match find_material(name) {
            Some(m) => Ok((m.formula, density.unwrap_or(m.density))),
            None => {
                let d = density.ok_or_else(|| {
                    XrayDbError::invalid_input(format!(
                        "unknown material '{name}': not in the materials table, \
                         so an explicit density is required"
                    ))
                })?;
                Ok((name, d))
            }
        }
    }

    /// Convert a measured ion-chamber current into photon fluxes.
    ///
    /// * `gases` — mixture as `(name, fraction)`; fractions are normalized.
    ///   Names may be materials-table names (`"nitrogen"`) or formulas (`"N2"`).
    /// * `volts` — measured voltage from the current amplifier
    /// * `length_cm` — active length of the chamber
    /// * `energy` — X-ray energy in eV
    /// * `sensitivity` — amplifier sensitivity in A/V
    /// * `with_compton` — include the Compton-electron energy contribution
    /// * `both_carriers` — count both electron and ion carriers (true for most chambers)
    #[allow(clippy::too_many_arguments)]
    pub fn ionchamber_fluxes(
        &self,
        gases: &[(&str, f64)],
        volts: f64,
        length_cm: f64,
        energy: f64,
        sensitivity: f64,
        with_compton: bool,
        both_carriers: bool,
    ) -> Result<IonChamberFluxes> {
        if gases.is_empty() {
            return Err(XrayDbError::invalid_input("gas mixture must not be empty"));
        }
        if !length_cm.is_finite() || length_cm <= 0.0 {
            return Err(XrayDbError::invalid_input(format!(
                "chamber length must be finite and > 0 cm, got {length_cm}"
            )));
        }

        let ncarriers: f64 = if both_carriers { 2.0 } else { 1.0 };

        let gas_total: f64 = gases.iter().map(|(_, f)| f).sum();
        if !gas_total.is_finite() || gas_total <= 0.0 {
            return Err(XrayDbError::invalid_input("gas fractions must sum to > 0"));
        }

        let energy_compton = if with_compton {
            self.compton_energies(energy).electron_mean
        } else {
            0.0
        };

        let mut mu_photo = 0.0;
        let mut mu_incoh = 0.0;
        let mut mu_total = 0.0;
        let mut mu_coh = 0.0;
        let mut ion_pot = 0.0;

        for &(gas_name, frac) in gases {
            let weight = frac / gas_total;
            let canonical = resolve_gas_name(gas_name);

            let ip = self
                .ionization_potential(gas_name)
                .or_else(|_| self.ionization_potential(canonical))
                .unwrap_or(DEFAULT_ION_POTENTIAL_EV);

            // Resolve and parse the material once, then evaluate all four
            // cross-sections against the same mass fractions.
            let (formula, density) = self.resolve_named_material(canonical, None)?;
            let fractions = self.mass_fractions(formula)?;

            mu_photo += weight
                * self.mu_from_fractions_at(
                    &fractions,
                    density,
                    energy,
                    CrossSectionKind::Photo,
                )?;
            mu_total += weight
                * self.mu_from_fractions_at(
                    &fractions,
                    density,
                    energy,
                    CrossSectionKind::Total,
                )?;
            mu_incoh += weight
                * self.mu_from_fractions_at(
                    &fractions,
                    density,
                    energy,
                    CrossSectionKind::Incoherent,
                )?;
            mu_coh += weight
                * self.mu_from_fractions_at(
                    &fractions,
                    density,
                    energy,
                    CrossSectionKind::Coherent,
                )?;
            ion_pot += ip * weight;
        }

        let atten_total = 1.0 - (-length_cm * mu_total).exp();
        let share = |mu: f64| {
            if mu_total > 0.0 {
                atten_total * mu / mu_total
            } else {
                0.0
            }
        };
        let atten_photo = share(mu_photo);
        let atten_incoh = share(mu_incoh);
        let atten_coh = share(mu_coh);

        let absorbed_energy = ncarriers * (energy * atten_photo + energy_compton * atten_incoh);

        let flux_in = if absorbed_energy > 0.0 {
            volts * sensitivity * ion_pot / (ELEMENTARY_CHARGE * absorbed_energy)
        } else {
            0.0
        };

        Ok(IonChamberFluxes {
            incident: flux_in,
            transmitted: flux_in * (1.0 - atten_total),
            photo: flux_in * atten_photo,
            incoherent: flux_in * atten_incoh,
            coherent: flux_in * atten_coh,
        })
    }
}
