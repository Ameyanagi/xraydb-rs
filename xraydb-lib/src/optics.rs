//! X-ray optics: crystal Darwin widths, mirror and multilayer reflectivity.
//!
//! Requires the `optics` feature.
//!
//! Every entry point takes a parameter struct with a `Default` impl, so call sites read
//! as named arguments rather than a row of anonymous floats and bools:
//!
//! ```
//! # #[cfg(feature = "optics")] {
//! use xraydb::{XrayDb, optics::DarwinParams};
//!
//! let db = XrayDb::try_new()?;
//! let dw = db.darwin_width(DarwinParams {
//!     energy: 10_000.0,
//!     crystal: "Si",
//!     hkl: (1, 1, 1),
//!     ..Default::default()
//! })?;
//! assert!(dw.is_some());
//! # }
//! # Ok::<(), xraydb::XrayDbError>(())
//! ```

use num_complex::Complex64;
use std::f64::consts::{PI, SQRT_2};

use crate::constants::{PLANCK_HC_ANGSTROM, R_ELECTRON_ANG};
use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

/// Polarization state for optics calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Polarization {
    /// Perpendicular to the scattering plane (σ).
    #[default]
    S,
    /// Parallel to the scattering plane (π).
    P,
    /// Average of S and P. Not supported for multilayers.
    Unpolarized,
}

/// Lattice constants in Å for the diamond-structure crystals supported here.
const CRYSTALS: &[(&str, f64)] = &[
    ("si", 5.4309),
    ("ge", 5.6578),
    ("c", 3.567),
    ("diamond", 3.567),
];

/// Lattice constant in Å for a supported crystal, by case-insensitive name.
///
/// ```
/// # #[cfg(feature = "optics")] {
/// assert_eq!(xraydb::optics::crystal_lattice_constant("Si"), Some(5.4309));
/// assert_eq!(xraydb::optics::crystal_lattice_constant("NaCl"), None);
/// # }
/// ```
pub fn crystal_lattice_constant(name: &str) -> Option<f64> {
    CRYSTALS
        .iter()
        .find(|(n, _)| n.eq_ignore_ascii_case(name))
        .map(|&(_, a)| a)
}

/// Names of the crystals [`crystal_lattice_constant`] recognises.
pub fn supported_crystals() -> Vec<&'static str> {
    CRYSTALS.iter().map(|&(n, _)| n).collect()
}

/// Upper bound on the sampled rocking-curve grid, so a pathological input cannot
/// request an unbounded allocation.
const MAX_ZETA_POINTS: usize = 100_000;

/// Arguments to [`XrayDb::darwin_width`].
#[derive(Debug, Clone, Copy)]
pub struct DarwinParams<'a> {
    /// X-ray energy in eV.
    pub energy: f64,
    /// Crystal name: `"Si"`, `"Ge"`, or `"C"`/`"diamond"`.
    pub crystal: &'a str,
    /// Miller indices, all even (with h+k+l divisible by 4) or all odd.
    pub hkl: (i32, i32, i32),
    /// Lattice constant in Å; `None` uses the built-in value for `crystal`.
    pub lattice_constant: Option<f64>,
    /// Polarization state.
    pub polarization: Polarization,
    /// Skip the f1 dispersion correction.
    pub ignore_f1: bool,
    /// Skip f2 absorption.
    pub ignore_f2: bool,
    /// Reflection order.
    pub order: i32,
}

impl Default for DarwinParams<'_> {
    fn default() -> Self {
        DarwinParams {
            energy: 10_000.0,
            crystal: "Si",
            hkl: (1, 1, 1),
            lattice_constant: None,
            polarization: Polarization::S,
            ignore_f1: false,
            ignore_f2: false,
            order: 1,
        }
    }
}

/// Arguments to [`XrayDb::mirror_reflectivity`].
#[derive(Debug, Clone, Copy)]
pub struct MirrorParams<'a> {
    /// Mirror material formula, e.g. `"Pt"` or `"SiO2"`.
    pub formula: &'a str,
    /// Grazing angles in radians.
    pub theta: &'a [f64],
    /// X-ray energy in eV.
    pub energy: f64,
    /// Density in g/cm³.
    pub density: f64,
    /// RMS surface roughness in Å; 0 for an ideal surface.
    pub roughness: f64,
    /// Polarization state.
    pub polarization: Polarization,
}

impl Default for MirrorParams<'_> {
    fn default() -> Self {
        MirrorParams {
            formula: "Si",
            theta: &[],
            energy: 10_000.0,
            density: 2.329,
            roughness: 0.0,
            polarization: Polarization::S,
        }
    }
}

/// Arguments to [`XrayDb::multilayer_reflectivity`].
#[derive(Debug, Clone, Copy)]
pub struct MultilayerParams<'a> {
    /// Layer formulas, surface first.
    pub stackup: &'a [&'a str],
    /// Layer thicknesses in Å, matching `stackup`.
    pub thickness: &'a [f64],
    /// Layer densities in g/cm³, matching `stackup`.
    pub density: &'a [f64],
    /// Substrate formula.
    pub substrate: &'a str,
    /// Substrate density in g/cm³.
    pub substrate_density: f64,
    /// Grazing angles in radians.
    pub theta: &'a [f64],
    /// X-ray energy in eV.
    pub energy: f64,
    /// How many times to repeat `stackup`.
    pub n_periods: usize,
    /// Substrate–layer interface roughness in Å.
    pub substrate_roughness: f64,
    /// Air–surface interface roughness in Å.
    pub surface_roughness: f64,
    /// Polarization state; [`Polarization::Unpolarized`] is rejected.
    pub polarization: Polarization,
}

impl Default for MultilayerParams<'_> {
    fn default() -> Self {
        MultilayerParams {
            stackup: &[],
            thickness: &[],
            density: &[],
            substrate: "Si",
            substrate_density: 2.329,
            theta: &[],
            energy: 10_000.0,
            n_periods: 1,
            substrate_roughness: 0.0,
            surface_roughness: 0.0,
            polarization: Polarization::S,
        }
    }
}

/// Arguments to [`XrayDb::coated_reflectivity`].
#[derive(Debug, Clone, Copy)]
pub struct CoatedMirrorParams<'a> {
    /// Coating formula.
    pub coating: &'a str,
    /// Coating thickness in Å.
    pub coating_thickness: f64,
    /// Coating density in g/cm³.
    pub coating_density: f64,
    /// Substrate formula.
    pub substrate: &'a str,
    /// Substrate density in g/cm³.
    pub substrate_density: f64,
    /// Grazing angles in radians.
    pub theta: &'a [f64],
    /// X-ray energy in eV.
    pub energy: f64,
    /// Air–coating roughness in Å.
    pub surface_roughness: f64,
    /// Substrate–coating roughness in Å.
    pub substrate_roughness: f64,
    /// Optional binder layer: `(formula, thickness_Å, density)`.
    pub binder: Option<(&'a str, f64, f64)>,
    /// Polarization state.
    pub polarization: Polarization,
}

impl Default for CoatedMirrorParams<'_> {
    fn default() -> Self {
        CoatedMirrorParams {
            coating: "Pt",
            coating_thickness: 500.0,
            coating_density: 21.45,
            substrate: "Si",
            substrate_density: 2.329,
            theta: &[],
            energy: 10_000.0,
            surface_roughness: 0.0,
            substrate_roughness: 0.0,
            binder: None,
            polarization: Polarization::S,
        }
    }
}

/// Result of a Darwin-width calculation for a crystal reflection.
#[derive(Debug, Clone)]
pub struct DarwinWidth {
    /// Bragg angle in radians.
    pub theta: f64,
    /// Angular offset of the reflection peak, radians.
    pub theta_offset: f64,
    /// Intrinsic angular width, radians.
    pub theta_width: f64,
    /// FWHM angular width, radians.
    pub theta_fwhm: f64,
    /// Rocking-curve FWHM in angle, radians.
    pub rocking_theta_fwhm: f64,
    /// Intrinsic energy width, eV.
    pub energy_width: f64,
    /// FWHM energy width, eV.
    pub energy_fwhm: f64,
    /// Rocking-curve FWHM in energy, eV.
    pub rocking_energy_fwhm: f64,
    /// Dimensionless deviation parameter grid.
    pub zeta: Vec<f64>,
    /// Angular deviation grid, radians.
    pub dtheta: Vec<f64>,
    /// Energy deviation grid, eV.
    pub denergy: Vec<f64>,
    /// Single-crystal reflectivity on the `zeta` grid.
    pub intensity: Vec<f64>,
    /// Two-crystal rocking curve (self-convolution of `intensity`).
    pub rocking_curve: Vec<f64>,
}

/// Real number as a complex value.
#[inline]
fn c(re: f64) -> Complex64 {
    Complex64::new(re, 0.0)
}

/// Discrete convolution with 'same' output size (centered).
fn convolve_same(a: &[f64], b: &[f64]) -> Vec<f64> {
    let na = a.len();
    let nb = b.len();
    if na == 0 || nb == 0 {
        return Vec::new();
    }
    let full_len = na + nb - 1;
    let mut full = vec![0.0; full_len];
    for (i, &ai) in a.iter().enumerate() {
        if ai == 0.0 {
            continue;
        }
        for (j, &bj) in b.iter().enumerate() {
            full[i + j] += ai * bj;
        }
    }
    let start = (full_len - na) / 2;
    full[start..start + na].to_vec()
}

impl XrayDb {
    /// Darwin width and rocking curve for a crystal reflection.
    ///
    /// Supports diamond-structure crystals (Si, Ge, C). Returns `Ok(None)` when the
    /// Bragg condition cannot be satisfied at this energy, or when the reflection is
    /// too weak to define a rocking-curve grid.
    pub fn darwin_width(&self, params: DarwinParams<'_>) -> Result<Option<DarwinWidth>> {
        let DarwinParams {
            energy,
            crystal,
            hkl,
            lattice_constant,
            polarization,
            ignore_f1,
            ignore_f2,
            order,
        } = params;

        if !energy.is_finite() || energy <= 0.0 {
            return Err(XrayDbError::invalid_input(format!(
                "energy must be finite and > 0 eV, got {energy}"
            )));
        }
        if order < 1 {
            return Err(XrayDbError::invalid_input(format!(
                "reflection order must be >= 1, got {order}"
            )));
        }

        let (h, k, l) = hkl;
        if (h, k, l) == (0, 0, 0) {
            return Err(XrayDbError::invalid_input("hkl must not be (0, 0, 0)"));
        }

        // Structure factor magnitude for the diamond lattice.
        let hkl_sum = h + k + l;
        let all_even = h % 2 == 0 && k % 2 == 0 && l % 2 == 0;
        let all_odd = h % 2 != 0 && k % 2 != 0 && l % 2 != 0;
        let eqr: f64 = if all_even && hkl_sum % 4 == 0 {
            8.0
        } else if all_odd {
            4.0 * SQRT_2
        } else {
            return Err(XrayDbError::invalid_input(
                "hkl must be all even with h+k+l divisible by 4, or all odd",
            ));
        };

        let a = match lattice_constant {
            Some(a) => a,
            None => crystal_lattice_constant(crystal).ok_or_else(|| {
                XrayDbError::invalid_input(format!(
                    "unsupported crystal '{crystal}'; supported: {}",
                    supported_crystals().join(", ")
                ))
            })?,
        };
        if !a.is_finite() || a <= 0.0 {
            return Err(XrayDbError::invalid_input(format!(
                "lattice constant must be finite and > 0 Å, got {a}"
            )));
        }

        let dspace = a / f64::from(h * h + k * k + l * l).sqrt();
        let lambd = PLANCK_HC_ANGSTROM / energy;

        if lambd > 2.0 * dspace {
            return Ok(None);
        }

        let theta = (lambd / (2.0 * dspace)).asin();
        let q = 0.5 / dspace;

        let f1_val = if ignore_f1 {
            0.0
        } else {
            self.f1_chantler_at(crystal, energy)?
        };
        let f2_val = if ignore_f2 {
            0.0
        } else {
            self.f2_chantler_at(crystal, energy)?
        };

        let mf = f64::from(order);
        let gscale = 2.0 * dspace * dspace * R_ELECTRON_ANG / (mf * a.powi(3));

        let eqr = match polarization {
            Polarization::Unpolarized => eqr * (1.0 + (2.0 * theta).cos().abs()) / 2.0,
            Polarization::P => eqr * (2.0 * theta).cos().abs(),
            Polarization::S => eqr,
        };

        let f0_0 = self.f0_at(crystal, 0.0)?;
        let f0_q = self.f0_at(crystal, q)?;

        let f_anom = Complex64::new(f1_val, -f2_val);
        // Forward scattering (q = 0), always unpolarized.
        let g0 = c(8.0 * gscale) * (c(f0_0) + f_anom);
        // Bragg scattering at q.
        let g = c(eqr * gscale) * (c(f0_q) + f_anom);

        let total = (c(2.0) * g / c(mf * PI)).norm();
        let fwhm = total * 3.0 / (2.0 * SQRT_2);

        let zeta_offset = g0.re / PI;
        let theta_offset = theta.tan() * zeta_offset;

        // The grid spans [-2.5ζ₀, 4.5ζ₀]; a non-positive or non-finite ζ₀ makes it
        // empty or inverted, which previously produced a `Some` full of empty vectors.
        let zeta_step = 0.01 * total;
        if !zeta_step.is_finite()
            || zeta_step <= 0.0
            || !zeta_offset.is_finite()
            || zeta_offset <= 0.0
        {
            return Ok(None);
        }

        let zeta_start = -2.5 * zeta_offset;
        let zeta_end = 4.5 * zeta_offset;
        let raw_points = ((zeta_end - zeta_start) / zeta_step).ceil();
        if !raw_points.is_finite() || raw_points < 1.0 {
            return Ok(None);
        }
        let n_points = (raw_points as usize).min(MAX_ZETA_POINTS);

        let zeta: Vec<f64> = (0..n_points)
            .map(|i| zeta_start + i as f64 * zeta_step)
            .collect();

        let one = c(1.0);
        let intensity: Vec<f64> = zeta
            .iter()
            .map(|&z| {
                let xc = (c(mf * PI * z) - g0) / g;
                let r = if xc.re > 1.0 {
                    xc - (xc * xc - one).sqrt()
                } else if xc.re < -1.0 {
                    xc + (xc * xc - one).sqrt()
                } else {
                    xc - Complex64::i() * (one - xc * xc).sqrt()
                };
                (r * r.conj()).re
            })
            .collect();

        let denergy: Vec<f64> = zeta.iter().map(|&z| -z * energy).collect();
        let dtheta: Vec<f64> = zeta.iter().map(|&z| z * theta.tan()).collect();

        let intensity_sum: f64 = intensity.iter().sum();
        let rocking_curve = if intensity_sum > 0.0 {
            convolve_same(&intensity, &intensity)
                .iter()
                .map(|&v| v / intensity_sum)
                .collect()
        } else {
            vec![0.0; intensity.len()]
        };

        let rc_max = rocking_curve
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let half_max = rc_max / 2.0;

        let first = rocking_curve.iter().position(|&v| v >= half_max);
        let last = rocking_curve.iter().rposition(|&v| v >= half_max);
        let (re_fwhm, rt_fwhm) = match (first, last) {
            (Some(f), Some(l)) if l > f => (
                (denergy[l] - denergy[f]).abs(),
                (dtheta[l] - dtheta[f]).abs(),
            ),
            _ => (0.0, 0.0),
        };

        Ok(Some(DarwinWidth {
            theta,
            theta_offset,
            theta_width: total * theta.tan(),
            theta_fwhm: fwhm * theta.tan(),
            rocking_theta_fwhm: rt_fwhm,
            energy_width: total * energy,
            energy_fwhm: fwhm * energy,
            rocking_energy_fwhm: re_fwhm,
            zeta,
            dtheta,
            denergy,
            intensity,
            rocking_curve,
        }))
    }

    /// Reflectivity of a thick, single-material mirror.
    pub fn mirror_reflectivity(&self, params: MirrorParams<'_>) -> Result<Vec<f64>> {
        let MirrorParams {
            formula,
            theta,
            energy,
            density,
            roughness,
            polarization,
        } = params;

        let n_index = self.xray_delta_beta(formula, density, energy)?;
        let n = Complex64::new(1.0 - n_index.delta, -n_index.beta);
        let qf = 2.0 * PI * energy / PLANCK_HC_ANGSTROM;

        Ok(theta
            .iter()
            .map(|&th| {
                let kiz = c(qf * th.sin());
                let cos_th = th.cos();
                let mut ktz = (n * n - c(cos_th * cos_th)).sqrt() * c(qf);

                if polarization == Polarization::P {
                    ktz /= n;
                }

                let mut r = (kiz - ktz) / (kiz + ktz);
                if roughness > 1e-12 {
                    r *= (c(-2.0 * roughness * roughness) * kiz * ktz).exp();
                }
                (r * r.conj()).re
            })
            .collect())
    }

    /// Multilayer reflectivity via the Parratt recursion.
    pub fn multilayer_reflectivity(&self, params: MultilayerParams<'_>) -> Result<Vec<f64>> {
        let MultilayerParams {
            stackup,
            thickness,
            density,
            substrate,
            substrate_density,
            theta,
            energy,
            n_periods,
            substrate_roughness,
            surface_roughness,
            polarization,
        } = params;

        if stackup.is_empty() {
            return Err(XrayDbError::invalid_input("stackup must not be empty"));
        }
        if n_periods == 0 {
            return Err(XrayDbError::invalid_input("n_periods must be >= 1"));
        }
        if stackup.len() != thickness.len() {
            return Err(XrayDbError::invalid_input(format!(
                "stackup ({}) and thickness ({}) lengths must match",
                stackup.len(),
                thickness.len()
            )));
        }
        if stackup.len() != density.len() {
            return Err(XrayDbError::invalid_input(format!(
                "stackup ({}) and density ({}) lengths must match",
                stackup.len(),
                density.len()
            )));
        }
        if polarization == Polarization::Unpolarized {
            return Err(XrayDbError::invalid_input(
                "multilayer reflectivity requires S or P polarization",
            ));
        }

        let k0 = 2.0 * PI * energy / PLANCK_HC_ANGSTROM;

        // Parratt convention: n = 1 - delta + i*beta.
        let mut n_vals = Vec::with_capacity(stackup.len());
        for (i, layer) in stackup.iter().enumerate() {
            let idx = self.xray_delta_beta(layer, density[i], energy)?;
            n_vals.push(Complex64::new(1.0 - idx.delta, idx.beta));
        }

        let total_layers = stackup.len() * n_periods;
        let mut t_all = Vec::with_capacity(total_layers);
        let mut n_all = Vec::with_capacity(total_layers);
        for _ in 0..n_periods {
            t_all.extend_from_slice(thickness);
            n_all.extend_from_slice(&n_vals);
        }

        let sub_idx = self.xray_delta_beta(substrate, substrate_density, energy)?;
        let n_sub = Complex64::new(1.0 - sub_idx.delta, sub_idx.beta);

        // `total_layers >= 1` is guaranteed by the checks above.
        let last = total_layers - 1;
        let one = c(1.0);
        let two_i = Complex64::new(0.0, 2.0);

        let mut result = Vec::with_capacity(theta.len());
        for &th in theta {
            let cos_th = th.cos();
            let cos2 = c(cos_th * cos_th);
            let kiz = c(k0 * th.sin());

            let kz: Vec<Complex64> = n_all
                .iter()
                .map(|&ni| (ni * ni - cos2).sqrt() * c(k0))
                .collect();
            let kz_sub = (n_sub * n_sub - cos2).sqrt() * c(k0);

            let mut r_amp = match polarization {
                Polarization::P => {
                    let a = kz[last] / n_all[last] * n_sub;
                    let b = kz_sub / n_sub * n_all[last];
                    (a - b) / (a + b)
                }
                _ => (kz[last] - kz_sub) / (kz[last] + kz_sub),
            };

            if substrate_roughness >= 1e-12 {
                r_amp *=
                    (c(-2.0 * substrate_roughness * substrate_roughness) * kz[last] * kz_sub).exp();
            }

            for i in (0..last).rev() {
                let fresnel_r = match polarization {
                    Polarization::P => {
                        let a = kz[i] / n_all[i] * n_all[i + 1];
                        let b = kz[i + 1] / n_all[i + 1] * n_all[i];
                        (a - b) / (a + b)
                    }
                    _ => (kz[i] - kz[i + 1]) / (kz[i] + kz[i + 1]),
                };
                let p2 = (two_i * c(t_all[i + 1]) * kz[i + 1]).exp();
                r_amp = (fresnel_r + r_amp * p2) / (one + fresnel_r * r_amp * p2);
            }

            let fresnel_r = match polarization {
                Polarization::P => (kiz - kz[0] / n_all[0]) / (kiz + kz[0] / n_all[0]),
                _ => (kiz - kz[0]) / (kiz + kz[0]),
            };
            let p2 = (two_i * c(t_all[0]) * kz[0]).exp();
            r_amp = (fresnel_r + r_amp * p2) / (one + fresnel_r * r_amp * p2);

            if surface_roughness >= 1e-12 {
                r_amp *= (c(-2.0 * surface_roughness * surface_roughness) * kiz * kz[0]).exp();
            }

            result.push((r_amp * r_amp.conj()).re);
        }

        Ok(result)
    }

    /// Reflectivity of a coated mirror, with an optional binder layer.
    ///
    /// A thin convenience wrapper over [`XrayDb::multilayer_reflectivity`].
    pub fn coated_reflectivity(&self, params: CoatedMirrorParams<'_>) -> Result<Vec<f64>> {
        let CoatedMirrorParams {
            coating,
            coating_thickness,
            coating_density,
            substrate,
            substrate_density,
            theta,
            energy,
            surface_roughness,
            substrate_roughness,
            binder,
            polarization,
        } = params;

        let (stackup, thickness, densities): (Vec<&str>, Vec<f64>, Vec<f64>) = match binder {
            Some((binder_mat, binder_thick, binder_dens)) => (
                vec![coating, binder_mat],
                vec![coating_thickness, binder_thick],
                vec![coating_density, binder_dens],
            ),
            None => (
                vec![coating],
                vec![coating_thickness],
                vec![coating_density],
            ),
        };

        self.multilayer_reflectivity(MultilayerParams {
            stackup: &stackup,
            thickness: &thickness,
            density: &densities,
            substrate,
            substrate_density,
            theta,
            energy,
            n_periods: 1,
            substrate_roughness,
            surface_roughness,
            polarization,
        })
    }
}
