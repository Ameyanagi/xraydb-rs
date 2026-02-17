//! X-ray optics calculations: mirror/multilayer reflectivity, Darwin width.
//!
//! Requires the `optics` feature.

use num_complex::Complex64;
use std::f64::consts::PI;

use crate::constants::{PLANCK_HC_ANGSTROM, R_ELECTRON_ANG};
use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

/// Polarization state for X-ray optics calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarization {
    S,
    P,
    Unpolarized,
}

/// Result of Darwin width calculation for crystal diffraction.
#[derive(Debug, Clone)]
pub struct DarwinWidth {
    /// Bragg angle (radians)
    pub theta: f64,
    /// Angular offset of reflection peak (radians)
    pub theta_offset: f64,
    /// Intrinsic angular width (radians)
    pub theta_width: f64,
    /// FWHM angular width (radians)
    pub theta_fwhm: f64,
    /// Rocking curve FWHM in angle (radians)
    pub rocking_theta_fwhm: f64,
    /// Intrinsic energy width (eV)
    pub energy_width: f64,
    /// FWHM energy width (eV)
    pub energy_fwhm: f64,
    /// Rocking curve FWHM in energy (eV)
    pub rocking_energy_fwhm: f64,
    /// Zeta parameter array
    pub zeta: Vec<f64>,
    /// Angular deviation array (radians)
    pub dtheta: Vec<f64>,
    /// Energy deviation array (eV)
    pub denergy: Vec<f64>,
    /// Reflectivity intensity array
    pub intensity: Vec<f64>,
    /// Rocking curve (self-convolution of intensity)
    pub rocking_curve: Vec<f64>,
}

/// Convert f64 to Complex64 (real part only).
#[inline]
fn c(re: f64) -> Complex64 {
    Complex64::new(re, 0.0)
}

/// Discrete convolution with 'same' output size (centered).
fn convolve_same(a: &[f64], b: &[f64]) -> Vec<f64> {
    let na = a.len();
    let nb = b.len();
    let full_len = na + nb - 1;
    let mut full = vec![0.0; full_len];
    for i in 0..na {
        let ai = a[i];
        if ai == 0.0 {
            continue;
        }
        for j in 0..nb {
            full[i + j] += ai * b[j];
        }
    }
    let start = (full_len - na) / 2;
    full[start..start + na].to_vec()
}

impl XrayDb {
    /// Calculate Darwin width for a crystal reflection.
    ///
    /// Supports diamond-structure crystals: Si, Ge, C (diamond).
    /// Returns `None` if the Bragg condition cannot be satisfied.
    ///
    /// # Arguments
    /// * `energy` - X-ray energy in eV
    /// * `crystal` - Crystal name: "Si", "Ge", or "C"
    /// * `hkl` - Miller indices (h, k, l)
    /// * `a` - Lattice constant in Å (None for built-in value)
    /// * `polarization` - S, P, or Unpolarized
    /// * `ignore_f1` - Ignore f1 dispersion correction
    /// * `ignore_f2` - Ignore f2 absorption
    /// * `m` - Reflection order
    pub fn darwin_width(
        &self,
        energy: f64,
        crystal: &str,
        hkl: (i32, i32, i32),
        a: Option<f64>,
        polarization: Polarization,
        ignore_f1: bool,
        ignore_f2: bool,
        m: i32,
    ) -> Result<Option<DarwinWidth>> {
        let (h, k, l) = hkl;
        let hkl_sum = h + k + l;

        // Structure factor for diamond-structure crystals
        let eqr: f64 = if hkl_sum % 4 == 0 && h % 2 == 0 && k % 2 == 0 && l % 2 == 0 {
            8.0
        } else if h % 2 != 0 && k % 2 != 0 && l % 2 != 0 {
            4.0 * std::f64::consts::SQRT_2
        } else {
            return Err(XrayDbError::DataError(
                "hkl must all be even (sum divisible by 4) or all odd".to_string(),
            ));
        };

        let lattice = match crystal.to_lowercase().as_str() {
            "si" => 5.4309,
            "ge" => 5.6578,
            "c" | "diamond" => 3.567,
            _ => {
                return Err(XrayDbError::DataError(format!(
                    "unsupported crystal '{crystal}', use Si, Ge, or C"
                )));
            }
        };
        let a = a.unwrap_or(lattice);

        let dspace = a / ((h * h + k * k + l * l) as f64).sqrt();
        let lambd = PLANCK_HC_ANGSTROM / energy;

        // Check Bragg condition
        if lambd > 2.0 * dspace {
            return Ok(None);
        }

        let theta = (lambd / (2.0 * dspace)).asin();
        let q = 0.5 / dspace;

        let mut f1_val = 0.0;
        let mut f2_val = 0.0;
        if !ignore_f1 {
            f1_val = self.f1_chantler(crystal, &[energy])?[0];
        }
        if !ignore_f2 {
            f2_val = self.f2_chantler(crystal, &[energy])?[0];
        }

        let mf = m as f64;
        let gscale = 2.0 * dspace * dspace * R_ELECTRON_ANG / (mf * a.powi(3));

        // Apply polarization factor
        let eqr = match polarization {
            Polarization::Unpolarized => eqr * (1.0 + (2.0 * theta).cos().abs()) / 2.0,
            Polarization::P => eqr * (2.0 * theta).cos().abs(),
            Polarization::S => eqr,
        };

        let f0_0 = self.f0(crystal, &[0.0])?[0];
        let f0_q = self.f0(crystal, &[q])?[0];

        let f_anom = Complex64::new(f1_val, -f2_val);
        // g0: forward scattering (q=0), always unpolarized
        let g0 = c(8.0 * gscale) * (c(f0_0) + f_anom);
        // g: Bragg scattering at q
        let g = c(eqr * gscale) * (c(f0_q) + f_anom);

        let total = (c(2.0) * g / c(mf * PI)).norm();
        let fwhm = total * 3.0 / (2.0 * std::f64::consts::SQRT_2);

        let zeta_offset = g0.re / PI;
        let theta_offset = theta.tan() * zeta_offset;

        // Build zeta array centered at zeta_offset
        let zeta_step = 0.01 * total;
        if zeta_step <= 0.0 {
            return Ok(None);
        }
        let zeta_start = -2.5 * zeta_offset;
        let zeta_end = 4.5 * zeta_offset;
        let n_points = ((zeta_end - zeta_start) / zeta_step).ceil() as usize;
        let zeta: Vec<f64> = (0..n_points)
            .map(|i| zeta_start + i as f64 * zeta_step)
            .collect();

        // Compute reflectivity at each zeta
        let mut intensity = Vec::with_capacity(zeta.len());
        let one = c(1.0);
        for &z in &zeta {
            let xc = (c(mf * PI * z) - g0) / g;

            let r = if xc.re > 1.0 {
                xc - (xc * xc - one).sqrt()
            } else if xc.re < -1.0 {
                xc + (xc * xc - one).sqrt()
            } else {
                xc - Complex64::i() * (one - xc * xc).sqrt()
            };

            intensity.push((r * r.conj()).re);
        }

        let denergy: Vec<f64> = zeta.iter().map(|&z| -z * energy).collect();
        let dtheta: Vec<f64> = zeta.iter().map(|&z| z * theta.tan()).collect();

        // Rocking curve = self-convolution normalized by intensity sum
        let intensity_sum: f64 = intensity.iter().sum();
        let rocking_curve = if intensity_sum > 0.0 {
            convolve_same(&intensity, &intensity)
                .iter()
                .map(|&v| v / intensity_sum)
                .collect()
        } else {
            vec![0.0; intensity.len()]
        };

        // Find FWHM of rocking curve
        let rc_max = rocking_curve
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let half_max = rc_max / 2.0;

        let big: Vec<usize> = rocking_curve
            .iter()
            .enumerate()
            .filter(|&(_, v)| *v >= half_max)
            .map(|(i, _)| i)
            .collect();

        let (re_fwhm, rt_fwhm) = if big.len() >= 2 {
            let first = big[0];
            let last = big[big.len() - 1];
            (
                (denergy[last] - denergy[first]).abs(),
                (dtheta[last] - dtheta[first]).abs(),
            )
        } else {
            (0.0, 0.0)
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

    /// Mirror reflectivity for a thick, single-layer mirror.
    ///
    /// # Arguments
    /// * `formula` - Chemical formula (e.g., "SiO2", "Pt")
    /// * `theta` - Grazing angles in radians
    /// * `energy` - X-ray energy in eV
    /// * `density` - Density in g/cm³
    /// * `roughness` - RMS surface roughness in Å (0 for ideal surface)
    /// * `polarization` - S or P polarization
    pub fn mirror_reflectivity(
        &self,
        formula: &str,
        theta: &[f64],
        energy: f64,
        density: f64,
        roughness: f64,
        polarization: Polarization,
    ) -> Result<Vec<f64>> {
        let (delta, beta, _) = self.xray_delta_beta(formula, density, energy)?;
        let n = Complex64::new(1.0 - delta, -beta);
        let qf = 2.0 * PI * energy / PLANCK_HC_ANGSTROM;

        let mut result = Vec::with_capacity(theta.len());
        for &th in theta {
            let sin_th = th.sin();
            let cos_th = th.cos();

            let kiz = c(qf * sin_th);
            let mut ktz = (n * n - c(cos_th * cos_th)).sqrt() * c(qf);

            if polarization == Polarization::P {
                ktz = ktz / n;
            }

            let mut r = (kiz - ktz) / (kiz + ktz);

            if roughness > 1e-12 {
                r = r * (c(-2.0 * roughness * roughness) * kiz * ktz).exp();
            }

            result.push((r * r.conj()).re);
        }

        Ok(result)
    }

    /// Multilayer reflectivity using Parratt recursion.
    ///
    /// # Arguments
    /// * `stackup` - Layer materials from surface to substrate (formulas)
    /// * `thickness` - Layer thicknesses in Å (matching stackup)
    /// * `substrate` - Substrate material formula
    /// * `theta` - Grazing angles in radians
    /// * `energy` - X-ray energy in eV
    /// * `n_periods` - Number of times to repeat the stackup
    /// * `density` - Densities in g/cm³ for each layer in stackup
    /// * `substrate_density` - Substrate density in g/cm³
    /// * `substrate_rough` - Substrate-layer interface roughness in Å
    /// * `surface_rough` - Air-surface interface roughness in Å
    /// * `polarization` - S or P polarization
    pub fn multilayer_reflectivity(
        &self,
        stackup: &[&str],
        thickness: &[f64],
        substrate: &str,
        theta: &[f64],
        energy: f64,
        n_periods: usize,
        density: &[f64],
        substrate_density: f64,
        substrate_rough: f64,
        surface_rough: f64,
        polarization: Polarization,
    ) -> Result<Vec<f64>> {
        if stackup.len() != thickness.len() {
            return Err(XrayDbError::DataError(format!(
                "stackup ({}) and thickness ({}) lengths must match",
                stackup.len(),
                thickness.len()
            )));
        }
        if stackup.len() != density.len() {
            return Err(XrayDbError::DataError(format!(
                "stackup ({}) and density ({}) lengths must match",
                stackup.len(),
                density.len()
            )));
        }

        let k0 = 2.0 * PI * energy / PLANCK_HC_ANGSTROM;

        // Compute n for each unique layer (Parratt convention: n = 1 - delta + i*beta)
        let n_layers = stackup.len();
        let mut n_vals = Vec::with_capacity(n_layers);
        for i in 0..n_layers {
            let (delta, beta, _) = self.xray_delta_beta(stackup[i], density[i], energy)?;
            n_vals.push(Complex64::new(1.0 - delta, beta));
        }

        // Repeat for n_periods
        let mut t_all = Vec::with_capacity(n_layers * n_periods);
        let mut n_all = Vec::with_capacity(n_layers * n_periods);
        for _ in 0..n_periods {
            t_all.extend_from_slice(thickness);
            n_all.extend_from_slice(&n_vals);
        }

        // Substrate
        let (delta_sub, beta_sub, _) =
            self.xray_delta_beta(substrate, substrate_density, energy)?;
        let n_sub = Complex64::new(1.0 - delta_sub, beta_sub);

        let total_layers = t_all.len();
        let mut result = Vec::with_capacity(theta.len());
        let one = c(1.0);
        let two_i = Complex64::new(0.0, 2.0);

        for &th in theta {
            let sin_th = th.sin();
            let cos_th = th.cos();
            let cos2 = c(cos_th * cos_th);

            let kiz = c(k0 * sin_th);

            // kz for each layer
            let kz: Vec<Complex64> = n_all
                .iter()
                .map(|&ni| (ni * ni - cos2).sqrt() * c(k0))
                .collect();
            let kz_sub = (n_sub * n_sub - cos2).sqrt() * c(k0);

            let last = total_layers - 1;

            // Parratt recursion starting from substrate
            let mut r_amp = match polarization {
                Polarization::S => (kz[last] - kz_sub) / (kz[last] + kz_sub),
                Polarization::P => {
                    let a = kz[last] / n_all[last] * n_sub;
                    let b = kz_sub / n_sub * n_all[last];
                    (a - b) / (a + b)
                }
                Polarization::Unpolarized => {
                    return Err(XrayDbError::DataError(
                        "use S or P polarization for multilayer".to_string(),
                    ));
                }
            };

            // Substrate roughness
            if substrate_rough >= 1e-12 {
                r_amp =
                    r_amp * (c(-2.0 * substrate_rough * substrate_rough) * kz[last] * kz_sub).exp();
            }

            // Recurse upward through layers
            for i in (0..last).rev() {
                let fresnel_r = match polarization {
                    Polarization::S => (kz[i] - kz[i + 1]) / (kz[i] + kz[i + 1]),
                    Polarization::P => {
                        let a = kz[i] / n_all[i] * n_all[i + 1];
                        let b = kz[i + 1] / n_all[i + 1] * n_all[i];
                        (a - b) / (a + b)
                    }
                    Polarization::Unpolarized => unreachable!(),
                };
                let p2 = (two_i * c(t_all[i + 1]) * kz[i + 1]).exp();
                r_amp = (fresnel_r + r_amp * p2) / (one + fresnel_r * r_amp * p2);
            }

            // Surface (air to first layer)
            let fresnel_r = match polarization {
                Polarization::S => (kiz - kz[0]) / (kiz + kz[0]),
                Polarization::P => (kiz - kz[0] / n_all[0]) / (kiz + kz[0] / n_all[0]),
                Polarization::Unpolarized => unreachable!(),
            };
            let p2 = (two_i * c(t_all[0]) * kz[0]).exp();
            r_amp = (fresnel_r + r_amp * p2) / (one + fresnel_r * r_amp * p2);

            // Surface roughness
            if surface_rough >= 1e-12 {
                r_amp =
                    r_amp * (c(-2.0 * surface_rough * surface_rough) * kiz * kz[0]).exp();
            }

            result.push((r_amp * r_amp.conj()).re);
        }

        Ok(result)
    }

    /// Reflectivity for a coated mirror (convenience wrapper around multilayer).
    ///
    /// # Arguments
    /// * `coating` - Coating material formula
    /// * `coating_thick` - Coating thickness in Å
    /// * `substrate` - Substrate material formula
    /// * `theta` - Grazing angles in radians
    /// * `energy` - X-ray energy in eV
    /// * `coating_density` - Coating density in g/cm³
    /// * `surface_roughness` - Air-coating roughness in Å
    /// * `substrate_density` - Substrate density in g/cm³
    /// * `substrate_roughness` - Substrate-coating roughness in Å
    /// * `binder` - Optional binder layer: (material, thickness_Å, density)
    /// * `polarization` - S or P polarization
    pub fn coated_reflectivity(
        &self,
        coating: &str,
        coating_thick: f64,
        substrate: &str,
        theta: &[f64],
        energy: f64,
        coating_density: f64,
        surface_roughness: f64,
        substrate_density: f64,
        substrate_roughness: f64,
        binder: Option<(&str, f64, f64)>,
        polarization: Polarization,
    ) -> Result<Vec<f64>> {
        let (stackup, thickness, densities): (Vec<&str>, Vec<f64>, Vec<f64>) = match binder {
            Some((binder_mat, binder_thick, binder_dens)) => (
                vec![coating, binder_mat],
                vec![coating_thick, binder_thick],
                vec![coating_density, binder_dens],
            ),
            None => (
                vec![coating],
                vec![coating_thick],
                vec![coating_density],
            ),
        };

        self.multilayer_reflectivity(
            &stackup,
            &thickness,
            substrate,
            theta,
            energy,
            1,
            &densities,
            substrate_density,
            substrate_roughness,
            surface_roughness,
            polarization,
        )
    }
}
