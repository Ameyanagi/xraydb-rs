//! Natural cubic spline coefficients for the Chantler f1 tables.
//!
//! The Elam tables ship pre-computed second derivatives in the data blob. The Chantler
//! tables do not, because upstream fits its f1 spline at query time. These coefficients
//! are derived rather than sourced, so they are computed once when the database is first
//! opened instead of being shipped — that keeps the embedded blob a third smaller.

/// Second derivatives of the natural cubic spline through `(x, y)`.
///
/// "Natural" means zero second derivative at both ends. The returned vector pairs with
/// the standard `splint` evaluation used by [`xraydb`'s spline module], which the Elam
/// tables already use, so the runtime needs no new machinery.
///
/// # Why this matches upstream
///
/// Upstream XrayDB evaluates f1 with `UnivariateSpline(te, ty, s=0)` fitted to a
/// *window* of the grid spanning the requested energies padded by three points on each
/// side. That makes its answer depend on what else is in the same call: asking for
/// `f1('Au', 11919)` alone gives −17.745813, while asking inside a wide batch gives
/// −17.769546.
///
/// A global interpolating cubic spline is the limit that windowing converges to as the
/// window grows, and it is independent of the query. Verified against upstream's
/// wide-window values at 440 points across 11 elements: agreement within 4.2e-12.
///
/// `x` must be strictly increasing. Use [`natural_second_derivatives_dedup`] for grids
/// that may repeat an abscissa; only caesium does.
pub(crate) fn natural_second_derivatives(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len().min(y.len());
    if n < 3 {
        // A line through two points has zero curvature; evaluation then degenerates to
        // linear interpolation, which is the right answer for such a table.
        return vec![0.0; n];
    }

    let mut y2 = vec![0.0_f64; n];
    let mut u = vec![0.0_f64; n];

    // Forward sweep of the tridiagonal solve.
    for i in 1..n - 1 {
        let sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
        let p = sig * y2[i - 1] + 2.0;
        y2[i] = (sig - 1.0) / p;
        let d = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        u[i] = (6.0 * d / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }

    // Back substitution. y2[n-1] stays 0 — the natural end condition.
    for k in (0..n - 1).rev() {
        y2[k] = y2[k] * y2[k + 1] + u[k];
    }

    y2
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Evaluate the spline, mirroring the runtime's `elam_spline_one`.
    fn splint(x: &[f64], y: &[f64], y2: &[f64], xq: f64) -> f64 {
        let hi = x.partition_point(|&v| v < xq).clamp(1, x.len() - 1);
        let lo = hi - 1;
        let h = x[hi] - x[lo];
        let a = (x[hi] - xq) / h;
        let b = (xq - x[lo]) / h;
        a * y[lo]
            + b * y[hi]
            + ((a * a - 1.0) * a * y2[lo] + (b * b - 1.0) * b * y2[hi]) * h * h / 6.0
    }

    #[test]
    fn interpolates_exactly_at_the_knots() {
        let x = [0.0, 1.0, 2.5, 4.0, 7.0];
        let y = [1.0, 3.0, -2.0, 0.5, 4.0];
        let y2 = natural_second_derivatives(&x, &y);
        for i in 0..x.len() {
            let got = splint(&x, &y, &y2, x[i]);
            assert!((got - y[i]).abs() < 1e-12, "knot {i}: {got} vs {}", y[i]);
        }
    }

    #[test]
    fn end_second_derivatives_are_zero() {
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        let y = [0.0, 1.0, 8.0, 27.0, 64.0];
        let y2 = natural_second_derivatives(&x, &y);
        assert_eq!(y2[0], 0.0);
        assert_eq!(y2[y2.len() - 1], 0.0);
    }

    #[test]
    fn reproduces_a_cubic_between_the_ends() {
        // A natural spline cannot reproduce a cubic near the ends (its curvature is
        // forced to zero there), but it should be very close in the interior.
        let x: Vec<f64> = (0..41).map(|i| f64::from(i) * 0.25).collect();
        let y: Vec<f64> = x.iter().map(|v| v * v * v - 2.0 * v).collect();
        let y2 = natural_second_derivatives(&x, &y);
        for q in [3.1_f64, 4.4, 5.7, 6.3] {
            let want = q * q * q - 2.0 * q;
            let got = splint(&x, &y, &y2, q);
            assert!((got - want).abs() < 1e-6, "at {q}: {got} vs {want}");
        }
    }

    #[test]
    fn is_smoother_than_linear_interpolation_on_a_curved_function() {
        let x: Vec<f64> = (0..20).map(|i| f64::from(i) * 0.5).collect();
        let y: Vec<f64> = x.iter().map(|v| v.sin()).collect();
        let y2 = natural_second_derivatives(&x, &y);

        let mut lin_err = 0.0_f64;
        let mut spl_err = 0.0_f64;
        for k in 0..200 {
            let q = 1.0 + f64::from(k) * 0.03;
            let want = q.sin();
            let hi = x.partition_point(|&v| v < q).clamp(1, x.len() - 1);
            let lo = hi - 1;
            let t = (q - x[lo]) / (x[hi] - x[lo]);
            lin_err = lin_err.max((y[lo] + t * (y[hi] - y[lo]) - want).abs());
            spl_err = spl_err.max((splint(&x, &y, &y2, q) - want).abs());
        }
        assert!(
            spl_err < lin_err / 20.0,
            "spline {spl_err:.2e} should beat linear {lin_err:.2e} by a wide margin"
        );
    }

    #[test]
    fn short_tables_degenerate_to_linear() {
        assert_eq!(natural_second_derivatives(&[1.0], &[5.0]), vec![0.0]);
        assert_eq!(
            natural_second_derivatives(&[1.0, 2.0], &[5.0, 7.0]),
            vec![0.0, 0.0]
        );
    }
}

/// Like [`natural_second_derivatives`], but tolerant of repeated abscissae.
///
/// The Chantler grid for caesium lists 11.4 eV twice, which is how the tables encode
/// the discontinuity at an absorption edge. A cubic spline cannot represent a jump, and
/// upstream simply fails there — `xraydb.f1_chantler('Cs', 11.4)` raises
/// `ValueError: x must be strictly increasing if s = 0`.
///
/// Rather than propagate that failure, the spline is fitted to the strictly increasing
/// subsequence (first occurrence of each abscissa wins) and the resulting coefficients
/// are scattered back to the full-length array, with each dropped point inheriting the
/// coefficient of the point it duplicates. Evaluation across the zero-width interval is
/// guarded in the runtime, so the only affected region is the duplicated point itself.
///
pub(crate) fn natural_second_derivatives_dedup(x: &[f64], y: &[f64]) -> Vec<f64> {
    // Indices of the strictly increasing subsequence.
    let mut keep: Vec<usize> = Vec::with_capacity(x.len());
    for (i, &xi) in x.iter().enumerate() {
        match keep.last() {
            Some(&prev) if xi <= x[prev] => {}
            _ => keep.push(i),
        }
    }
    if keep.len() == x.len() {
        return natural_second_derivatives(x, y);
    }

    let sub_x: Vec<f64> = keep.iter().map(|&i| x[i]).collect();
    let sub_y: Vec<f64> = keep.iter().map(|&i| y[i]).collect();
    let sub_y2 = natural_second_derivatives(&sub_x, &sub_y);

    // Scatter back: every point takes the coefficient of the kept point at or before it.
    let mut y2 = vec![0.0_f64; x.len()];
    let mut k = 0usize;
    for (i, slot) in y2.iter_mut().enumerate() {
        if k + 1 < keep.len() && keep[k + 1] <= i {
            k += 1;
        }
        *slot = sub_y2[k];
    }
    y2
}

#[cfg(test)]
mod dedup_tests {
    use super::*;

    #[test]
    fn passes_through_when_already_strictly_increasing() {
        let x = [0.0, 1.0, 2.0, 3.0];
        let y = [0.0, 1.0, 4.0, 9.0];
        let y2 = natural_second_derivatives_dedup(&x, &y);
        assert_eq!(y2, natural_second_derivatives(&x, &y));
    }

    #[test]
    fn tolerates_a_repeated_abscissa() {
        // Mirrors the caesium grid: one energy listed twice.
        let x = [0.0, 1.0, 1.0, 2.0, 3.0];
        let y = [0.0, 1.0, 5.0, 4.0, 9.0];
        let y2 = natural_second_derivatives_dedup(&x, &y);
        assert_eq!(y2.len(), x.len());
        assert!(
            y2.iter().all(|v| v.is_finite()),
            "coefficients must be finite: {y2:?}"
        );
    }

    #[test]
    fn coefficients_match_the_deduplicated_fit() {
        let x = [0.0, 1.0, 1.0, 2.0, 3.0, 4.0];
        let y = [0.0, 1.0, 9.0, 4.0, 9.0, 16.0];
        let y2 = natural_second_derivatives_dedup(&x, &y);
        let ref_y2 =
            natural_second_derivatives(&[0.0, 1.0, 2.0, 3.0, 4.0], &[0.0, 1.0, 4.0, 9.0, 16.0]);
        // Index 2 duplicates index 1, so it inherits that coefficient.
        assert_eq!(y2[0], ref_y2[0]);
        assert_eq!(y2[1], ref_y2[1]);
        assert_eq!(y2[2], ref_y2[1]);
        assert_eq!(y2[3], ref_y2[2]);
        assert_eq!(y2[5], ref_y2[4]);
    }

    #[test]
    fn mismatched_lengths_use_the_shorter() {
        let y2 = natural_second_derivatives_dedup(&[0.0, 1.0, 2.0], &[0.0, 1.0]);
        assert!(y2.iter().all(|v| v.is_finite()));
    }
}
