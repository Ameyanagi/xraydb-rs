//! Internal cubic-spline evaluation for the Elam tables.

/// Evaluate a cubic spline at one point using pre-computed second derivatives.
///
/// This is the core interpolation behind every Elam photoabsorption and scattering
/// cross-section lookup; it is a direct port of upstream XrayDB's `elam_spline`.
///
/// * `xin` — knot abscissae, strictly increasing
/// * `yin` — knot ordinates
/// * `yspl` — pre-computed second derivatives at the knots
///
/// Returns `NaN` for a degenerate (empty or single-knot) table rather than panicking.
pub(crate) fn elam_spline_one(xin: &[f64], yin: &[f64], yspl: &[f64], x: f64) -> f64 {
    let n = xin.len().min(yin.len()).min(yspl.len());
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 {
        return yin[0];
    }

    // Bracket [lo, hi] with 1 <= hi <= n-1.
    let hi = xin[..n].partition_point(|&v| v < x).clamp(1, n - 1);
    let lo = hi - 1;

    let diff = xin[hi] - xin[lo];
    if diff.is_nan() || diff <= 0.0 {
        // Non-increasing knots: fall back to the nearest ordinate instead of
        // dividing by zero.
        return yin[hi];
    }

    let a = (xin[hi] - x) / diff;
    let b = (x - xin[lo]) / diff;

    a * yin[lo]
        + b * yin[hi]
        + (diff * diff / 6.0) * ((a * a - 1.0) * a * yspl[lo] + (b * b - 1.0) * b * yspl[hi])
}

/// Evaluate a cubic spline at each point of `xout`.
#[cfg(test)]
pub(crate) fn elam_spline(xin: &[f64], yin: &[f64], yspl: &[f64], xout: &[f64]) -> Vec<f64> {
    xout.iter()
        .map(|&x| elam_spline_one(xin, yin, yspl, x))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spline_at_knot_points() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0];
        let y = [1.0, 4.0, 9.0, 16.0, 25.0];
        let spl = [0.0; 5]; // zero second derivatives = linear

        for (&xi, &yi) in x.iter().zip(y.iter()) {
            let result = elam_spline(&x, &y, &spl, &[xi]);
            assert!(
                (result[0] - yi).abs() < 1e-10,
                "at x={xi}: got {} expected {yi}",
                result[0]
            );
        }
    }

    #[test]
    fn test_spline_linear_interpolation() {
        // With zero spline coefficients, should do linear interpolation
        let x = [0.0, 1.0, 2.0];
        let y = [0.0, 1.0, 2.0];
        let spl = [0.0; 3];

        let result = elam_spline(&x, &y, &spl, &[0.5, 1.5]);
        assert!((result[0] - 0.5).abs() < 1e-10);
        assert!((result[1] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn empty_table_yields_nan_not_panic() {
        let out = elam_spline(&[], &[], &[], &[1.0]);
        assert!(out[0].is_nan());
    }

    #[test]
    fn single_knot_returns_that_knot() {
        assert_eq!(elam_spline_one(&[2.0], &[9.0], &[0.0], 100.0), 9.0);
    }

    #[test]
    fn mismatched_lengths_do_not_read_out_of_bounds() {
        // yspl shorter than xin: n collapses to 2, no panic.
        let v = elam_spline_one(&[0.0, 1.0, 2.0], &[0.0, 1.0, 2.0], &[0.0, 0.0], 0.5);
        assert!((v - 0.5).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn extrapolates_beyond_the_table_without_panicking() {
        let x = [0.0, 1.0, 2.0];
        let y = [0.0, 1.0, 2.0];
        let spl = [0.0; 3];
        assert!(elam_spline_one(&x, &y, &spl, -10.0).is_finite());
        assert!(elam_spline_one(&x, &y, &spl, 10.0).is_finite());
    }
}
