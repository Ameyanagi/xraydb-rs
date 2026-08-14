//! Internal interpolation primitives.
//!
//! These are total functions: degenerate inputs (empty tables, mismatched lengths)
//! yield `NaN` rather than panicking. The shipped database never produces such inputs,
//! so this is defence in depth against a corrupted blob.

/// Linear interpolation of a single value, equivalent to `numpy.interp`.
///
/// Values outside `[xp[0], xp[len-1]]` are clamped to the boundary values.
/// Returns `NaN` if the table is empty.
pub(crate) fn interp_one(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    let n = xp.len().min(fp.len());
    if n == 0 {
        return f64::NAN;
    }
    if n == 1 || x <= xp[0] {
        return fp[0];
    }
    if x >= xp[n - 1] {
        return fp[n - 1];
    }

    // Binary search for the bracket. `x` is strictly inside the table here, so
    // `idx` lands in `1..n`.
    let idx = xp[..n].partition_point(|&v| v < x);
    if idx == 0 {
        return fp[0];
    }
    if idx >= n {
        return fp[n - 1];
    }

    // Exact hit on a knot.
    if (xp[idx] - x).abs() < f64::EPSILON * xp[idx].abs() {
        return fp[idx];
    }

    let lo = idx - 1;
    let span = xp[idx] - xp[lo];
    if span == 0.0 {
        return fp[idx];
    }
    let t = (x - xp[lo]) / span;
    fp[lo] + t * (fp[idx] - fp[lo])
}

/// Linear interpolation over a slice of query points.
#[cfg(test)]
pub(crate) fn interp(x: &[f64], xp: &[f64], fp: &[f64]) -> Vec<f64> {
    x.iter().map(|&xi| interp_one(xi, xp, fp)).collect()
}

/// Log-log linear interpolation against **pre-computed** log-space tables.
///
/// Equivalent to `exp(interp(ln(x), log_xp, log_fp))`, but without transforming the
/// tables on every call. `log_xp` and `log_fp` must already hold `ln()` values; see
/// [`crate::db`] for where they are built.
pub(crate) fn interp_loglog_pre(x: f64, log_xp: &[f64], log_fp: &[f64]) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    interp_one(x.ln(), log_xp, log_fp).exp()
}

/// Natural log of a value, floored away from zero so `ln` stays finite.
///
/// Mirrors the upstream XrayDB treatment of tabulated zeros in log-space data.
pub(crate) fn safe_ln(v: f64) -> f64 {
    const FLOOR: f64 = 1e-99;
    if v.abs() < FLOOR {
        FLOOR.ln()
    } else {
        v.abs().ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interp_basic() {
        let xp = [0.0, 1.0, 2.0];
        let fp = [0.0, 10.0, 20.0];

        let result = interp(&[0.5, 1.5], &xp, &fp);
        assert!((result[0] - 5.0).abs() < 1e-10);
        assert!((result[1] - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_interp_clamping() {
        let xp = [1.0, 2.0, 3.0];
        let fp = [10.0, 20.0, 30.0];

        let result = interp(&[0.0, 4.0], &xp, &fp);
        assert!((result[0] - 10.0).abs() < 1e-10);
        assert!((result[1] - 30.0).abs() < 1e-10);
    }

    #[test]
    fn empty_table_yields_nan_not_panic() {
        assert!(interp_one(1.0, &[], &[]).is_nan());
        assert!(interp_one(1.0, &[1.0], &[]).is_nan());
    }

    #[test]
    fn single_point_table_returns_that_point() {
        assert_eq!(interp_one(-5.0, &[1.0], &[7.0]), 7.0);
        assert_eq!(interp_one(500.0, &[1.0], &[7.0]), 7.0);
    }

    #[test]
    fn mismatched_lengths_use_the_shorter() {
        // Longer fp than xp must not read past the end of xp.
        assert_eq!(interp_one(10.0, &[0.0, 1.0], &[0.0, 5.0, 9.0]), 5.0);
    }

    #[test]
    fn duplicate_knots_do_not_divide_by_zero() {
        let v = interp_one(1.0, &[0.0, 1.0, 1.0, 2.0], &[0.0, 3.0, 4.0, 5.0]);
        assert!(v.is_finite(), "got {v}");
    }

    #[test]
    fn loglog_matches_naive_transform() {
        let xp: [f64; 3] = [1.0, 10.0, 100.0];
        let fp: [f64; 3] = [1.0, 100.0, 10000.0];
        let log_xp: Vec<f64> = xp.iter().map(|v| v.ln()).collect();
        let log_fp: Vec<f64> = fp.iter().map(|v| v.ln()).collect();

        // Exactly a power law, so log-log interpolation is exact.
        let got = interp_loglog_pre(50.0, &log_xp, &log_fp);
        assert!((got - 2500.0).abs() < 1e-6, "got {got}");
    }

    #[test]
    fn loglog_rejects_non_positive_x() {
        assert!(interp_loglog_pre(0.0, &[0.0, 1.0], &[0.0, 1.0]).is_nan());
        assert!(interp_loglog_pre(-1.0, &[0.0, 1.0], &[0.0, 1.0]).is_nan());
    }

    #[test]
    fn safe_ln_floors_zero() {
        assert!(safe_ln(0.0).is_finite());
        assert_eq!(safe_ln(0.0), 1e-99_f64.ln());
        assert!((safe_ln(std::f64::consts::E) - 1.0).abs() < 1e-12);
    }
}
