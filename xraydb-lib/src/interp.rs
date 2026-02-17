/// Linear interpolation (equivalent to numpy.interp).
///
/// Interpolates values from `(xp, fp)` at points `x`.
/// Values outside the range are clamped to the boundary values.
pub fn interp(x: &[f64], xp: &[f64], fp: &[f64]) -> Vec<f64> {
    x.iter()
        .map(|&xi| interp_one(xi, xp, fp))
        .collect()
}

/// Interpolate a single value.
pub fn interp_one(x: f64, xp: &[f64], fp: &[f64]) -> f64 {
    if x <= xp[0] {
        return fp[0];
    }
    if x >= xp[xp.len() - 1] {
        return fp[fp.len() - 1];
    }

    // Binary search for the bracket
    let idx = xp.partition_point(|&v| v < x);
    if idx == 0 {
        return fp[0];
    }

    // Check for exact match
    if (xp[idx] - x).abs() < f64::EPSILON * xp[idx].abs() {
        return fp[idx];
    }

    let lo = idx - 1;
    let t = (x - xp[lo]) / (xp[idx] - xp[lo]);
    fp[lo] + t * (fp[idx] - fp[lo])
}

/// Log-log linear interpolation.
///
/// Equivalent to `exp(interp(log(x), log(xp), log(fp)))`.
/// Used for Chantler f2 and mu data.
pub fn interp_loglog(x: &[f64], xp: &[f64], fp: &[f64]) -> Vec<f64> {
    let log_xp: Vec<f64> = xp.iter().map(|v| v.ln()).collect();
    let log_fp: Vec<f64> = fp.iter().map(|v| v.ln()).collect();
    let log_x: Vec<f64> = x.iter().map(|v| v.ln()).collect();
    interp(&log_x, &log_xp, &log_fp)
        .into_iter()
        .map(|v| v.exp())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interp_basic() {
        let xp = vec![0.0, 1.0, 2.0];
        let fp = vec![0.0, 10.0, 20.0];

        let result = interp(&[0.5, 1.5], &xp, &fp);
        assert!((result[0] - 5.0).abs() < 1e-10);
        assert!((result[1] - 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_interp_clamping() {
        let xp = vec![1.0, 2.0, 3.0];
        let fp = vec![10.0, 20.0, 30.0];

        let result = interp(&[0.0, 4.0], &xp, &fp);
        assert!((result[0] - 10.0).abs() < 1e-10);
        assert!((result[1] - 30.0).abs() < 1e-10);
    }
}
