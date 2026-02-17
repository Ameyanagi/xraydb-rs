/// Cubic spline interpolation using pre-computed second derivatives (Elam method).
///
/// This is the core interpolation used for all Elam photoabsorption and
/// scattering cross-section lookups. Ported from the Python `elam_spline` function.
///
/// # Arguments
/// * `xin` - Input x values (must be strictly increasing)
/// * `yin` - Input y values
/// * `yspl` - Pre-computed spline coefficients (second derivatives of y)
/// * `xout` - Output x values to interpolate at
///
/// # Returns
/// Interpolated y values at each xout point
pub fn elam_spline(xin: &[f64], yin: &[f64], yspl: &[f64], xout: &[f64]) -> Vec<f64> {
    xout.iter()
        .map(|&x| {
            // Find bracket indices using binary search
            let hi = match xin.partition_point(|&v| v < x) {
                i if i >= xin.len() => xin.len() - 1,
                0 => 1.min(xin.len() - 1),
                i => i,
            };
            let lo = hi - 1;

            let diff = xin[hi] - xin[lo];
            debug_assert!(diff > 0.0, "xin must be strictly increasing");

            let a = (xin[hi] - x) / diff;
            let b = (x - xin[lo]) / diff;

            a * yin[lo]
                + b * yin[hi]
                + (diff * diff / 6.0)
                    * ((a * a - 1.0) * a * yspl[lo] + (b * b - 1.0) * b * yspl[hi])
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spline_at_knot_points() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let y = vec![1.0, 4.0, 9.0, 16.0, 25.0];
        let spl = vec![0.0; 5]; // zero second derivatives = linear

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
        let x = vec![0.0, 1.0, 2.0];
        let y = vec![0.0, 1.0, 2.0];
        let spl = vec![0.0; 3];

        let result = elam_spline(&x, &y, &spl, &[0.5, 1.5]);
        assert!((result[0] - 0.5).abs() < 1e-10);
        assert!((result[1] - 1.5).abs() < 1e-10);
    }
}
