//! Parsing for the `--energy` / `--q` specification syntax.

use anyhow::{Context, Result, bail, ensure};

/// Parse an energy (or q) specification into a list of values.
///
/// Four forms are accepted:
///
/// | Syntax | Meaning |
/// |---|---|
/// | `10000` | a single value |
/// | `5000,7112,10000` | an explicit list |
/// | `5000:15000:100` | `start:stop:step`, inclusive of `stop` when it lands on the grid |
/// | `5000:15000/50` | `start:stop/count`, log-spaced (linear if `start <= 0`) |
pub fn parse_spec(spec: &str) -> Result<Vec<f64>> {
    let spec = spec.trim();
    ensure!(!spec.is_empty(), "empty energy specification");

    if let Some((range, count)) = spec.split_once('/') {
        let (start, stop) = parse_range(range)?;
        let count: usize = count
            .trim()
            .parse()
            .with_context(|| format!("invalid point count '{count}'"))?;
        ensure!(count >= 1, "point count must be at least 1");
        return Ok(spaced(start, stop, count));
    }

    if spec.contains(':') {
        let parts: Vec<&str> = spec.split(':').collect();
        ensure!(
            parts.len() == 3,
            "range must be start:stop:step or start:stop/count, got '{spec}'"
        );
        let start = number(parts[0])?;
        let stop = number(parts[1])?;
        let step = number(parts[2])?;
        ensure!(step > 0.0, "step must be > 0, got {step}");
        ensure!(stop >= start, "stop ({stop}) must be >= start ({start})");

        let n = ((stop - start) / step).floor() as usize + 1;
        ensure!(n <= 1_000_000, "range would produce {n} points; narrow it");
        return Ok((0..n).map(|i| start + step * i as f64).collect());
    }

    let values: Result<Vec<f64>> = spec.split(',').map(number).collect();
    let values = values?;
    ensure!(!values.is_empty(), "no values in '{spec}'");
    Ok(values)
}

fn parse_range(range: &str) -> Result<(f64, f64)> {
    let (start, stop) = range
        .split_once(':')
        .with_context(|| format!("expected start:stop before '/', got '{range}'"))?;
    let start = number(start)?;
    let stop = number(stop)?;
    ensure!(stop >= start, "stop ({stop}) must be >= start ({start})");
    Ok((start, stop))
}

/// `count` points from `start` to `stop`, log-spaced when both are positive.
fn spaced(start: f64, stop: f64, count: usize) -> Vec<f64> {
    if count == 1 {
        return vec![start];
    }
    let n = (count - 1) as f64;
    if start > 0.0 && stop > 0.0 {
        let (ls, le) = (start.ln(), stop.ln());
        (0..count)
            .map(|i| (ls + (le - ls) * i as f64 / n).exp())
            .collect()
    } else {
        (0..count)
            .map(|i| start + (stop - start) * i as f64 / n)
            .collect()
    }
}

fn number(text: &str) -> Result<f64> {
    let text = text.trim();
    let value: f64 = text
        .parse()
        .with_context(|| format!("'{text}' is not a number"))?;
    if !value.is_finite() {
        bail!("'{text}' is not a finite number");
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_value() {
        assert_eq!(parse_spec("10000").unwrap(), vec![10000.0]);
        assert_eq!(parse_spec(" 7112 ").unwrap(), vec![7112.0]);
    }

    #[test]
    fn comma_list() {
        assert_eq!(
            parse_spec("5000,7112,10000").unwrap(),
            vec![5000.0, 7112.0, 10000.0]
        );
        assert_eq!(parse_spec("1, 2 ,3").unwrap(), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn start_stop_step() {
        assert_eq!(
            parse_spec("1000:1400:100").unwrap(),
            vec![1000.0, 1100.0, 1200.0, 1300.0, 1400.0]
        );
    }

    #[test]
    fn start_stop_step_handles_non_divisible_ranges() {
        let v = parse_spec("0:10:3").unwrap();
        assert_eq!(v, vec![0.0, 3.0, 6.0, 9.0]);
    }

    #[test]
    fn log_spaced_count() {
        let v = parse_spec("100:10000/3").unwrap();
        assert_eq!(v.len(), 3);
        assert!((v[0] - 100.0).abs() < 1e-9);
        assert!((v[1] - 1000.0).abs() < 1e-6, "got {}", v[1]);
        assert!((v[2] - 10000.0).abs() < 1e-6);
    }

    #[test]
    fn log_spacing_falls_back_to_linear_at_zero() {
        let v = parse_spec("0:10/3").unwrap();
        assert_eq!(v, vec![0.0, 5.0, 10.0]);
    }

    #[test]
    fn single_point_range() {
        assert_eq!(parse_spec("500:5000/1").unwrap(), vec![500.0]);
    }

    #[test]
    fn scientific_notation_is_accepted() {
        assert_eq!(parse_spec("1e4").unwrap(), vec![10000.0]);
        assert_eq!(parse_spec("1.5e3,2e3").unwrap(), vec![1500.0, 2000.0]);
    }

    #[test]
    fn rejects_malformed_input() {
        for bad in [
            "",
            "  ",
            "abc",
            "1000:2000",
            "1000:2000:3000:4000",
            "1000:2000:0",
            "1000:2000:-5",
            "2000:1000:10",
            "1000:2000/0",
            "1000:2000/abc",
            "1,,2",
        ] {
            assert!(parse_spec(bad).is_err(), "'{bad}' should be rejected");
        }
    }

    #[test]
    fn rejects_absurdly_large_ranges() {
        assert!(parse_spec("0:1e12:1").is_err());
    }

    #[test]
    fn rejects_non_finite() {
        assert!(parse_spec("inf").is_err());
        assert!(parse_spec("NaN").is_err());
    }
}
