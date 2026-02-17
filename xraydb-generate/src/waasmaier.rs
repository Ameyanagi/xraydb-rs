use std::path::Path;

use xraydb_data::WaasmaierRecord;

/// Parse `waasmaeir_kirfel.dat` for f0 elastic scattering factors.
///
/// Format: sections starting with `#S <Z> <ion>`, followed by
/// `#N 11`, `#L ...`, then a data line with 11 floats:
///   a1 a2 a3 a4 a5 c b1 b2 b3 b4 b5
/// Stored as: offset=c, scale=[a1..a5], exponents=[b1..b5]
pub fn parse_waasmaier(path: &Path) -> Vec<WaasmaierRecord> {
    let content = std::fs::read_to_string(path).expect("failed to read waasmaier data");
    let lines: Vec<&str> = content.lines().collect();

    // Validate header
    if lines.len() < 2 || !lines[1].contains("Elastic Photon-Atom Scatt") {
        panic!("Source file not recognized for f0_WaasKirf data");
    }

    let mut records = Vec::new();
    let mut i = 0;

    while i < lines.len() {
        let line = lines[i];
        if let Some(stripped) = line.strip_prefix("#S ") {
            let parts: Vec<&str> = stripped.split_whitespace().collect();
            if parts.len() >= 2 {
                let atno: u16 = parts[0].parse().unwrap();
                let ion = parts[1].to_string();

                // Skip 3 lines (#N, #L, then data)
                i += 3;
                if i >= lines.len() {
                    break;
                }

                let data_line = lines[i];
                let words: Vec<f64> = data_line
                    .split_whitespace()
                    .filter_map(|w| w.parse().ok())
                    .collect();

                if words.len() >= 11 {
                    let scale = words[0..5].to_vec();
                    let offset = words[5];
                    let exponents = words[6..11].to_vec();

                    // Extract element symbol from ion name:
                    // strip digits, +, - then trim suffixes "va", "val"
                    let elem = extract_element(&ion);

                    records.push(WaasmaierRecord {
                        atomic_number: atno,
                        element: elem,
                        ion,
                        offset,
                        scale,
                        exponents,
                    });
                }
            }
        }
        i += 1;
    }

    records
}

fn extract_element(ion: &str) -> String {
    // Strip digits, +, - from ion name to get element symbol
    let stripped: String = ion
        .chars()
        .filter(|c| !c.is_ascii_digit() && *c != '+' && *c != '-')
        .collect();
    let mut elem = stripped.trim().to_string();

    // Strip common suffixes
    for suffix in &["val", "va"] {
        if elem.ends_with(suffix) {
            elem.truncate(elem.len() - suffix.len());
        }
    }

    elem
}
