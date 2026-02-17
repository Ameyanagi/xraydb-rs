use std::path::Path;

use xraydb_data::ChantlerRecord;

/// Parse Chantler data from per-element .dat files in a directory.
///
/// Each file is named `01.dat` through `92.dat` (Z=1 to 92).
/// Format per file:
///   - Line 1 (comment): `# <Elem>:  Z = <z>;  ...  density ... = <density> g/cm3`
///   - Line 2 (comment): `# ... sigma_mu = <sigma_mu>`
///   - Line 3 (comment): `# ... f2 = <mue_f2>`
///   - Further comment lines may contain relativistic/nuclear corrections
///   - Data lines: energy(keV)  f1(e/atom)  f2(e/atom)  mu_photo  mu_incoh  [mu_total]
pub fn parse_chantler(dir: &Path) -> Vec<ChantlerRecord> {
    let mut records = Vec::new();

    for z in 1..=92u16 {
        let fname = dir.join(format!("{:02}.dat", z));
        let content = std::fs::read_to_string(&fname)
            .unwrap_or_else(|_| panic!("failed to read Chantler file {:?}", fname));
        let lines: Vec<&str> = content.lines().collect();

        // Line 1: extract element symbol and density
        let line1 = lines[0].trim_start_matches('#').trim();
        let (elem, density) = parse_chantler_line1(line1);

        // Line 2: extract sigma_mu (last value on the line)
        let sigma_mu = parse_last_float(lines[1]);

        // Line 3: extract mue_f2 (last value on the line)
        let mue_f2 = parse_last_float(lines[2]);

        // Scan comment lines for relativistic and nuclear corrections
        let mut corr_henke: f64 = 0.0;
        let mut corr_cl35: f64 = 0.0;
        let mut corr_nucl: f64 = 0.0;

        for line in &lines {
            if !line.starts_with('#') {
                continue;
            }
            if line.contains("Relativistic") || line.contains("Nuclear Thomson") {
                // Strip # and parse "label = value"
                let stripped = line
                    .trim_start_matches('#')
                    .replace('#', " ")
                    .trim()
                    .to_string();

                if let Some(eq_pos) = stripped.find('=') {
                    let val_str = stripped[eq_pos + 1..]
                        .replace(',', "")
                        .replace("e/atom", "")
                        .trim()
                        .to_string();

                    if stripped.contains("Relativistic") {
                        let parts: Vec<f64> = val_str
                            .split_whitespace()
                            .filter_map(|w| w.parse().ok())
                            .collect();
                        if parts.len() >= 2 {
                            corr_henke = parts[0];
                            corr_cl35 = parts[1];
                        }
                    } else {
                        corr_nucl = val_str
                            .split_whitespace()
                            .next()
                            .and_then(|w| w.parse().ok())
                            .unwrap_or(0.0);
                    }
                }
            }
        }

        // Parse data lines
        let mut energy = Vec::new();
        let mut f1 = Vec::new();
        let mut f2 = Vec::new();
        let mut mu_photo = Vec::new();
        let mut mu_incoh = Vec::new();
        let mut mu_total = Vec::new();

        for line in &lines {
            if line.starts_with('#') {
                continue;
            }
            let words: Vec<f64> = line
                .split_whitespace()
                .filter_map(|w| w.parse().ok())
                .collect();
            if words.len() >= 5 {
                // Energy: convert keV to eV
                energy.push(1000.0 * words[0]);
                // f1: apply corrections (matching Python: words[1] - z + corr_cl35 + corr_nucl)
                f1.push(words[1] - z as f64 + corr_cl35 + corr_nucl);
                f2.push(words[2]);
                mu_photo.push(words[3]);
                mu_incoh.push(words[4]);
                // mu_total = mu_photo + mu_incoh (matching Python)
                mu_total.push(words[3] + words[4]);
            }
        }

        records.push(ChantlerRecord {
            element: elem,
            sigma_mu,
            mue_f2,
            density,
            corr_henke,
            corr_cl35,
            corr_nucl,
            energy,
            f1,
            f2,
            mu_photo,
            mu_incoh,
            mu_total,
        });
    }

    records
}

fn parse_chantler_line1(line: &str) -> (String, f64) {
    // Format: "Fe:    Z = 26;   Atomic weight = 55.84700 g/mol ;    nominal density:   [INLINE] = 7.8600E+00 g/cm3"
    // We need the element symbol (before ':') and density (second-to-last word before "g/cm3")
    let elem = line.split(':').next().unwrap().trim().to_string();

    // Find density: look for the pattern "= <number> g/cm3" near the end
    let mut words: Vec<&str> = line.split_whitespace().collect();
    // Remove trailing "g/cm3" if present
    if words.last().is_some_and(|w| w.contains("g/cm")) {
        words.pop();
    }
    // The density value should now be the last word
    let density: f64 = words.last().and_then(|w| w.parse().ok()).unwrap_or(0.0);

    (elem, density)
}

fn parse_last_float(line: &str) -> f64 {
    line.split_whitespace()
        .rev()
        .find_map(|w| w.parse::<f64>().ok())
        .unwrap_or(0.0)
}
