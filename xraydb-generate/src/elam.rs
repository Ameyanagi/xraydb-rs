use std::path::Path;

use xraydb_data::{
    CosterKronigRecord, PhotoabsorptionRecord, ScatteringRecord, XrayLevelRecord,
    XrayTransitionRecord,
};

type ElamParseResult = (
    Vec<XrayLevelRecord>,
    Vec<XrayTransitionRecord>,
    Vec<CosterKronigRecord>,
    Vec<PhotoabsorptionRecord>,
    Vec<ScatteringRecord>,
);

/// Parse `elam.dat` — the main Elam, Ravel, Sieber data file.
///
/// Uses a state-machine parser matching the Python `add_Elam` function.
/// Returns (xray_levels, xray_transitions, coster_kronig, photoabsorption, scattering).
pub fn parse_elam(path: &Path) -> ElamParseResult {
    let content = std::fs::read_to_string(path).expect("failed to read elam.dat");
    let mut lines: Vec<&str> = content.lines().collect();

    // Validate header
    if lines.is_empty() || !lines[0].contains("Elam, Ravel, Sieber") {
        panic!("Source file not recognized as Elam data");
    }

    // Skip comment lines at the start
    while !lines.is_empty() && lines[0].starts_with('/') {
        lines.remove(0);
    }

    let mut xray_levels = Vec::new();
    let mut xray_transitions = Vec::new();
    let mut coster_kronig = Vec::new();
    let mut photoabsorption = Vec::new();
    let mut scattering = Vec::new();

    let mut current_element = String::new();
    let mut current_edge = String::new();

    let mut idx = 0;
    while idx < lines.len() {
        let line = lines[idx];

        if line.starts_with("Element") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            // Element sym num mw rho
            current_element = parts[1].to_string();
        } else if line.starts_with("Edge") {
            let parts: Vec<&str> = line.split_whitespace().collect();
            // Edge label energy yield jump
            let label = parts[1];
            let energy: f64 = parts[2].parse().unwrap();
            let yield_: f64 = parts[3].parse().unwrap();
            let jump: f64 = parts[4].parse().unwrap();

            current_edge = label.to_string();

            xray_levels.push(XrayLevelRecord {
                element: current_element.clone(),
                iupac_symbol: label.to_string(),
                absorption_edge: energy,
                fluorescence_yield: yield_,
                jump_ratio: jump,
            });
        } else if line.starts_with("  Lines") {
            // Read indented lines that follow
            idx += 1;
            while idx < lines.len() && lines[idx].starts_with("    ") {
                let parts: Vec<&str> = lines[idx].split_whitespace().collect();
                // iupac siegbahn energy intensity
                if parts.len() >= 4 {
                    let iupac = parts[0];
                    let siegbahn = parts[1];
                    let energy: f64 = parts[2].parse().unwrap();
                    let intensity: f64 = parts[3].parse().unwrap();

                    // Split iupac on '-' to get initial/final levels
                    let (start, end) = if let Some(pos) = iupac.find('-') {
                        (iupac[..pos].to_string(), iupac[pos + 1..].to_string())
                    } else {
                        (iupac.to_string(), String::new())
                    };

                    xray_transitions.push(XrayTransitionRecord {
                        element: current_element.clone(),
                        iupac_symbol: iupac.to_string(),
                        siegbahn_symbol: siegbahn.to_string(),
                        initial_level: start,
                        final_level: end,
                        emission_energy: energy,
                        intensity,
                    });
                }
                idx += 1;
            }
            continue; // Don't increment idx again
        } else if line.starts_with("  CK ") {
            // Parse CK line: pairs of (final_level, probability)
            let temp: Vec<&str> = line.split_whitespace().skip(1).collect();
            let mut ck_pairs: Vec<(String, f64)> = Vec::new();
            let mut j = 0;
            while j + 1 < temp.len() {
                let final_level = temp[j].to_string();
                let prob: f64 = temp[j + 1].parse().unwrap();
                ck_pairs.push((final_level, prob));
                j += 2;
            }

            // Check if next line is CKtotal
            let mut total_pairs: Vec<(String, f64)> = Vec::new();
            if idx + 1 < lines.len() && lines[idx + 1].starts_with("  CKtotal") {
                idx += 1;
                let temp: Vec<&str> = lines[idx].split_whitespace().skip(1).collect();
                let mut j = 0;
                while j + 1 < temp.len() {
                    let final_level = temp[j].to_string();
                    let total_prob: f64 = temp[j + 1].parse().unwrap();
                    total_pairs.push((final_level, total_prob));
                    j += 2;
                }
            } else {
                // If no CKtotal line, use CK values as total
                total_pairs = ck_pairs.clone();
            }

            for (ck, total) in ck_pairs.iter().zip(total_pairs.iter()) {
                coster_kronig.push(CosterKronigRecord {
                    element: current_element.clone(),
                    initial_level: current_edge.clone(),
                    final_level: ck.0.clone(),
                    transition_probability: ck.1,
                    total_transition_probability: total.1,
                });
            }
        } else if line.starts_with("Photo") {
            let mut energy = Vec::new();
            let mut photo = Vec::new();
            let mut spline = Vec::new();

            idx += 1;
            while idx < lines.len() && lines[idx].starts_with("    ") {
                let parts: Vec<f64> = lines[idx]
                    .split_whitespace()
                    .filter_map(|w| w.parse().ok())
                    .collect();
                if parts.len() >= 3 {
                    energy.push(parts[0]);
                    photo.push(parts[1]);
                    spline.push(parts[2]);
                }
                idx += 1;
            }

            photoabsorption.push(PhotoabsorptionRecord {
                element: current_element.clone(),
                log_energy: energy,
                log_photoabsorption: photo,
                log_photoabsorption_spline: spline,
            });
            continue; // Don't increment idx again
        } else if line.starts_with("Scatter") {
            let mut energy = Vec::new();
            let mut cs = Vec::new();
            let mut css = Vec::new();
            let mut ics = Vec::new();
            let mut icss = Vec::new();

            idx += 1;
            while idx < lines.len() && lines[idx].starts_with("    ") {
                let parts: Vec<f64> = lines[idx]
                    .split_whitespace()
                    .filter_map(|w| w.parse().ok())
                    .collect();
                if parts.len() >= 5 {
                    energy.push(parts[0]);
                    cs.push(parts[1]);
                    css.push(parts[2]);
                    ics.push(parts[3]);
                    icss.push(parts[4]);
                }
                idx += 1;
            }

            scattering.push(ScatteringRecord {
                element: current_element.clone(),
                log_energy: energy,
                log_coherent_scatter: cs,
                log_coherent_scatter_spline: css,
                log_incoherent_scatter: ics,
                log_incoherent_scatter_spline: icss,
            });
            continue; // Don't increment idx again
        }

        idx += 1;
    }

    (
        xray_levels,
        xray_transitions,
        coster_kronig,
        photoabsorption,
        scattering,
    )
}
