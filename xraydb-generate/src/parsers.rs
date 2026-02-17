use std::collections::HashMap;
use std::path::Path;

use xraydb_data::{
    ComptonEnergiesRecord, CoreWidthRecord, ElementRecord, IonizationPotentialRecord,
    VersionRecord,
};

pub fn parse_version(path: &Path) -> Vec<VersionRecord> {
    let content = std::fs::read_to_string(path).expect("failed to read Version.dat");
    let mut records = Vec::new();
    for line in content.lines() {
        if line.starts_with('#') || line.trim().len() < 3 {
            continue;
        }
        let parts: Vec<&str> = line.splitn(3, "//").collect();
        if parts.len() == 3 {
            records.push(VersionRecord {
                tag: parts[0].trim().to_string(),
                date: parts[1].trim().to_string(),
                notes: parts[2].trim().to_string(),
            });
        }
    }
    records
}

pub fn parse_elements(path: &Path) -> Vec<ElementRecord> {
    let content = std::fs::read_to_string(path).expect("failed to read elemental_data.txt");
    let mut records = Vec::new();
    for line in content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 5 {
            records.push(ElementRecord {
                atomic_number: parts[0].parse().unwrap(),
                symbol: parts[1].to_string(),
                name: parts[2].to_string(),
                molar_mass: parts[3].parse().unwrap(),
                density: parts[4].parse().unwrap(),
            });
        }
    }
    records
}

pub fn parse_ionization_potentials(path: &Path) -> Vec<IonizationPotentialRecord> {
    let content =
        std::fs::read_to_string(path).expect("failed to read ion_chamber_potentials.txt");
    let mut records = Vec::new();
    for line in content.lines() {
        if line.starts_with('#') || line.trim().len() < 2 {
            continue;
        }
        let line = line.trim();
        let mut words: Vec<&str> = line.split_whitespace().collect();
        if words.len() >= 2 {
            let potential: f64 = words.pop().unwrap().parse().unwrap();
            let gas = words.join(" ");
            records.push(IonizationPotentialRecord { gas, potential });
        }
    }
    records
}

pub fn parse_compton_energies(path: &Path) -> ComptonEnergiesRecord {
    let content = std::fs::read_to_string(path).expect("failed to read Compton_energies.txt");
    let mut incident = Vec::new();
    let mut xray_90deg = Vec::new();
    let mut xray_mean = Vec::new();
    let mut electron_mean = Vec::new();

    for line in content.lines() {
        if line.starts_with('#') || line.trim().len() < 2 {
            continue;
        }
        let parts: Vec<f64> = line
            .split_whitespace()
            .filter_map(|w| w.parse().ok())
            .collect();
        if parts.len() >= 4 {
            incident.push(parts[0]);
            xray_90deg.push(parts[1]);
            xray_mean.push(parts[2]);
            electron_mean.push(parts[3]);
        }
    }

    ComptonEnergiesRecord {
        incident,
        xray_90deg,
        xray_mean,
        electron_mean,
    }
}

pub fn parse_corehole_data(
    kk_path: &Path,
    ko_path: &Path,
) -> (
    Vec<CoreWidthRecord>,
    Vec<CoreWidthRecord>,
    Vec<CoreWidthRecord>,
) {
    // Parse Keski-Rahkonen and Krause
    let kk_content = std::fs::read_to_string(kk_path).expect("failed to read KK data");
    let mut kk_records = Vec::new();
    for line in kk_content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 4 {
            kk_records.push(CoreWidthRecord {
                atomic_number: parts[0].parse().unwrap(),
                element: parts[1].to_string(),
                edge: parts[2].to_string(),
                width: parts[3].parse().unwrap(),
            });
        }
    }

    // Start corelevel_widths with KK data
    let mut corelevel: HashMap<(u16, String), CoreWidthRecord> = HashMap::new();
    for r in &kk_records {
        corelevel.insert(
            (r.atomic_number, r.edge.clone()),
            CoreWidthRecord {
                atomic_number: r.atomic_number,
                element: r.element.clone(),
                edge: r.edge.clone(),
                width: r.width,
            },
        );
    }

    // Parse Krause and Oliver
    let ko_content = std::fs::read_to_string(ko_path).expect("failed to read KO data");
    let mut ko_records = Vec::new();
    for line in ko_content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 6 {
            let atno: u16 = parts[0].parse().unwrap();
            let elem = parts[1].to_string();
            let kwid: f64 = parts[2].parse().unwrap();
            let l1wid: f64 = parts[3].parse().unwrap();
            let l2wid: f64 = parts[4].parse().unwrap();
            let l3wid: f64 = parts[5].parse().unwrap();

            for (edge, width) in [("K", kwid), ("L1", l1wid), ("L2", l2wid), ("L3", l3wid)] {
                ko_records.push(CoreWidthRecord {
                    atomic_number: atno,
                    element: elem.clone(),
                    edge: edge.to_string(),
                    width,
                });
                // Update corelevel_widths (KO overrides KK for K, L1, L2, L3)
                corelevel.insert(
                    (atno, edge.to_string()),
                    CoreWidthRecord {
                        atomic_number: atno,
                        element: elem.clone(),
                        edge: edge.to_string(),
                        width,
                    },
                );
            }
        }
    }

    // Sort corelevel_widths by (atomic_number, edge)
    let mut corelevel_vec: Vec<CoreWidthRecord> = corelevel.into_values().collect();
    corelevel_vec.sort_by(|a, b| {
        a.atomic_number
            .cmp(&b.atomic_number)
            .then_with(|| a.edge.cmp(&b.edge))
    });

    (kk_records, ko_records, corelevel_vec)
}
