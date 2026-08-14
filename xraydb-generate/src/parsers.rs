use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result, anyhow};

use xraydb_data::{
    ComptonEnergiesRecord, CoreWidthRecord, ElementRecord, IonizationPotentialRecord, VersionRecord,
};

/// Parse a whitespace-delimited field, reporting which file and line failed.
pub fn field<T: std::str::FromStr>(parts: &[&str], idx: usize, what: &str, line: &str) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    let raw = parts
        .get(idx)
        .with_context(|| format!("missing {what} (field {idx}) in line: {line}"))?;
    raw.parse::<T>()
        .map_err(|e| anyhow!("invalid {what} '{raw}' in line '{line}': {e}"))
}

/// Read a file, naming it in the error.
pub fn read(path: &Path) -> Result<String> {
    std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))
}

pub fn parse_version(path: &Path) -> Result<Vec<VersionRecord>> {
    let content = read(path)?;
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
    Ok(records)
}

pub fn parse_elements(path: &Path) -> Result<Vec<ElementRecord>> {
    let content = read(path)?;
    let mut records = Vec::new();
    for line in content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 5 {
            records.push(ElementRecord {
                atomic_number: field(&parts, 0, "atomic number", line)?,
                symbol: parts[1].to_string(),
                name: parts[2].to_string(),
                molar_mass: field(&parts, 3, "molar mass", line)?,
                density: field(&parts, 4, "density", line)?,
            });
        }
    }
    Ok(records)
}

pub fn parse_ionization_potentials(path: &Path) -> Result<Vec<IonizationPotentialRecord>> {
    let content = read(path)?;
    let mut records = Vec::new();
    for line in content.lines() {
        if line.starts_with('#') || line.trim().len() < 2 {
            continue;
        }
        let line = line.trim();
        let mut words: Vec<&str> = line.split_whitespace().collect();
        if words.len() >= 2 {
            let raw = words.pop().unwrap_or_default();
            let potential: f64 = raw
                .parse()
                .map_err(|e| anyhow!("invalid ionization potential '{raw}' in '{line}': {e}"))?;
            let gas = words.join(" ");
            records.push(IonizationPotentialRecord { gas, potential });
        }
    }
    Ok(records)
}

pub fn parse_compton_energies(path: &Path) -> Result<ComptonEnergiesRecord> {
    let content = read(path)?;
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

    Ok(ComptonEnergiesRecord {
        incident,
        xray_90deg,
        xray_mean,
        electron_mean,
    })
}

pub fn parse_corehole_data(
    kk_path: &Path,
    ko_path: &Path,
) -> Result<(
    Vec<CoreWidthRecord>,
    Vec<CoreWidthRecord>,
    Vec<CoreWidthRecord>,
)> {
    // Parse Keski-Rahkonen and Krause
    let kk_content = read(kk_path)?;
    let mut kk_records = Vec::new();
    for line in kk_content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 4 {
            kk_records.push(CoreWidthRecord {
                atomic_number: field(&parts, 0, "atomic number", line)?,
                element: parts[1].to_string(),
                edge: parts[2].to_string(),
                width: field(&parts, 3, "core width", line)?,
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
    let ko_content = read(ko_path)?;
    let mut ko_records = Vec::new();
    for line in ko_content.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 6 {
            let atno: u16 = field(&parts, 0, "atomic number", line)?;
            let elem = parts[1].to_string();
            let kwid: f64 = field(&parts, 2, "K width", line)?;
            let l1wid: f64 = field(&parts, 3, "L1 width", line)?;
            let l2wid: f64 = field(&parts, 4, "L2 width", line)?;
            let l3wid: f64 = field(&parts, 5, "L3 width", line)?;

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

    Ok((kk_records, ko_records, corelevel_vec))
}
