//! Absorption edges, emission lines, and edge identification.

use std::collections::HashMap;

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

/// One X-ray absorption edge.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XrayEdge {
    /// Absorption edge energy in eV.
    pub energy: f64,
    /// Fluorescence yield, 0–1.
    pub fluorescence_yield: f64,
    /// Absorption jump ratio across the edge.
    pub jump_ratio: f64,
}

/// One X-ray emission line.
#[derive(Debug, Clone, PartialEq)]
pub struct XrayLine {
    /// Emission energy in eV.
    pub energy: f64,
    /// Relative intensity within the initial level's lines.
    pub intensity: f64,
    /// IUPAC label of the initial (core-hole) level, e.g. `"K"`.
    pub initial_level: String,
    /// IUPAC label of the final level, e.g. `"L3"`.
    pub final_level: String,
}

/// A candidate match from [`XrayDb::guess_edge`].
#[derive(Debug, Clone, PartialEq)]
pub struct EdgeGuess {
    /// Element symbol.
    pub element: String,
    /// IUPAC edge label.
    pub edge: String,
    /// Tabulated edge energy in eV.
    pub energy: f64,
    /// Signed difference `energy - query` in eV.
    pub difference: f64,
}

/// Emission lines from one initial level, strongest first.
pub type LevelLines = (String, Vec<(String, XrayLine)>);

/// Edge labels searched by [`XrayDb::guess_edge`] when none are specified.
pub const DEFAULT_GUESS_EDGES: &[&str] = &["K", "L3", "L2", "L1", "M5"];

impl XrayDb {
    /// All X-ray absorption edges for an element, keyed by IUPAC label
    /// (`K`, `L1`, `L2`, `L3`, `M1`, …).
    pub fn xray_edges(&self, element: &str) -> Result<HashMap<String, XrayEdge>> {
        let sym = self.symbol(element)?;
        Ok(self
            .level_indices(sym)
            .iter()
            .filter_map(|&idx| self.raw().xray_levels.get(idx))
            .map(|l| {
                (
                    l.iupac_symbol.clone(),
                    XrayEdge {
                        energy: l.absorption_edge,
                        fluorescence_yield: l.fluorescence_yield,
                        jump_ratio: l.jump_ratio,
                    },
                )
            })
            .collect())
    }

    /// A single X-ray absorption edge.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let k = db.xray_edge("Fe", "K")?;
    /// assert_eq!(k.energy, 7112.0);
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn xray_edge(&self, element: &str, edge: &str) -> Result<XrayEdge> {
        let sym = self.symbol(element)?;
        self.level_by_edge(sym, edge)
            .map(|l| XrayEdge {
                energy: l.absorption_edge,
                fluorescence_yield: l.fluorescence_yield,
                jump_ratio: l.jump_ratio,
            })
            .ok_or_else(|| XrayDbError::UnknownEdge {
                element: element.to_string(),
                edge: edge.to_string(),
                available: self.edge_labels(sym),
            })
    }

    /// IUPAC edge labels available for an already-resolved symbol.
    fn edge_labels(&self, sym: &str) -> Vec<String> {
        self.level_indices(sym)
            .iter()
            .filter_map(|&idx| self.raw().xray_levels.get(idx))
            .map(|l| l.iupac_symbol.clone())
            .collect()
    }

    /// X-ray emission lines for an element, keyed by Siegbahn label (`Ka1`, `Kb1`, `La1`, …).
    ///
    /// * `initial_level` — restrict to lines from one core-hole level.
    /// * `excitation_energy` — drop lines whose initial level cannot be reached at that
    ///   incident energy.
    pub fn xray_lines(
        &self,
        element: &str,
        initial_level: Option<&str>,
        excitation_energy: Option<f64>,
    ) -> Result<HashMap<String, XrayLine>> {
        let sym = self.symbol(element)?;

        // Resolve each initial level's edge energy once, instead of re-scanning the
        // 1400-row level table for every one of the ~1800 transitions.
        let excitable: Option<HashMap<&str, bool>> = excitation_energy.map(|max_energy| {
            self.level_indices(sym)
                .iter()
                .filter_map(|&idx| self.raw().xray_levels.get(idx))
                .map(|l| (l.iupac_symbol.as_str(), l.absorption_edge <= max_energy))
                .collect()
        });

        let mut lines = HashMap::new();
        for &idx in self.transition_indices(sym) {
            let Some(trans) = self.raw().xray_transitions.get(idx) else {
                continue;
            };
            // Not a let-chain: those need Rust 1.88, and the MSRV is lower on purpose.
            #[allow(clippy::collapsible_if)]
            if let Some(level) = initial_level {
                if trans.initial_level != level {
                    continue;
                }
            }
            if let Some(ref reachable) = excitable {
                // A transition from a level with no tabulated edge is kept, matching
                // the previous behaviour of "no matching level found -> no filtering".
                if reachable.get(trans.initial_level.as_str()) == Some(&false) {
                    continue;
                }
            }

            lines.insert(
                trans.siegbahn_symbol.clone(),
                XrayLine {
                    energy: trans.emission_energy,
                    intensity: trans.intensity,
                    initial_level: trans.initial_level.clone(),
                    final_level: trans.final_level.clone(),
                },
            );
        }
        Ok(lines)
    }

    /// Emission lines sorted by descending intensity.
    ///
    /// Convenience over [`XrayDb::xray_lines`], whose `HashMap` has no useful order.
    ///
    /// # Comparing intensities across levels
    ///
    /// Elam intensities are normalized **within each initial level**, not across the
    /// whole element. Sorting globally therefore ranks a strong L line above Kα1 even
    /// though Kα1 dominates a real fluorescence spectrum. Pass `initial_level` to
    /// compare like with like, or use [`XrayDb::xray_lines_by_level`] to keep the
    /// grouping.
    pub fn xray_lines_by_intensity(
        &self,
        element: &str,
        initial_level: Option<&str>,
        excitation_energy: Option<f64>,
    ) -> Result<Vec<(String, XrayLine)>> {
        let mut lines: Vec<(String, XrayLine)> = self
            .xray_lines(element, initial_level, excitation_energy)?
            .into_iter()
            .collect();
        lines.sort_by(|a, b| {
            b.1.intensity
                .total_cmp(&a.1.intensity)
                .then_with(|| a.0.cmp(&b.0))
        });
        Ok(lines)
    }

    /// Emission lines grouped by initial level, each group sorted by descending
    /// intensity.
    ///
    /// Groups are ordered K, L1, L2, L3, M1…, matching the usual spectroscopic
    /// presentation. Because Elam intensities are relative within a level, this is
    /// the honest way to rank them; see [`XrayDb::xray_lines_by_intensity`].
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let groups = db.xray_lines_by_level("Fe", None)?;
    /// assert_eq!(groups.first().map(|(level, _)| level.as_str()), Some("K"));
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn xray_lines_by_level(
        &self,
        element: &str,
        excitation_energy: Option<f64>,
    ) -> Result<Vec<LevelLines>> {
        let lines = self.xray_lines_by_intensity(element, None, excitation_energy)?;

        let mut groups: Vec<LevelLines> = Vec::new();
        for (label, line) in lines {
            let level = line.initial_level.clone();
            match groups.iter_mut().find(|(l, _)| *l == level) {
                Some((_, items)) => items.push((label, line)),
                None => groups.push((level, vec![(label, line)])),
            }
        }
        groups.sort_by(|a, b| {
            level_order(&a.0)
                .cmp(&level_order(&b.0))
                .then(a.0.cmp(&b.0))
        });
        Ok(groups)
    }

    /// Identify the absorption edge closest to a given energy.
    ///
    /// * `edges` — labels to consider; defaults to [`DEFAULT_GUESS_EDGES`].
    /// * `max_difference` — reject matches further than this many eV away. Pass `None`
    ///   to accept the nearest edge at any distance.
    ///
    /// Returns `None` when nothing is within tolerance.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let guess = db.guess_edge(7100.0, None, Some(100.0)).expect("Fe K is nearby");
    /// assert_eq!(guess.element, "Fe");
    /// assert_eq!(guess.edge, "K");
    /// assert_eq!(guess.difference, 12.0); // Fe K is at 7112 eV
    ///
    /// // 200 keV is above every tabulated K edge; without a tolerance the nearest
    /// // match is 65 keV away, which is not a useful answer.
    /// assert!(db.guess_edge(200_000.0, None, None).is_some());
    /// assert!(db.guess_edge(200_000.0, None, Some(500.0)).is_none());
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn guess_edge(
        &self,
        energy: f64,
        edges: Option<&[&str]>,
        max_difference: Option<f64>,
    ) -> Option<EdgeGuess> {
        let edge_filter = edges.unwrap_or(DEFAULT_GUESS_EDGES);
        let levels = &self.raw().xray_levels;
        let sorted = self.levels_sorted_by_energy();

        // Binary search to the insertion point, then walk outward. The first candidate
        // found in each direction bounds the search: anything further out is worse.
        let start = sorted.partition_point(|&i| levels[i].absorption_edge < energy);

        let mut best: Option<EdgeGuess> = None;
        let consider = |idx: usize, best: &mut Option<EdgeGuess>| {
            let Some(level) = levels.get(idx) else {
                return;
            };
            if level.absorption_edge <= 0.0 || !edge_filter.contains(&level.iupac_symbol.as_str()) {
                return;
            }
            let difference = level.absorption_edge - energy;
            if best
                .as_ref()
                .is_none_or(|b| difference.abs() < b.difference.abs())
            {
                *best = Some(EdgeGuess {
                    element: level.element.clone(),
                    edge: level.iupac_symbol.clone(),
                    energy: level.absorption_edge,
                    difference,
                });
            }
        };

        // Upward from the insertion point.
        for &idx in &sorted[start..] {
            // Not a let-chain: those need Rust 1.88, and the MSRV is lower on purpose.
            #[allow(clippy::collapsible_if)]
            if let Some(ref b) = best {
                if levels[idx].absorption_edge - energy > b.difference.abs() {
                    break;
                }
            }
            consider(idx, &mut best);
        }
        // Downward from the insertion point.
        for &idx in sorted[..start].iter().rev() {
            // Not a let-chain: those need Rust 1.88, and the MSRV is lower on purpose.
            #[allow(clippy::collapsible_if)]
            if let Some(ref b) = best {
                if energy - levels[idx].absorption_edge > b.difference.abs() {
                    break;
                }
            }
            consider(idx, &mut best);
        }

        match (best, max_difference) {
            (Some(b), Some(tol)) if b.difference.abs() > tol => None,
            (best, _) => best,
        }
    }

    /// Ionization potential for a gas, in eV per ion pair.
    pub fn ionization_potential(&self, gas: &str) -> Result<f64> {
        self.gas_potential(gas)
            .ok_or_else(|| XrayDbError::UnknownGas {
                gas: gas.to_string(),
                available: self.known_gases(),
            })
    }
}

/// Sort key placing shells in spectroscopic order: K, then L, M, N, ...
fn level_order(level: &str) -> (u8, u8) {
    let shell = match level.chars().next() {
        Some('K') => 0,
        Some('L') => 1,
        Some('M') => 2,
        Some('N') => 3,
        Some('O') => 4,
        Some('P') => 5,
        _ => 9,
    };
    let index = level[1..].parse::<u8>().unwrap_or(0);
    (shell, index)
}
