use std::collections::HashMap;

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

/// X-ray absorption edge data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct XrayEdge {
    pub energy: f64,
    pub fluorescence_yield: f64,
    pub jump_ratio: f64,
}

/// X-ray emission line data.
#[derive(Debug, Clone, PartialEq)]
pub struct XrayLine {
    pub energy: f64,
    pub intensity: f64,
    pub initial_level: String,
    pub final_level: String,
}

impl XrayDb {
    /// Returns a map of all X-ray absorption edges for an element.
    ///
    /// Keys are IUPAC edge labels (K, L1, L2, L3, M1, ...).
    pub fn xray_edges(&self, element: &str) -> Result<HashMap<String, XrayEdge>> {
        let sym = self.symbol(element)?;
        let mut edges = HashMap::new();
        for level in &self.raw().xray_levels {
            if level.element == sym {
                edges.insert(
                    level.iupac_symbol.clone(),
                    XrayEdge {
                        energy: level.absorption_edge,
                        fluorescence_yield: level.fluorescence_yield,
                        jump_ratio: level.jump_ratio,
                    },
                );
            }
        }
        Ok(edges)
    }

    /// Returns data for a specific X-ray edge.
    pub fn xray_edge(&self, element: &str, edge: &str) -> Result<XrayEdge> {
        let sym = self.symbol(element)?;
        self.raw()
            .xray_levels
            .iter()
            .find(|l| l.element == sym && l.iupac_symbol == edge)
            .map(|l| XrayEdge {
                energy: l.absorption_edge,
                fluorescence_yield: l.fluorescence_yield,
                jump_ratio: l.jump_ratio,
            })
            .ok_or_else(|| XrayDbError::UnknownEdge {
                element: element.to_string(),
                edge: edge.to_string(),
            })
    }

    /// Returns a map of X-ray emission lines for an element.
    ///
    /// Keys are Siegbahn notation (Ka1, Ka2, Kb1, La1, ...).
    /// If `initial_level` is provided, returns only lines from that level.
    /// If `excitation_energy` is provided, returns only lines with excitation
    /// energy below that value.
    pub fn xray_lines(
        &self,
        element: &str,
        initial_level: Option<&str>,
        excitation_energy: Option<f64>,
    ) -> Result<HashMap<String, XrayLine>> {
        let sym = self.symbol(element)?;
        let mut lines = HashMap::new();

        for trans in &self.raw().xray_transitions {
            if trans.element != sym {
                continue;
            }
            if let Some(level) = initial_level {
                if trans.initial_level != level {
                    continue;
                }
            }
            if let Some(max_energy) = excitation_energy {
                // Find the edge energy for this transition's initial level
                if let Some(level) = self
                    .raw()
                    .xray_levels
                    .iter()
                    .find(|l| l.element == sym && l.iupac_symbol == trans.initial_level)
                {
                    if level.absorption_edge > max_energy {
                        continue;
                    }
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

    /// Guess the element and edge from an X-ray energy.
    ///
    /// Returns (element_symbol, edge_label) for the closest match.
    pub fn guess_edge(
        &self,
        energy: f64,
        edges: Option<&[&str]>,
    ) -> Option<(String, String)> {
        let default_edges = ["K", "L3", "L2", "L1", "M5"];
        let edge_filter = edges.unwrap_or(&default_edges);

        let mut best: Option<(String, String, f64)> = None;

        for level in &self.raw().xray_levels {
            if level.absorption_edge <= 0.0 {
                continue;
            }
            if !edge_filter.contains(&level.iupac_symbol.as_str()) {
                continue;
            }
            let diff = (level.absorption_edge - energy).abs();
            if best.is_none() || diff < best.as_ref().unwrap().2 {
                best = Some((level.element.clone(), level.iupac_symbol.clone(), diff));
            }
        }

        best.map(|(elem, edge, _)| (elem, edge))
    }

    /// Returns the ionization potential for a gas (in eV per ion pair).
    pub fn ionization_potential(&self, gas: &str) -> Result<f64> {
        let gas_lower = gas.to_lowercase();
        self.raw()
            .ionization_potentials
            .iter()
            .find(|ip| ip.gas.to_lowercase() == gas_lower)
            .map(|ip| ip.potential)
            .ok_or_else(|| XrayDbError::UnknownGas(gas.to_string()))
    }
}
