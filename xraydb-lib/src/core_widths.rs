//! Core-hole level widths.

use std::collections::HashMap;

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

impl XrayDb {
    /// Core-hole width in eV for a single edge.
    ///
    /// Uses the merged `corelevel_widths` table (Keski-Rahkonen & Krause, updated with
    /// Krause & Oliver for K, L1, L2, L3).
    pub fn core_width_at(&self, element: &str, edge: &str) -> Result<f64> {
        let z = self.resolve_element(element)?;
        let indices = self.width_indices(z);
        for &idx in indices {
            if let Some(w) = self.raw().corelevel_widths.get(idx)
                && w.edge == edge
            {
                return Ok(w.width);
            }
        }
        let available = indices
            .iter()
            .filter_map(|&i| self.raw().corelevel_widths.get(i))
            .map(|w| w.edge.clone())
            .collect();
        Err(XrayDbError::UnknownEdge {
            element: element.to_string(),
            edge: edge.to_string(),
            available,
        })
    }

    /// Core-hole width(s) for an element, in eV.
    ///
    /// With `edge = Some(label)` the map holds that one edge; with `None` it holds every
    /// edge available for the element.
    pub fn core_width(&self, element: &str, edge: Option<&str>) -> Result<HashMap<String, f64>> {
        if let Some(e) = edge {
            let width = self.core_width_at(element, e)?;
            return Ok(HashMap::from([(e.to_string(), width)]));
        }

        let z = self.resolve_element(element)?;
        let widths: HashMap<String, f64> = self
            .width_indices(z)
            .iter()
            .filter_map(|&idx| self.raw().corelevel_widths.get(idx))
            .map(|w| (w.edge.clone(), w.width))
            .collect();

        if widths.is_empty() {
            return Err(XrayDbError::UnknownElement(element.to_string()));
        }
        Ok(widths)
    }
}
