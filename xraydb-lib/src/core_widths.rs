use std::collections::HashMap;

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

impl XrayDb {
    /// Returns core hole width(s) for an element in eV.
    ///
    /// Uses the merged corelevel_widths table (Keski-Rahkonen & Krause
    /// updated with Krause & Oliver for K, L1, L2, L3).
    ///
    /// If `edge` is provided, returns the width for that specific edge.
    /// Otherwise returns a map of edge → width for all available edges.
    pub fn core_width(&self, element: &str, edge: Option<&str>) -> Result<HashMap<String, f64>> {
        let z = self.resolve_element(element)?;
        let mut widths = HashMap::new();

        for w in &self.raw().corelevel_widths {
            if w.atomic_number == z {
                if let Some(e) = edge {
                    if w.edge == e {
                        widths.insert(w.edge.clone(), w.width);
                        return Ok(widths);
                    }
                } else {
                    widths.insert(w.edge.clone(), w.width);
                }
            }
        }

        if widths.is_empty() {
            if let Some(e) = edge {
                return Err(XrayDbError::UnknownEdge {
                    element: element.to_string(),
                    edge: e.to_string(),
                });
            }
            return Err(XrayDbError::UnknownElement(element.to_string()));
        }

        Ok(widths)
    }
}
