//! Output formatting: aligned text tables, CSV, and JSON.

use std::io::{IsTerminal, Write};

use anyhow::Result;
use clap::ValueEnum;
use serde::Serialize;

/// How to render results.
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, Default)]
pub enum Format {
    /// Aligned plain-text columns, suitable for reading and for `awk`.
    #[default]
    Text,
    /// JSON, for scripting.
    Json,
    /// Comma-separated values, for plotting tools.
    Csv,
}

/// A simple column-aligned table.
pub struct Table {
    headers: Vec<String>,
    rows: Vec<Vec<String>>,
    /// Columns rendered right-aligned (numeric columns).
    numeric: Vec<bool>,
}

impl Table {
    /// Build a table with the given column headers.
    pub fn new(headers: &[&str]) -> Self {
        Table {
            headers: headers.iter().map(|h| (*h).to_string()).collect(),
            rows: Vec::new(),
            numeric: vec![true; headers.len()],
        }
    }

    /// Append a row. Cells are compared against the header count only loosely, so a
    /// short row simply renders with blanks.
    pub fn row(&mut self, cells: Vec<String>) {
        for (i, cell) in cells.iter().enumerate() {
            if i < self.numeric.len() && cell.parse::<f64>().is_err() && !cell.is_empty() {
                self.numeric[i] = false;
            }
        }
        self.rows.push(cells);
    }

    /// Render to `out` in the requested format.
    pub fn write(&self, out: &mut impl Write, format: Format) -> Result<()> {
        match format {
            Format::Csv => {
                writeln!(out, "{}", self.headers.join(","))?;
                for row in &self.rows {
                    let cells: Vec<String> = row.iter().map(|c| escape_csv(c)).collect();
                    writeln!(out, "{}", cells.join(","))?;
                }
            }
            Format::Json => {
                let objects: Vec<serde_json::Map<String, serde_json::Value>> = self
                    .rows
                    .iter()
                    .map(|row| {
                        self.headers
                            .iter()
                            .zip(row.iter())
                            .map(|(h, c)| (h.clone(), json_cell(c)))
                            .collect()
                    })
                    .collect();
                writeln!(out, "{}", serde_json::to_string_pretty(&objects)?)?;
            }
            Format::Text => {
                let mut widths: Vec<usize> =
                    self.headers.iter().map(|h| h.chars().count()).collect();
                for row in &self.rows {
                    for (i, cell) in row.iter().enumerate() {
                        if i < widths.len() {
                            widths[i] = widths[i].max(cell.chars().count());
                        }
                    }
                }

                let line: Vec<String> = self
                    .headers
                    .iter()
                    .enumerate()
                    .map(|(i, h)| pad(h, widths[i], self.numeric[i]))
                    .collect();
                writeln!(out, "{}", bold(line.join("  ").trim_end()))?;

                let rule: Vec<String> = widths.iter().map(|w| "-".repeat(*w)).collect();
                writeln!(out, "{}", rule.join("  "))?;

                for row in &self.rows {
                    let line: Vec<String> = row
                        .iter()
                        .enumerate()
                        .map(|(i, c)| pad(c, widths.get(i).copied().unwrap_or(0), self.numeric[i]))
                        .collect();
                    writeln!(out, "{}", line.join("  ").trim_end())?;
                }
            }
        }
        Ok(())
    }
}

fn pad(text: &str, width: usize, right_align: bool) -> String {
    let len = text.chars().count();
    let fill = width.saturating_sub(len);
    if right_align {
        format!("{}{}", " ".repeat(fill), text)
    } else {
        format!("{}{}", text, " ".repeat(fill))
    }
}

fn escape_csv(cell: &str) -> String {
    if cell.contains([',', '"', '\n']) {
        format!("\"{}\"", cell.replace('"', "\"\""))
    } else {
        cell.to_string()
    }
}

fn json_cell(cell: &str) -> serde_json::Value {
    match cell.parse::<f64>() {
        Ok(v) if v.is_finite() => serde_json::Number::from_f64(v).map_or_else(
            || serde_json::Value::String(cell.to_string()),
            serde_json::Value::Number,
        ),
        _ => serde_json::Value::String(cell.to_string()),
    }
}

/// Write a serializable record as JSON, or fall back to a two-column key/value table.
pub fn write_record<T: Serialize>(
    out: &mut impl Write,
    format: Format,
    pairs: &[(&str, String)],
    json_value: &T,
) -> Result<()> {
    match format {
        Format::Json => {
            writeln!(out, "{}", serde_json::to_string_pretty(json_value)?)?;
        }
        Format::Csv => {
            writeln!(out, "field,value")?;
            for (k, v) in pairs {
                writeln!(out, "{},{}", escape_csv(k), escape_csv(v))?;
            }
        }
        Format::Text => {
            let width = pairs
                .iter()
                .map(|(k, _)| k.chars().count())
                .max()
                .unwrap_or(0);
            for (k, v) in pairs {
                writeln!(out, "{}  {}", bold(&pad(k, width, false)), v)?;
            }
        }
    }
    Ok(())
}

/// Format a number with significant figures appropriate to its magnitude.
///
/// Edge energies want integers, µ spans six decades, and delta is ~1e-6, so a blanket
/// `{:.4}` is wrong for at least one of them.
pub fn num(value: f64) -> String {
    if !value.is_finite() {
        return if value.is_nan() {
            "NaN".to_string()
        } else if value > 0.0 {
            "inf".to_string()
        } else {
            "-inf".to_string()
        };
    }
    if value == 0.0 {
        return "0".to_string();
    }
    let mag = value.abs();
    if !(1e-4..1e6).contains(&mag) {
        format!("{value:.5e}")
    } else if value.fract() == 0.0 && mag < 1e9 {
        format!("{value:.0}")
    } else {
        // Six significant figures, trailing zeros trimmed.
        let s = format!("{value:.6}");
        let s = s.trim_end_matches('0').trim_end_matches('.');
        s.to_string()
    }
}

/// True when colour should be used: stdout is a terminal and `NO_COLOR` is unset.
pub fn use_color() -> bool {
    std::io::stdout().is_terminal() && std::env::var_os("NO_COLOR").is_none()
}

/// Wrap `text` in a bold ANSI sequence when colour is enabled.
pub fn bold(text: &str) -> String {
    if use_color() {
        format!("\u{1b}[1m{text}\u{1b}[0m")
    } else {
        text.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn render(table: &Table, format: Format) -> String {
        let mut buf = Vec::new();
        table.write(&mut buf, format).unwrap();
        String::from_utf8(buf).unwrap()
    }

    #[test]
    fn text_table_aligns_columns() {
        let mut t = Table::new(&["edge", "energy"]);
        t.row(vec!["K".into(), "7112".into()]);
        t.row(vec!["L3".into(), "706.8".into()]);
        let out = render(&t, Format::Text);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 4);
        // The numeric column is right-aligned, so both values end at the same column.
        assert!(lines[2].ends_with("7112"), "{:?}", lines[2]);
        assert!(lines[3].ends_with("706.8"), "{:?}", lines[3]);
    }

    #[test]
    fn csv_escapes_separators() {
        let mut t = Table::new(&["name", "formula"]);
        t.row(vec!["a,b".into(), "H2O".into()]);
        let out = render(&t, Format::Csv);
        assert!(out.contains("\"a,b\",H2O"), "{out}");
    }

    #[test]
    fn json_emits_numbers_as_numbers() {
        let mut t = Table::new(&["edge", "energy"]);
        t.row(vec!["K".into(), "7112".into()]);
        let out = render(&t, Format::Json);
        let parsed: serde_json::Value = serde_json::from_str(&out).unwrap();
        assert_eq!(parsed[0]["edge"], "K");
        assert_eq!(parsed[0]["energy"], 7112.0);
    }

    #[test]
    fn num_picks_a_sensible_precision() {
        assert_eq!(num(7112.0), "7112");
        assert_eq!(num(0.0), "0");
        assert_eq!(num(55.845), "55.845");
        assert!(num(4.599e-6).contains('e'), "{}", num(4.599e-6));
        assert!(num(1.5e7).contains('e'));
        assert_eq!(num(f64::INFINITY), "inf");
        assert_eq!(num(f64::NAN), "NaN");
    }
}
