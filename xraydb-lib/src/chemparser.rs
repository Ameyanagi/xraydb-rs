//! Chemical formula parsing.
//!
//! This is the single formula parser used by the whole crate — both
//! [`XrayDb::material_mu`](crate::XrayDb::material_mu) and
//! [`XrayDb::xray_delta_beta`](crate::XrayDb::xray_delta_beta) go through it, so any
//! formula accepted by one is accepted by the other.
//!
//! Supported notation:
//!
//! | Form | Example |
//! |---|---|
//! | plain | `H2O` |
//! | nested groups | `Mn(SO4)2(H2O)7` |
//! | fractional stoichiometry | `Fe0.7Mg0.3O`, `Fe.7Mg.3O` |
//! | scientific notation | `Zn1.e-5Fe3O4` |
//! | fractional group multipliers | `(N2)0.7808(O2)0.2095` |
//! | deuterium alias | `D2O` (counted as hydrogen) |

use std::ops::Deref;

use crate::db::XrayDb;
use crate::error::{Result, XrayDbError};

/// A parsed formula: element symbols with their (possibly fractional) counts.
///
/// Entries are sorted by symbol, which makes downstream floating-point sums
/// reproducible run to run — a `HashMap` would not be.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Composition(Vec<(String, f64)>);

impl Composition {
    /// Count for one element symbol, or `None` if absent.
    pub fn get(&self, symbol: &str) -> Option<f64> {
        self.0
            .iter()
            .find(|(s, _)| s == symbol)
            .map(|&(_, count)| count)
    }

    /// Iterate over `(symbol, count)` pairs in symbol order.
    pub fn iter(&self) -> std::slice::Iter<'_, (String, f64)> {
        self.0.iter()
    }

    /// Consume into the underlying sorted vector.
    pub fn into_vec(self) -> Vec<(String, f64)> {
        self.0
    }

    /// Number of distinct elements.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// True when the formula contained no elements.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

impl Deref for Composition {
    type Target = [(String, f64)];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> IntoIterator for &'a Composition {
    type Item = &'a (String, f64);
    type IntoIter = std::slice::Iter<'a, (String, f64)>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl IntoIterator for Composition {
    type Item = (String, f64);
    type IntoIter = std::vec::IntoIter<(String, f64)>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

/// Deuterium is written `D` but is chemically hydrogen; the database has no `D` row.
const DEUTERIUM: &str = "D";

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Name(String),
    Num(f64),
    LParen,
    RParen,
    Eos,
}

struct Tokenizer {
    chars: Vec<char>,
    pos: usize,
}

impl Tokenizer {
    fn new(input: &str) -> Self {
        Tokenizer {
            chars: input.chars().collect(),
            pos: 0,
        }
    }

    fn next_token(&mut self) -> core::result::Result<Token, String> {
        let Some(&ch) = self.chars.get(self.pos) else {
            return Ok(Token::Eos);
        };

        if ch == '(' {
            self.pos += 1;
            return Ok(Token::LParen);
        }
        if ch == ')' {
            self.pos += 1;
            return Ok(Token::RParen);
        }
        if ch.is_ascii_digit() || ch == '.' {
            return self.read_number();
        }
        if ch.is_ascii_uppercase() {
            let start = self.pos;
            self.pos += 1;
            while self
                .chars
                .get(self.pos)
                .is_some_and(char::is_ascii_lowercase)
            {
                self.pos += 1;
            }
            let name: String = self.chars[start..self.pos].iter().collect();
            return Ok(Token::Name(name));
        }

        Err(format!(
            "unrecognized character '{ch}' at position {}",
            self.pos
        ))
    }

    fn read_number(&mut self) -> core::result::Result<Token, String> {
        let start = self.pos;

        while self.chars.get(self.pos).is_some_and(char::is_ascii_digit) {
            self.pos += 1;
        }
        if self.chars.get(self.pos) == Some(&'.') {
            self.pos += 1;
            while self.chars.get(self.pos).is_some_and(char::is_ascii_digit) {
                self.pos += 1;
            }
        }
        if matches!(self.chars.get(self.pos), Some('e' | 'E')) {
            self.pos += 1;
            if matches!(self.chars.get(self.pos), Some('+' | '-')) {
                self.pos += 1;
            }
            while self.chars.get(self.pos).is_some_and(char::is_ascii_digit) {
                self.pos += 1;
            }
        }

        let s: String = self.chars[start..self.pos].iter().collect();
        s.parse::<f64>()
            .map(Token::Num)
            .map_err(|_| format!("invalid number '{s}'"))
    }
}

impl XrayDb {
    /// Parse a chemical formula, validating element symbols against the database.
    ///
    /// ```
    /// # use xraydb::XrayDb;
    /// let db = XrayDb::try_new()?;
    /// let c = db.parse_formula("Mn(SO4)2(H2O)7")?;
    /// assert_eq!(c.get("O"), Some(15.0));
    /// assert_eq!(c.get("H"), Some(14.0));
    /// # Ok::<(), xraydb::XrayDbError>(())
    /// ```
    pub fn parse_formula(&self, formula: &str) -> Result<Composition> {
        let prepared = preprocess_formula(formula);

        let mut tokenizer = Tokenizer::new(&prepared);
        let current = tokenizer
            .next_token()
            .map_err(|e| XrayDbError::invalid_formula(formula, e))?;

        let (node, next) = self.parse_sequence(formula, &mut tokenizer, current)?;
        if next != Token::Eos {
            return Err(XrayDbError::invalid_formula(
                formula,
                "unexpected trailing token",
            ));
        }

        let mut accumulated: Vec<(String, f64)> = Vec::new();
        accumulate(&node, 1.0, &mut accumulated);
        accumulated.sort_by(|a, b| a.0.cmp(&b.0));

        // Merge duplicate symbols (e.g. the two carbons in "C6H4(CH3)2").
        let mut merged: Vec<(String, f64)> = Vec::with_capacity(accumulated.len());
        for (symbol, count) in accumulated {
            match merged.last_mut() {
                Some((prev, total)) if *prev == symbol => *total += count,
                _ => merged.push((symbol, count)),
            }
        }

        if merged.is_empty() {
            return Err(XrayDbError::invalid_formula(formula, "no elements found"));
        }
        Ok(Composition(merged))
    }

    /// True if [`XrayDb::parse_formula`] would succeed.
    pub fn validate_formula(&self, formula: &str) -> bool {
        self.parse_formula(formula).is_ok()
    }

    fn resolve_symbol(&self, formula: &str, name: &str) -> Result<String> {
        if name == DEUTERIUM {
            return Ok("H".to_string());
        }
        if self.contains_symbol(name) {
            return Ok(name.to_string());
        }
        Err(XrayDbError::invalid_formula(
            formula,
            format!("'{name}' is not an element symbol"),
        ))
    }

    fn parse_sequence(
        &self,
        formula: &str,
        tokenizer: &mut Tokenizer,
        mut current: Token,
    ) -> Result<(Node, Token)> {
        let mut items: Vec<(Node, f64)> = Vec::new();
        let bad = |e: String| XrayDbError::invalid_formula(formula, e);

        loop {
            match &current {
                Token::LParen => {
                    current = tokenizer.next_token().map_err(bad)?;
                    let (inner, next) = self.parse_sequence(formula, tokenizer, current)?;
                    if next != Token::RParen {
                        return Err(XrayDbError::invalid_formula(
                            formula,
                            "expected closing parenthesis",
                        ));
                    }
                    current = tokenizer.next_token().map_err(bad)?;

                    let count = if let Token::Num(n) = current {
                        current = tokenizer.next_token().map_err(bad)?;
                        n
                    } else {
                        1.0
                    };
                    items.push((inner, count));
                }
                Token::Name(name) => {
                    let symbol = self.resolve_symbol(formula, name)?;
                    current = tokenizer.next_token().map_err(bad)?;

                    let count = if let Token::Num(n) = current {
                        current = tokenizer.next_token().map_err(bad)?;
                        n
                    } else {
                        1.0
                    };
                    items.push((Node::Element(symbol), count));
                }
                _ => break,
            }
        }

        Ok((Node::Sequence(items), current))
    }
}

/// Parse a chemical formula using the embedded database for symbol validation.
///
/// Equivalent to [`XrayDb::parse_formula`]; provided as a free function for one-off use.
/// The database is a lazily-initialised static, so this is cheap after first call.
///
/// Uses [`XrayDb::current`], so it works whether the data is compiled in or was
/// supplied at runtime.
///
/// ```
/// # #[cfg(feature = "embedded-data")] {
/// let water = xraydb::chemparser::chemparse("H2O")?;
/// assert_eq!(water.get("H"), Some(2.0));
/// assert_eq!(water.get("O"), Some(1.0));
/// # }
/// # Ok::<(), xraydb::XrayDbError>(())
/// ```
pub fn chemparse(formula: &str) -> Result<Composition> {
    XrayDb::current()?.parse_formula(formula)
}

/// True if [`chemparse`] would succeed.
pub fn validate_formula(formula: &str) -> bool {
    chemparse(formula).is_ok()
}

/// Strip spaces and give bare leading decimal points an explicit zero,
/// so `Fe.7Mg.3O` tokenizes the same as `Fe0.7Mg0.3O`.
fn preprocess_formula(formula: &str) -> String {
    let chars: Vec<char> = formula.chars().filter(|c| !c.is_whitespace()).collect();
    let mut result = String::with_capacity(chars.len() + 8);

    for (i, &ch) in chars.iter().enumerate() {
        if ch == '.' && (i == 0 || !chars[i - 1].is_ascii_digit()) {
            result.push('0');
        }
        result.push(ch);
    }
    result
}

#[derive(Debug)]
enum Node {
    Element(String),
    Sequence(Vec<(Node, f64)>),
}

fn accumulate(node: &Node, weight: f64, out: &mut Vec<(String, f64)>) {
    match node {
        Node::Element(sym) => out.push((sym.clone(), weight)),
        Node::Sequence(items) => {
            for (child, count) in items {
                accumulate(child, weight * count, out);
            }
        }
    }
}

#[cfg(all(test, feature = "embedded-data"))]
mod tests {
    use super::*;

    fn get(formula: &str, symbol: &str) -> f64 {
        chemparse(formula).unwrap().get(symbol).unwrap()
    }

    #[test]
    fn test_water() {
        assert_eq!(get("H2O", "H"), 2.0);
        assert_eq!(get("H2O", "O"), 1.0);
    }

    #[test]
    fn test_nested_parens() {
        assert_eq!(get("Mn(SO4)2(H2O)7", "Mn"), 1.0);
        assert_eq!(get("Mn(SO4)2(H2O)7", "S"), 2.0);
        assert_eq!(get("Mn(SO4)2(H2O)7", "O"), 15.0);
        assert_eq!(get("Mn(SO4)2(H2O)7", "H"), 14.0);
    }

    #[test]
    fn test_scientific_notation() {
        assert!((get("Zn1.e-5Fe3O4", "Zn") - 1e-5).abs() < 1e-10);
        assert_eq!(get("Zn1.e-5Fe3O4", "Fe"), 3.0);
        assert_eq!(get("Zn1.e-5Fe3O4", "O"), 4.0);
    }

    #[test]
    fn test_co_vs_co() {
        assert_eq!(get("CO", "C"), 1.0);
        assert_eq!(get("CO", "O"), 1.0);
        assert_eq!(get("Co", "Co"), 1.0);
    }

    #[test]
    fn test_decimal_stoichiometry() {
        assert!((get("Fe0.7Mg0.3O", "Fe") - 0.7).abs() < 1e-10);
        assert!((get("Fe0.7Mg0.3O", "Mg") - 0.3).abs() < 1e-10);
        assert_eq!(get("Fe0.7Mg0.3O", "O"), 1.0);
    }

    #[test]
    fn test_decimal_starting_with_dot() {
        assert!((get("Fe.7Mg.3O", "Fe") - 0.7).abs() < 1e-10);
        assert!((get("Fe.7Mg.3O", "Mg") - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_formula() {
        assert!(chemparse("co").is_err()); // lowercase
        assert!(chemparse("Xx").is_err()); // not an element
        assert!(chemparse("H2(O").is_err()); // unbalanced
        assert!(chemparse("H2O)").is_err()); // unbalanced
        assert!(chemparse("").is_err()); // empty
        assert!(chemparse("2H").is_err()); // leading count
    }

    #[test]
    fn test_validate() {
        assert!(validate_formula("H2O"));
        assert!(validate_formula("Mn(SO4)2(H2O)7"));
        assert!(!validate_formula("co"));
        assert!(!validate_formula("Xx"));
    }

    #[test]
    fn test_deuterium() {
        assert_eq!(get("D2O", "H"), 2.0);
        assert_eq!(get("D2O", "O"), 1.0);
    }

    #[test]
    fn duplicate_symbols_merge() {
        // Two separate carbon terms must sum, not overwrite.
        assert_eq!(get("C6H4(CH3)2", "C"), 8.0);
        assert_eq!(get("C6H4(CH3)2", "H"), 10.0);
    }

    #[test]
    fn composition_is_sorted_for_reproducible_sums() {
        let c = chemparse("Mn(SO4)2(H2O)7").unwrap();
        let symbols: Vec<&str> = c.iter().map(|(s, _)| s.as_str()).collect();
        let mut sorted = symbols.clone();
        sorted.sort_unstable();
        assert_eq!(symbols, sorted);
    }

    #[test]
    fn the_shipped_materials_table_all_parses() {
        // Regression guard: the crate must be able to parse its own materials data.
        let db = XrayDb::new();
        for material in db.materials() {
            assert!(
                db.validate_formula(material.formula),
                "material '{}' has unparseable formula '{}'",
                material.name,
                material.formula
            );
        }
    }

    #[test]
    fn whitespace_is_ignored() {
        assert_eq!(get("H2 O", "O"), 1.0);
        assert_eq!(get(" Fe 2 O 3 ", "Fe"), 2.0);
    }
}
