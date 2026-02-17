use std::collections::HashMap;

use crate::error::{Result, XrayDbError};

const ELEMENTS: &[&str] = &[
    "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au", "B", "Ba", "Be", "Bi", "Bk", "Br", "C", "Ca",
    "Cd", "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu", "Dy", "Er", "Es", "Eu", "F", "Fe", "Fm",
    "Fr", "Ga", "Gd", "Ge", "H", "He", "Hf", "Hg", "Ho", "I", "In", "Ir", "K", "Kr", "La", "Li",
    "Lr", "Lu", "Md", "Mg", "Mn", "Mo", "N", "Na", "Nb", "Nd", "Ne", "Ni", "No", "Np", "O", "Os",
    "P", "Pa", "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu", "Ra", "Rb", "Re", "Rh", "Rn", "Ru", "S",
    "Sb", "Sc", "Se", "Si", "Sm", "Sn", "Sr", "Ta", "Tb", "Tc", "Te", "Th", "Ti", "Tl", "Tm", "U",
    "Unh", "Unp", "Unq", "Uns", "V", "W", "Xe", "Y", "Yb", "Zn", "Zr",
];

fn is_element(sym: &str) -> bool {
    // D is an alias for H
    sym == "D" || ELEMENTS.contains(&sym)
}

fn resolve_element(sym: &str) -> &str {
    if sym == "D" { "H" } else { sym }
}

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

    fn next_token(&mut self) -> std::result::Result<Token, String> {
        if self.pos >= self.chars.len() {
            return Ok(Token::Eos);
        }

        let ch = self.chars[self.pos];

        if ch == '(' {
            self.pos += 1;
            return Ok(Token::LParen);
        }
        if ch == ')' {
            self.pos += 1;
            return Ok(Token::RParen);
        }

        // Number: starts with digit or '.'
        if ch.is_ascii_digit() || ch == '.' {
            return self.read_number();
        }

        // Element name: starts with uppercase letter
        if ch.is_ascii_uppercase() {
            let start = self.pos;
            self.pos += 1;
            while self.pos < self.chars.len() && self.chars[self.pos].is_ascii_lowercase() {
                self.pos += 1;
            }
            let name: String = self.chars[start..self.pos].iter().collect();
            return Ok(Token::Name(name));
        }

        Err(format!(
            "unrecognized character '{}' at position {}",
            ch, self.pos
        ))
    }

    fn read_number(&mut self) -> std::result::Result<Token, String> {
        let start = self.pos;

        // Integer part
        while self.pos < self.chars.len() && self.chars[self.pos].is_ascii_digit() {
            self.pos += 1;
        }

        // Decimal part
        if self.pos < self.chars.len() && self.chars[self.pos] == '.' {
            self.pos += 1;
            while self.pos < self.chars.len() && self.chars[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
        }

        // Exponent part
        if self.pos < self.chars.len()
            && (self.chars[self.pos] == 'e' || self.chars[self.pos] == 'E')
        {
            self.pos += 1;
            if self.pos < self.chars.len()
                && (self.chars[self.pos] == '+' || self.chars[self.pos] == '-')
            {
                self.pos += 1;
            }
            while self.pos < self.chars.len() && self.chars[self.pos].is_ascii_digit() {
                self.pos += 1;
            }
        }

        let s: String = self.chars[start..self.pos].iter().collect();
        s.parse::<f64>()
            .map(Token::Num)
            .map_err(|_| format!("invalid number '{s}'"))
    }
}

/// Parse a chemical formula into a map of element symbol to count.
///
/// Supports nested parentheses, floating-point stoichiometries, and
/// scientific notation.
///
/// # Examples
/// ```
/// let result = xraydb::chemparser::chemparse("H2O").unwrap();
/// assert_eq!(*result.get("H").unwrap(), 2.0);
/// assert_eq!(*result.get("O").unwrap(), 1.0);
/// ```
pub fn chemparse(formula: &str) -> Result<HashMap<String, f64>> {
    // Handle numbers that start with '.' by inserting '0':
    //   "Fe.7Mg.3O" -> "Fe0.7Mg0.3O"
    let formula = preprocess_formula(formula);

    let mut tokenizer = Tokenizer::new(&formula);
    let current = tokenizer
        .next_token()
        .map_err(XrayDbError::InvalidFormula)?;

    let (result, next) = parse_sequence(&mut tokenizer, current)?;

    if next != Token::Eos {
        return Err(XrayDbError::InvalidFormula(format!(
            "unexpected token after formula: {formula}"
        )));
    }

    let mut out = HashMap::new();
    add_to_result(&result, 1.0, &mut out);
    Ok(out)
}

/// Returns true if the formula can be successfully parsed.
pub fn validate_formula(formula: &str) -> bool {
    chemparse(formula).is_ok()
}

fn preprocess_formula(formula: &str) -> String {
    let formula = formula.replace(' ', "");
    let chars: Vec<char> = formula.chars().collect();
    let mut result = String::with_capacity(formula.len() + 10);

    for (i, &ch) in chars.iter().enumerate() {
        if ch == '.' && (i == 0 || !chars[i - 1].is_ascii_digit()) {
            result.push('0');
        }
        result.push(ch);
    }

    result
}

#[derive(Debug)]
enum FormulaNode {
    Element(String),
    Sequence(Vec<(FormulaNode, f64)>),
}

fn add_to_result(node: &FormulaNode, weight: f64, result: &mut HashMap<String, f64>) {
    match node {
        FormulaNode::Element(sym) => {
            *result.entry(sym.clone()).or_insert(0.0) += weight;
        }
        FormulaNode::Sequence(items) => {
            for (child, count) in items {
                add_to_result(child, weight * count, result);
            }
        }
    }
}

fn parse_sequence(tokenizer: &mut Tokenizer, mut current: Token) -> Result<(FormulaNode, Token)> {
    let mut items: Vec<(FormulaNode, f64)> = Vec::new();

    loop {
        match &current {
            Token::LParen => {
                // Parse nested sequence
                current = tokenizer
                    .next_token()
                    .map_err(XrayDbError::InvalidFormula)?;
                let (inner, next) = parse_sequence(tokenizer, current)?;
                if next != Token::RParen {
                    return Err(XrayDbError::InvalidFormula(
                        "expected closing parenthesis".to_string(),
                    ));
                }
                current = tokenizer
                    .next_token()
                    .map_err(XrayDbError::InvalidFormula)?;

                // Optional count after ')'
                let count = if let Token::Num(n) = current {
                    current = tokenizer
                        .next_token()
                        .map_err(XrayDbError::InvalidFormula)?;
                    n
                } else {
                    1.0
                };
                items.push((inner, count));
            }
            Token::Name(name) => {
                let sym = name.clone();
                if !is_element(&sym) {
                    return Err(XrayDbError::InvalidFormula(format!(
                        "'{sym}' is not an element symbol"
                    )));
                }
                let resolved = resolve_element(&sym).to_string();
                current = tokenizer
                    .next_token()
                    .map_err(XrayDbError::InvalidFormula)?;

                // Optional count after element
                let count = if let Token::Num(n) = current {
                    current = tokenizer
                        .next_token()
                        .map_err(XrayDbError::InvalidFormula)?;
                    n
                } else {
                    1.0
                };
                items.push((FormulaNode::Element(resolved), count));
            }
            _ => break,
        }
    }

    Ok((FormulaNode::Sequence(items), current))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_water() {
        let result = chemparse("H2O").unwrap();
        assert_eq!(result["H"], 2.0);
        assert_eq!(result["O"], 1.0);
    }

    #[test]
    fn test_nested_parens() {
        let result = chemparse("Mn(SO4)2(H2O)7").unwrap();
        assert_eq!(result["Mn"], 1.0);
        assert_eq!(result["S"], 2.0);
        assert_eq!(result["O"], 15.0);
        assert_eq!(result["H"], 14.0);
    }

    #[test]
    fn test_scientific_notation() {
        let result = chemparse("Zn1.e-5Fe3O4").unwrap();
        assert!((result["Zn"] - 1e-5).abs() < 1e-10);
        assert_eq!(result["Fe"], 3.0);
        assert_eq!(result["O"], 4.0);
    }

    #[test]
    fn test_co_vs_co() {
        // CO = carbon monoxide
        let co = chemparse("CO").unwrap();
        assert_eq!(co["C"], 1.0);
        assert_eq!(co["O"], 1.0);

        // Co = cobalt
        let cobalt = chemparse("Co").unwrap();
        assert_eq!(cobalt["Co"], 1.0);
    }

    #[test]
    fn test_decimal_stoichiometry() {
        let result = chemparse("Fe0.7Mg0.3O").unwrap();
        assert!((result["Fe"] - 0.7).abs() < 1e-10);
        assert!((result["Mg"] - 0.3).abs() < 1e-10);
        assert_eq!(result["O"], 1.0);
    }

    #[test]
    fn test_decimal_starting_with_dot() {
        let result = chemparse("Fe.7Mg.3O").unwrap();
        assert!((result["Fe"] - 0.7).abs() < 1e-10);
        assert!((result["Mg"] - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_formula() {
        assert!(chemparse("co").is_err()); // lowercase
        assert!(chemparse("Xx").is_err()); // not an element
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
        let result = chemparse("D2O").unwrap();
        assert_eq!(result["H"], 2.0); // D maps to H
        assert_eq!(result["O"], 1.0);
    }
}
