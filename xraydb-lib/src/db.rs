use std::collections::HashMap;
use std::sync::OnceLock;

use xraydb_data::XrayDatabase;

use crate::error::{Result, XrayDbError};

const COMPRESSED_DATA: &[u8] = include_bytes!("../../data/xraydb.bin.zst");

struct InitializedDb {
    data: XrayDatabase,
    symbol_to_z: HashMap<String, u16>,
    name_to_z: HashMap<String, u16>,
}

static DATABASE: OnceLock<InitializedDb> = OnceLock::new();

fn db() -> &'static InitializedDb {
    DATABASE.get_or_init(|| {
        // Decompress with ruzstd
        let mut decoder = ruzstd::decoding::StreamingDecoder::new(COMPRESSED_DATA)
            .expect("failed to create zstd decoder");
        let mut decompressed = Vec::new();
        std::io::Read::read_to_end(&mut decoder, &mut decompressed)
            .expect("failed to decompress data");

        // Deserialize with postcard
        let data: XrayDatabase =
            postcard::from_bytes(&decompressed).expect("failed to deserialize data");

        // Build lookup indices
        let mut symbol_to_z = HashMap::new();
        let mut name_to_z = HashMap::new();
        for elem in &data.elements {
            symbol_to_z.insert(elem.symbol.clone(), elem.atomic_number);
            symbol_to_z.insert(elem.symbol.to_lowercase(), elem.atomic_number);
            name_to_z.insert(elem.name.clone(), elem.atomic_number);
            name_to_z.insert(elem.name.to_lowercase(), elem.atomic_number);
        }

        InitializedDb {
            data,
            symbol_to_z,
            name_to_z,
        }
    })
}

/// The main interface to the X-ray database.
///
/// Cheap to create — holds a reference to statically-allocated data
/// that is decompressed on first use.
pub struct XrayDb {
    db: &'static InitializedDb,
}

impl XrayDb {
    pub fn new() -> Self {
        XrayDb { db: db() }
    }

    /// Access the raw database.
    pub fn raw(&self) -> &XrayDatabase {
        &self.db.data
    }

    /// Resolve an element identifier (symbol, name, or atomic number) to Z.
    pub fn resolve_element(&self, element: &str) -> Result<u16> {
        // Try as atomic number first
        if let Ok(z) = element.parse::<u16>() {
            if self
                .db
                .data
                .elements
                .iter()
                .any(|e| e.atomic_number == z)
            {
                return Ok(z);
            }
        }
        // Try as symbol
        if let Some(&z) = self.db.symbol_to_z.get(element) {
            return Ok(z);
        }
        // Try as name (case-insensitive)
        if let Some(&z) = self.db.name_to_z.get(&element.to_lowercase()) {
            return Ok(z);
        }
        Err(XrayDbError::UnknownElement(element.to_string()))
    }

    fn element_record(&self, element: &str) -> Result<&xraydb_data::ElementRecord> {
        let z = self.resolve_element(element)?;
        self.db
            .data
            .elements
            .iter()
            .find(|e| e.atomic_number == z)
            .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))
    }

    pub fn atomic_number(&self, element: &str) -> Result<u16> {
        self.resolve_element(element)
    }

    pub fn symbol(&self, element: &str) -> Result<&str> {
        Ok(&self.element_record(element)?.symbol)
    }

    pub fn atomic_name(&self, element: &str) -> Result<&str> {
        Ok(&self.element_record(element)?.name)
    }

    pub fn molar_mass(&self, element: &str) -> Result<f64> {
        Ok(self.element_record(element)?.molar_mass)
    }

    pub fn density(&self, element: &str) -> Result<f64> {
        Ok(self.element_record(element)?.density)
    }
}

impl Default for XrayDb {
    fn default() -> Self {
        Self::new()
    }
}
