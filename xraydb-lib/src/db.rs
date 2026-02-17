use std::collections::HashMap;
use std::sync::OnceLock;

use xraydb_data::XrayDatabase;

use crate::error::{Result, XrayDbError};

const COMPRESSED_DATA: &[u8] = include_bytes!("../data/xraydb.bin.zst");

struct InitializedDb {
    data: XrayDatabase,
    symbol_to_z: HashMap<String, u16>,
    name_to_z: HashMap<String, u16>,
    z_to_element_idx: HashMap<u16, usize>,
    symbol_to_chantler_idx: HashMap<String, usize>,
    symbol_to_photo_idx: HashMap<String, usize>,
    symbol_to_scatter_idx: HashMap<String, usize>,
    ion_to_waasmaier_idx: HashMap<String, usize>,
    symbol_to_waasmaier_idxs: HashMap<String, Vec<usize>>,
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
        let mut symbol_to_z = HashMap::with_capacity(data.elements.len() * 2);
        let mut name_to_z = HashMap::with_capacity(data.elements.len() * 2);
        let mut z_to_element_idx = HashMap::with_capacity(data.elements.len());
        for (idx, elem) in data.elements.iter().enumerate() {
            symbol_to_z.insert(elem.symbol.clone(), elem.atomic_number);
            symbol_to_z.insert(elem.symbol.to_lowercase(), elem.atomic_number);
            name_to_z.insert(elem.name.clone(), elem.atomic_number);
            name_to_z.insert(elem.name.to_lowercase(), elem.atomic_number);
            z_to_element_idx.insert(elem.atomic_number, idx);
        }

        let mut symbol_to_chantler_idx = HashMap::with_capacity(data.chantler.len());
        for (idx, row) in data.chantler.iter().enumerate() {
            symbol_to_chantler_idx.insert(row.element.clone(), idx);
        }

        let mut symbol_to_photo_idx = HashMap::with_capacity(data.photoabsorption.len());
        for (idx, row) in data.photoabsorption.iter().enumerate() {
            symbol_to_photo_idx.insert(row.element.clone(), idx);
        }

        let mut symbol_to_scatter_idx = HashMap::with_capacity(data.scattering.len());
        for (idx, row) in data.scattering.iter().enumerate() {
            symbol_to_scatter_idx.insert(row.element.clone(), idx);
        }

        let mut ion_to_waasmaier_idx = HashMap::with_capacity(data.waasmaier.len());
        let mut symbol_to_waasmaier_idxs = HashMap::with_capacity(data.elements.len());
        for (idx, row) in data.waasmaier.iter().enumerate() {
            ion_to_waasmaier_idx.insert(row.ion.clone(), idx);
            symbol_to_waasmaier_idxs
                .entry(row.element.clone())
                .or_insert_with(Vec::new)
                .push(idx);
        }

        InitializedDb {
            data,
            symbol_to_z,
            name_to_z,
            z_to_element_idx,
            symbol_to_chantler_idx,
            symbol_to_photo_idx,
            symbol_to_scatter_idx,
            ion_to_waasmaier_idx,
            symbol_to_waasmaier_idxs,
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
        if let Ok(z) = element.parse::<u16>()
            && self.db.z_to_element_idx.contains_key(&z)
        {
            return Ok(z);
        }

        // Try as symbol
        if let Some(&z) = self.db.symbol_to_z.get(element) {
            return Ok(z);
        }

        let lower = element.to_lowercase();
        if let Some(&z) = self.db.symbol_to_z.get(&lower) {
            return Ok(z);
        }

        // Try as name (case-insensitive)
        if let Some(&z) = self.db.name_to_z.get(element) {
            return Ok(z);
        }
        if let Some(&z) = self.db.name_to_z.get(&lower) {
            return Ok(z);
        }

        Err(XrayDbError::UnknownElement(element.to_string()))
    }

    fn element_record(&self, element: &str) -> Result<&xraydb_data::ElementRecord> {
        let z = self.resolve_element(element)?;
        self.element_by_z(z)
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

    pub(crate) fn element_by_z(&self, z: u16) -> Option<&xraydb_data::ElementRecord> {
        self.db
            .z_to_element_idx
            .get(&z)
            .and_then(|&idx| self.db.data.elements.get(idx))
    }

    pub(crate) fn chantler_by_symbol(&self, symbol: &str) -> Option<&xraydb_data::ChantlerRecord> {
        self.db
            .symbol_to_chantler_idx
            .get(symbol)
            .and_then(|&idx| self.db.data.chantler.get(idx))
    }

    pub(crate) fn photo_by_symbol(
        &self,
        symbol: &str,
    ) -> Option<&xraydb_data::PhotoabsorptionRecord> {
        self.db
            .symbol_to_photo_idx
            .get(symbol)
            .and_then(|&idx| self.db.data.photoabsorption.get(idx))
    }

    pub(crate) fn scatter_by_symbol(&self, symbol: &str) -> Option<&xraydb_data::ScatteringRecord> {
        self.db
            .symbol_to_scatter_idx
            .get(symbol)
            .and_then(|&idx| self.db.data.scattering.get(idx))
    }

    pub(crate) fn waasmaier_by_ion(&self, ion: &str) -> Option<&xraydb_data::WaasmaierRecord> {
        self.db
            .ion_to_waasmaier_idx
            .get(ion)
            .and_then(|&idx| self.db.data.waasmaier.get(idx))
    }

    pub(crate) fn waasmaier_indices_by_symbol(&self, symbol: &str) -> Option<&[usize]> {
        self.db
            .symbol_to_waasmaier_idxs
            .get(symbol)
            .map(Vec::as_slice)
    }
}

impl Default for XrayDb {
    fn default() -> Self {
        Self::new()
    }
}
