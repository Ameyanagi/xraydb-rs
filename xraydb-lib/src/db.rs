//! Database loading, indexing, and element resolution.

use std::collections::HashMap;
use std::sync::OnceLock;

use xraydb_data::XrayDatabase;

use crate::error::{Result, XrayDbError};
use crate::interp::safe_ln;

const COMPRESSED_DATA: &[u8] = include_bytes!("../data/xraydb.bin.zst");

/// Chantler tables pre-transformed into log space.
///
/// Building these once at load time removes the dominant cost from every Chantler
/// query: the naive implementation took the natural log of ~1425 abscissae *and*
/// ~1425 ordinates on every single call, which cost ~20 µs regardless of how many
/// energies were requested.
pub(crate) struct ChantlerLogs {
    pub(crate) log_energy: Vec<f64>,
    pub(crate) log_f2: Vec<f64>,
    pub(crate) log_photo: Vec<f64>,
    pub(crate) log_incoh: Vec<f64>,
    pub(crate) log_total: Vec<f64>,
}

impl ChantlerLogs {
    fn build(row: &xraydb_data::ChantlerRecord) -> Self {
        let ln_all = |v: &Vec<f64>| v.iter().copied().map(safe_ln).collect();
        ChantlerLogs {
            log_energy: row.energy.iter().copied().map(safe_ln).collect(),
            log_f2: ln_all(&row.f2),
            log_photo: ln_all(&row.mu_photo),
            log_incoh: ln_all(&row.mu_incoh),
            log_total: ln_all(&row.mu_total),
        }
    }
}

pub(crate) struct InitializedDb {
    pub(crate) data: XrayDatabase,

    // Element resolution
    symbol_to_z: HashMap<String, u16>,
    name_to_z: HashMap<String, u16>,
    z_to_element_idx: HashMap<u16, usize>,

    // Per-element table rows
    symbol_to_chantler_idx: HashMap<String, usize>,
    symbol_to_photo_idx: HashMap<String, usize>,
    symbol_to_scatter_idx: HashMap<String, usize>,
    ion_to_waasmaier_idx: HashMap<String, usize>,
    symbol_to_waasmaier_idxs: HashMap<String, Vec<usize>>,

    // Indices added to replace full-table linear scans
    levels_by_symbol: HashMap<String, Vec<usize>>,
    level_by_symbol_edge: HashMap<(String, String), usize>,
    transitions_by_symbol: HashMap<String, Vec<usize>>,
    widths_by_z: HashMap<u16, Vec<usize>>,
    ck_by_key: HashMap<(String, String, String), usize>,
    gas_to_potential: HashMap<String, f64>,
    /// Indices into `data.xray_levels`, sorted by absorption edge energy.
    levels_by_energy: Vec<usize>,

    // Pre-computed log-space Chantler tables, parallel to `data.chantler`.
    chantler_logs: Vec<ChantlerLogs>,
}

static DATABASE: OnceLock<core::result::Result<InitializedDb, String>> = OnceLock::new();

/// Decode the embedded blob. Separated out so the error path stays `?`-shaped.
fn decode_embedded() -> core::result::Result<XrayDatabase, String> {
    let mut decoder = ruzstd::decoding::StreamingDecoder::new(COMPRESSED_DATA)
        .map_err(|e| format!("failed to create zstd decoder: {e}"))?;
    let mut decompressed = Vec::new();
    std::io::Read::read_to_end(&mut decoder, &mut decompressed)
        .map_err(|e| format!("failed to decompress embedded database: {e}"))?;
    postcard::from_bytes(&decompressed)
        .map_err(|e| format!("failed to deserialize embedded database: {e}"))
}

fn build() -> core::result::Result<InitializedDb, String> {
    let data = decode_embedded()?;

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
    let chantler_logs: Vec<ChantlerLogs> = data.chantler.iter().map(ChantlerLogs::build).collect();

    let mut symbol_to_photo_idx = HashMap::with_capacity(data.photoabsorption.len());
    for (idx, row) in data.photoabsorption.iter().enumerate() {
        symbol_to_photo_idx.insert(row.element.clone(), idx);
    }

    let mut symbol_to_scatter_idx = HashMap::with_capacity(data.scattering.len());
    for (idx, row) in data.scattering.iter().enumerate() {
        symbol_to_scatter_idx.insert(row.element.clone(), idx);
    }

    let mut ion_to_waasmaier_idx = HashMap::with_capacity(data.waasmaier.len());
    let mut symbol_to_waasmaier_idxs: HashMap<String, Vec<usize>> =
        HashMap::with_capacity(data.elements.len());
    for (idx, row) in data.waasmaier.iter().enumerate() {
        ion_to_waasmaier_idx.insert(row.ion.clone(), idx);
        symbol_to_waasmaier_idxs
            .entry(row.element.clone())
            .or_default()
            .push(idx);
    }

    let mut levels_by_symbol: HashMap<String, Vec<usize>> = HashMap::new();
    let mut level_by_symbol_edge = HashMap::with_capacity(data.xray_levels.len());
    for (idx, level) in data.xray_levels.iter().enumerate() {
        levels_by_symbol
            .entry(level.element.clone())
            .or_default()
            .push(idx);
        level_by_symbol_edge.insert((level.element.clone(), level.iupac_symbol.clone()), idx);
    }

    let mut levels_by_energy: Vec<usize> = (0..data.xray_levels.len()).collect();
    levels_by_energy.sort_by(|&a, &b| {
        data.xray_levels[a]
            .absorption_edge
            .total_cmp(&data.xray_levels[b].absorption_edge)
    });

    let mut transitions_by_symbol: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, trans) in data.xray_transitions.iter().enumerate() {
        transitions_by_symbol
            .entry(trans.element.clone())
            .or_default()
            .push(idx);
    }

    let mut widths_by_z: HashMap<u16, Vec<usize>> = HashMap::new();
    for (idx, w) in data.corelevel_widths.iter().enumerate() {
        widths_by_z.entry(w.atomic_number).or_default().push(idx);
    }

    let mut ck_by_key = HashMap::with_capacity(data.coster_kronig.len());
    for (idx, ck) in data.coster_kronig.iter().enumerate() {
        ck_by_key.insert(
            (
                ck.element.clone(),
                ck.initial_level.clone(),
                ck.final_level.clone(),
            ),
            idx,
        );
    }

    let mut gas_to_potential = HashMap::with_capacity(data.ionization_potentials.len());
    for ip in &data.ionization_potentials {
        gas_to_potential.insert(ip.gas.to_lowercase(), ip.potential);
    }

    Ok(InitializedDb {
        data,
        symbol_to_z,
        name_to_z,
        z_to_element_idx,
        symbol_to_chantler_idx,
        symbol_to_photo_idx,
        symbol_to_scatter_idx,
        ion_to_waasmaier_idx,
        symbol_to_waasmaier_idxs,
        levels_by_symbol,
        level_by_symbol_edge,
        transitions_by_symbol,
        widths_by_z,
        ck_by_key,
        gas_to_potential,
        levels_by_energy,
        chantler_logs,
    })
}

fn db() -> core::result::Result<&'static InitializedDb, &'static str> {
    match DATABASE.get_or_init(build) {
        Ok(inner) => Ok(inner),
        Err(msg) => Err(msg.as_str()),
    }
}

/// The main interface to the X-ray database.
///
/// `XrayDb` is a zero-sized handle over a lazily-initialised `static`: constructing one
/// is free after the first call, and it is `Copy`, so pass it by value.
///
/// ```
/// use xraydb::XrayDb;
/// let db = XrayDb::try_new()?;
/// assert_eq!(db.atomic_number("Fe")?, 26);
/// # Ok::<(), xraydb::XrayDbError>(())
/// ```
#[derive(Clone, Copy)]
pub struct XrayDb {
    db: &'static InitializedDb,
}

impl std::fmt::Debug for XrayDb {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("XrayDb")
            .field("elements", &self.db.data.elements.len())
            .field("version", &self.data_version())
            .finish()
    }
}

impl XrayDb {
    /// Open the database, decoding the embedded blob on first use.
    ///
    /// Returns [`XrayDbError::DataError`] if the embedded data cannot be decoded, which
    /// indicates a corrupted build rather than a caller mistake.
    pub fn try_new() -> Result<Self> {
        match db() {
            Ok(inner) => Ok(XrayDb { db: inner }),
            Err(msg) => Err(XrayDbError::DataError(msg.to_string())),
        }
    }

    /// Open the database, panicking if the embedded blob is corrupt.
    ///
    /// # Panics
    ///
    /// Panics if the compiled-in database fails to decompress or deserialize. This can
    /// only happen if the build is corrupt; use [`XrayDb::try_new`] to handle it.
    pub fn new() -> Self {
        match Self::try_new() {
            Ok(db) => db,
            Err(e) => panic!("xraydb: embedded database is corrupt: {e}"),
        }
    }

    /// Access the raw deserialized database.
    pub fn raw(&self) -> &'static XrayDatabase {
        &self.db.data
    }

    /// Resolve an element identifier — symbol, name, or atomic number — to its Z.
    ///
    /// Matching is case-insensitive: `"Fe"`, `"fe"`, `"iron"`, `"Iron"`, and `"26"` all
    /// resolve to 26.
    pub fn resolve_element(&self, element: &str) -> Result<u16> {
        // Atomic number, e.g. "26".
        if let Ok(z) = element.parse::<u16>()
            && self.db.z_to_element_idx.contains_key(&z)
        {
            return Ok(z);
        }

        // Exact symbol/name first, so the common path allocates nothing.
        if let Some(&z) = self.db.symbol_to_z.get(element) {
            return Ok(z);
        }
        if let Some(&z) = self.db.name_to_z.get(element) {
            return Ok(z);
        }

        let lower = element.to_lowercase();
        if let Some(&z) = self.db.symbol_to_z.get(&lower) {
            return Ok(z);
        }
        if let Some(&z) = self.db.name_to_z.get(&lower) {
            return Ok(z);
        }

        Err(XrayDbError::UnknownElement(element.to_string()))
    }

    fn element_record(&self, element: &str) -> Result<&'static xraydb_data::ElementRecord> {
        let z = self.resolve_element(element)?;
        self.element_by_z(z)
            .ok_or_else(|| XrayDbError::UnknownElement(element.to_string()))
    }

    /// Atomic number for an element identifier.
    pub fn atomic_number(&self, element: &str) -> Result<u16> {
        self.resolve_element(element)
    }

    /// Canonical element symbol (e.g. `"Fe"`) for any accepted identifier.
    pub fn symbol(&self, element: &str) -> Result<&'static str> {
        Ok(&self.element_record(element)?.symbol)
    }

    /// Full element name (e.g. `"iron"`) for any accepted identifier.
    pub fn atomic_name(&self, element: &str) -> Result<&'static str> {
        Ok(&self.element_record(element)?.name)
    }

    /// Molar mass in g/mol.
    pub fn molar_mass(&self, element: &str) -> Result<f64> {
        Ok(self.element_record(element)?.molar_mass)
    }

    /// Elemental density in g/cm³ at standard conditions.
    pub fn density(&self, element: &str) -> Result<f64> {
        Ok(self.element_record(element)?.density)
    }

    /// Number of elements in the database.
    pub fn element_count(&self) -> usize {
        self.db.data.elements.len()
    }

    /// Upstream XrayDB version tag the embedded data was generated from.
    pub fn data_version(&self) -> Option<&'static str> {
        self.db.data.version.first().map(|v| v.tag.as_str())
    }

    // ── Internal accessors ────────────────────────────────────────────────

    pub(crate) fn element_by_z(&self, z: u16) -> Option<&'static xraydb_data::ElementRecord> {
        self.db
            .z_to_element_idx
            .get(&z)
            .and_then(|&idx| self.db.data.elements.get(idx))
    }

    pub(crate) fn chantler_by_symbol(
        &self,
        symbol: &str,
    ) -> Option<(&'static xraydb_data::ChantlerRecord, &'static ChantlerLogs)> {
        let &idx = self.db.symbol_to_chantler_idx.get(symbol)?;
        let row = self.db.data.chantler.get(idx)?;
        let logs = self.db.chantler_logs.get(idx)?;
        Some((row, logs))
    }

    pub(crate) fn photo_by_symbol(
        &self,
        symbol: &str,
    ) -> Option<&'static xraydb_data::PhotoabsorptionRecord> {
        self.db
            .symbol_to_photo_idx
            .get(symbol)
            .and_then(|&idx| self.db.data.photoabsorption.get(idx))
    }

    pub(crate) fn scatter_by_symbol(
        &self,
        symbol: &str,
    ) -> Option<&'static xraydb_data::ScatteringRecord> {
        self.db
            .symbol_to_scatter_idx
            .get(symbol)
            .and_then(|&idx| self.db.data.scattering.get(idx))
    }

    pub(crate) fn waasmaier_by_ion(
        &self,
        ion: &str,
    ) -> Option<&'static xraydb_data::WaasmaierRecord> {
        self.db
            .ion_to_waasmaier_idx
            .get(ion)
            .and_then(|&idx| self.db.data.waasmaier.get(idx))
    }

    pub(crate) fn waasmaier_indices_by_symbol(&self, symbol: &str) -> &'static [usize] {
        self.db
            .symbol_to_waasmaier_idxs
            .get(symbol)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub(crate) fn level_indices(&self, symbol: &str) -> &'static [usize] {
        self.db
            .levels_by_symbol
            .get(symbol)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub(crate) fn level_by_edge(
        &self,
        symbol: &str,
        edge: &str,
    ) -> Option<&'static xraydb_data::XrayLevelRecord> {
        // The tuple key would otherwise force two allocations per lookup; the borrowed
        // form is only usable with an owned key, so allocate once and reuse.
        let key = (symbol.to_string(), edge.to_string());
        let &idx = self.db.level_by_symbol_edge.get(&key)?;
        self.db.data.xray_levels.get(idx)
    }

    pub(crate) fn transition_indices(&self, symbol: &str) -> &'static [usize] {
        self.db
            .transitions_by_symbol
            .get(symbol)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub(crate) fn width_indices(&self, z: u16) -> &'static [usize] {
        self.db
            .widths_by_z
            .get(&z)
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    pub(crate) fn coster_kronig_record(
        &self,
        symbol: &str,
        initial: &str,
        final_level: &str,
    ) -> Option<&'static xraydb_data::CosterKronigRecord> {
        let key = (
            symbol.to_string(),
            initial.to_string(),
            final_level.to_string(),
        );
        let &idx = self.db.ck_by_key.get(&key)?;
        self.db.data.coster_kronig.get(idx)
    }

    pub(crate) fn gas_potential(&self, gas: &str) -> Option<f64> {
        self.db.gas_to_potential.get(&gas.to_lowercase()).copied()
    }

    pub(crate) fn known_gases(&self) -> Vec<String> {
        let mut gases: Vec<String> = self
            .db
            .data
            .ionization_potentials
            .iter()
            .map(|ip| ip.gas.clone())
            .collect();
        gases.sort();
        gases
    }

    pub(crate) fn levels_sorted_by_energy(&self) -> &'static [usize] {
        &self.db.levels_by_energy
    }

    pub(crate) fn contains_symbol(&self, symbol: &str) -> bool {
        self.db.symbol_to_z.contains_key(symbol)
    }
}

impl Default for XrayDb {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_new_succeeds_on_the_shipped_blob() {
        assert!(XrayDb::try_new().is_ok());
    }

    #[test]
    fn indices_agree_with_a_linear_scan() {
        let db = XrayDb::new();
        let raw = db.raw();

        // xray_levels index
        for sym in ["Fe", "U", "H", "Cu"] {
            let scanned: Vec<usize> = raw
                .xray_levels
                .iter()
                .enumerate()
                .filter(|(_, l)| l.element == sym)
                .map(|(i, _)| i)
                .collect();
            assert_eq!(db.level_indices(sym), scanned.as_slice(), "levels {sym}");

            let scanned: Vec<usize> = raw
                .xray_transitions
                .iter()
                .enumerate()
                .filter(|(_, t)| t.element == sym)
                .map(|(i, _)| i)
                .collect();
            assert_eq!(
                db.transition_indices(sym),
                scanned.as_slice(),
                "transitions {sym}"
            );
        }

        // corelevel_widths index
        for z in [1u16, 26, 92] {
            let scanned: Vec<usize> = raw
                .corelevel_widths
                .iter()
                .enumerate()
                .filter(|(_, w)| w.atomic_number == z)
                .map(|(i, _)| i)
                .collect();
            assert_eq!(db.width_indices(z), scanned.as_slice(), "widths Z={z}");
        }
    }

    #[test]
    fn levels_by_energy_is_sorted() {
        let db = XrayDb::new();
        let raw = db.raw();
        let sorted = db.levels_sorted_by_energy();
        assert_eq!(sorted.len(), raw.xray_levels.len());
        for w in sorted.windows(2) {
            assert!(raw.xray_levels[w[0]].absorption_edge <= raw.xray_levels[w[1]].absorption_edge);
        }
    }

    #[test]
    fn chantler_logs_are_parallel_to_chantler_rows() {
        let db = XrayDb::new();
        assert_eq!(db.db.chantler_logs.len(), db.raw().chantler.len());
        for (row, logs) in db.raw().chantler.iter().zip(db.db.chantler_logs.iter()) {
            assert_eq!(row.energy.len(), logs.log_energy.len(), "{}", row.element);
            assert_eq!(row.f2.len(), logs.log_f2.len(), "{}", row.element);
            assert_eq!(row.mu_total.len(), logs.log_total.len(), "{}", row.element);
        }
    }

    #[test]
    fn resolution_is_case_insensitive() {
        let db = XrayDb::new();
        for id in ["Fe", "fe", "FE", "iron", "Iron", "IRON", "26"] {
            assert_eq!(db.atomic_number(id).ok(), Some(26), "{id}");
        }
    }
}
