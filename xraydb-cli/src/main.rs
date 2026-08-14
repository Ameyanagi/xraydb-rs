//! `xraydb` — command-line access to the X-ray reference database.

#![forbid(unsafe_code)]
#![deny(clippy::unwrap_used, clippy::expect_used)]
#![cfg_attr(test, allow(clippy::unwrap_used, clippy::expect_used))]

mod energy;
mod output;

use std::io::{self, Write};

use anyhow::{Context, Result, bail};
use clap::{Parser, Subcommand, ValueEnum};
use serde::Serialize;

use output::{Format, Table, num, write_record};
use xraydb::{ChantlerKind, CrossSectionKind, XrayDb};

/// X-ray reference data for the elements.
///
/// Element identifiers accept a symbol (`Fe`), a name (`iron`), or an atomic
/// number (`26`), case-insensitively.
#[derive(Parser)]
#[command(name = "xraydb", version, about, long_about = None)]
struct Cli {
    /// Output format.
    #[arg(long, short = 'f', value_enum, default_value_t = Format::Text, global = true)]
    format: Format,

    /// Shorthand for `--format json`.
    #[arg(long, global = true, conflicts_with_all = ["format", "csv"])]
    json: bool,

    /// Shorthand for `--format csv`.
    #[arg(long, global = true, conflicts_with_all = ["format", "json"])]
    csv: bool,

    #[command(subcommand)]
    command: Command,
}

impl Cli {
    fn resolved_format(&self) -> Format {
        if self.json {
            Format::Json
        } else if self.csv {
            Format::Csv
        } else {
            self.format
        }
    }
}

/// Which Elam cross-section to report.
#[derive(Debug, Clone, Copy, ValueEnum, Default)]
enum Kind {
    /// Photoelectric + coherent + incoherent.
    #[default]
    Total,
    /// Photoelectric absorption only.
    Photo,
    /// Coherent (Rayleigh) scattering only.
    Coherent,
    /// Incoherent (Compton) scattering only.
    Incoherent,
}

impl From<Kind> for CrossSectionKind {
    fn from(k: Kind) -> Self {
        match k {
            Kind::Total => CrossSectionKind::Total,
            Kind::Photo => CrossSectionKind::Photo,
            Kind::Coherent => CrossSectionKind::Coherent,
            Kind::Incoherent => CrossSectionKind::Incoherent,
        }
    }
}

/// Which Chantler attenuation coefficient to report.
#[derive(Debug, Clone, Copy, ValueEnum, Default)]
enum CKind {
    /// Total attenuation.
    #[default]
    Total,
    /// Photoelectric only.
    Photo,
    /// Incoherent only.
    Incoherent,
}

impl From<CKind> for ChantlerKind {
    fn from(k: CKind) -> Self {
        match k {
            CKind::Total => ChantlerKind::Total,
            CKind::Photo => ChantlerKind::Photo,
            CKind::Incoherent => ChantlerKind::Incoherent,
        }
    }
}

const ENERGY_HELP: &str = "Energies in eV. Accepts a single value (10000), a list \
(5000,7112,10000), a range (5000:15000:100 as start:stop:step), or a log-spaced \
range (5000:15000/50 as start:stop/count).";

#[derive(Subcommand)]
enum Command {
    /// Element facts: Z, symbol, name, molar mass, density.
    #[command(after_help = "EXAMPLE:\n  xraydb element Fe\n  xraydb element 26 --json")]
    Element {
        /// Element symbol, name, or atomic number.
        element: String,
    },

    /// All absorption edges for an element.
    #[command(after_help = "EXAMPLE:\n  xraydb edges Fe\n  xraydb edges uranium --csv")]
    Edges {
        /// Element symbol, name, or atomic number.
        element: String,
    },

    /// One absorption edge.
    #[command(after_help = "EXAMPLE:\n  xraydb edge Fe K")]
    Edge {
        /// Element symbol, name, or atomic number.
        element: String,
        /// IUPAC edge label (K, L1, L2, L3, M1, ...).
        edge: String,
    },

    /// Emission lines, strongest first.
    #[command(
        after_help = "EXAMPLE:\n  xraydb lines Cu --level K\n  xraydb lines Au --excitation 12000"
    )]
    Lines {
        /// Element symbol, name, or atomic number.
        element: String,
        /// Restrict to lines from this initial level.
        #[arg(long)]
        level: Option<String>,
        /// Drop lines that cannot be excited at this energy (eV).
        #[arg(long)]
        excitation: Option<f64>,
        /// Show at most this many lines.
        #[arg(long, default_value_t = 20)]
        limit: usize,
    },

    /// Mass or linear attenuation coefficient.
    ///
    /// With `--density`, reports linear attenuation µ in 1/cm for a compound.
    /// Without it, reports mass attenuation µ/ρ in cm²/g for a single element.
    #[command(after_help = "EXAMPLE:\n  xraydb mu Fe --energy 10000\n  \
xraydb mu H2O --density 1.0 --energy 5000:15000/11\n  \
xraydb mu kapton --energy 10000    # density from the materials table")]
    Mu {
        /// Element identifier, chemical formula, or material name.
        target: String,
        #[arg(long, short = 'e', help = ENERGY_HELP)]
        energy: String,
        /// Density in g/cm³. Required for formulas not in the materials table.
        #[arg(long, short = 'd')]
        density: Option<f64>,
        /// Which cross-section to report.
        #[arg(long, short = 'k', value_enum, default_value_t = Kind::Total)]
        kind: Kind,
    },

    /// Chantler anomalous scattering factors f1 and f2.
    #[command(after_help = "EXAMPLE:\n  xraydb f1f2 Fe --energy 7000:7200:20")]
    F1f2 {
        /// Element symbol, name, or atomic number.
        element: String,
        #[arg(long, short = 'e', help = ENERGY_HELP)]
        energy: String,
    },

    /// Chantler mass attenuation coefficient.
    #[command(after_help = "EXAMPLE:\n  xraydb mu-chantler Fe --energy 10000")]
    MuChantler {
        /// Element symbol, name, or atomic number.
        element: String,
        #[arg(long, short = 'e', help = ENERGY_HELP)]
        energy: String,
        /// Which coefficient to report.
        #[arg(long, short = 'k', value_enum, default_value_t = CKind::Total)]
        kind: CKind,
    },

    /// Waasmaier-Kirfel elastic scattering factor f0.
    #[command(after_help = "EXAMPLE:\n  xraydb f0 Fe2+ --q 0:6/13\n  xraydb f0 --list Fe")]
    F0 {
        /// Ion name, e.g. `Fe`, `Fe2+`, `O1-`.
        ion: Option<String>,
        /// Momentum transfer sin(theta)/lambda in 1/Angstrom. Same syntax as --energy.
        #[arg(long, short = 'q', default_value = "0:6/13")]
        q: String,
        /// List available ions instead (optionally for one element).
        #[arg(long)]
        list: bool,
    },

    /// Refractive index decrements delta and beta.
    #[command(after_help = "EXAMPLE:\n  xraydb delta-beta SiO2 --density 2.2 --energy 10000")]
    DeltaBeta {
        /// Chemical formula or material name.
        formula: String,
        #[arg(long, short = 'e', help = ENERGY_HELP)]
        energy: String,
        /// Density in g/cm³. Taken from the materials table when omitted.
        #[arg(long, short = 'd')]
        density: Option<f64>,
    },

    /// Look up one material by name or formula.
    #[command(after_help = "EXAMPLE:\n  xraydb material kapton")]
    Material {
        /// Material name or chemical formula.
        name: String,
    },

    /// List the built-in materials table.
    #[command(after_help = "EXAMPLE:\n  xraydb materials --filter poly")]
    Materials {
        /// Only show materials whose name or formula contains this text.
        #[arg(long)]
        filter: Option<String>,
    },

    /// Identify the absorption edge nearest an energy.
    #[command(after_help = "EXAMPLE:\n  xraydb guess-edge --energy 7100\n  \
xraydb guess-edge --energy 7100 --edges K,L3 --tolerance 50")]
    GuessEdge {
        /// Energy in eV.
        #[arg(long, short = 'e')]
        energy: f64,
        /// Comma-separated edge labels to consider.
        #[arg(long, default_value = "K,L3,L2,L1,M5")]
        edges: String,
        /// Reject matches further than this many eV away.
        #[arg(long, short = 't')]
        tolerance: Option<f64>,
    },

    /// Core-hole level widths.
    #[command(after_help = "EXAMPLE:\n  xraydb core-width Au\n  xraydb core-width Au --edge L3")]
    CoreWidth {
        /// Element symbol, name, or atomic number.
        element: String,
        /// Restrict to one edge.
        #[arg(long)]
        edge: Option<String>,
    },

    /// Ion-chamber photon fluxes from a measured current.
    #[command(
        after_help = "EXAMPLE:\n  xraydb ionchamber --gas nitrogen --volts 1.0 \
--length 30 --energy 10000 --sensitivity 1e-6"
    )]
    Ionchamber {
        /// Gas mixture as name:fraction, repeatable (e.g. --gas nitrogen:0.8 --gas argon:0.2).
        #[arg(long = "gas", required = true)]
        gases: Vec<String>,
        /// Measured voltage.
        #[arg(long)]
        volts: f64,
        /// Active chamber length in cm.
        #[arg(long, short = 'l')]
        length: f64,
        /// X-ray energy in eV.
        #[arg(long, short = 'e')]
        energy: f64,
        /// Amplifier sensitivity in A/V.
        #[arg(long, short = 's')]
        sensitivity: f64,
        /// Exclude the Compton-electron contribution.
        #[arg(long)]
        no_compton: bool,
        /// Count only one carrier species instead of two.
        #[arg(long)]
        single_carrier: bool,
    },

    /// Compton scattering energies.
    #[command(after_help = "EXAMPLE:\n  xraydb compton --energy 20000")]
    Compton {
        /// Incident energy in eV.
        #[arg(long, short = 'e')]
        energy: f64,
    },

    /// Database provenance and size.
    Info,
}

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    match run(&cli) {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(err) => {
            // Broken pipe (e.g. `| head`) is a normal way for a CLI to end.
            if let Some(io_err) = err.downcast_ref::<io::Error>()
                && io_err.kind() == io::ErrorKind::BrokenPipe
            {
                return std::process::ExitCode::SUCCESS;
            }
            eprintln!("error: {err:#}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run(cli: &Cli) -> Result<()> {
    let db = XrayDb::try_new().context("opening the embedded database")?;
    let format = cli.resolved_format();
    let stdout = io::stdout();
    let mut out = stdout.lock();

    match &cli.command {
        Command::Element { element } => cmd_element(&db, &mut out, format, element),
        Command::Edges { element } => cmd_edges(&db, &mut out, format, element),
        Command::Edge { element, edge } => cmd_edge(&db, &mut out, format, element, edge),
        Command::Lines {
            element,
            level,
            excitation,
            limit,
        } => cmd_lines(
            &db,
            &mut out,
            format,
            element,
            level.as_deref(),
            *excitation,
            *limit,
        ),
        Command::Mu {
            target,
            energy,
            density,
            kind,
        } => cmd_mu(&db, &mut out, format, target, energy, *density, *kind),
        Command::F1f2 { element, energy } => cmd_f1f2(&db, &mut out, format, element, energy),
        Command::MuChantler {
            element,
            energy,
            kind,
        } => cmd_mu_chantler(&db, &mut out, format, element, energy, *kind),
        Command::F0 { ion, q, list } => cmd_f0(&db, &mut out, format, ion.as_deref(), q, *list),
        Command::DeltaBeta {
            formula,
            energy,
            density,
        } => cmd_delta_beta(&db, &mut out, format, formula, energy, *density),
        Command::Material { name } => cmd_material(&db, &mut out, format, name),
        Command::Materials { filter } => cmd_materials(&db, &mut out, format, filter.as_deref()),
        Command::GuessEdge {
            energy,
            edges,
            tolerance,
        } => cmd_guess_edge(&db, &mut out, format, *energy, edges, *tolerance),
        Command::CoreWidth { element, edge } => {
            cmd_core_width(&db, &mut out, format, element, edge.as_deref())
        }
        Command::Ionchamber {
            gases,
            volts,
            length,
            energy,
            sensitivity,
            no_compton,
            single_carrier,
        } => cmd_ionchamber(
            &db,
            &mut out,
            format,
            gases,
            *volts,
            *length,
            *energy,
            *sensitivity,
            !*no_compton,
            !*single_carrier,
        ),
        Command::Compton { energy } => cmd_compton(&db, &mut out, format, *energy),
        Command::Info => cmd_info(&db, &mut out, format),
    }
}

// ── Commands ────────────────────────────────────────────────────────────────

#[derive(Serialize)]
struct ElementRecord<'a> {
    z: u16,
    symbol: &'a str,
    name: &'a str,
    molar_mass: f64,
    density: f64,
}

fn cmd_element(db: &XrayDb, out: &mut impl Write, format: Format, element: &str) -> Result<()> {
    let z = db.atomic_number(element)?;
    let symbol = db.symbol(element)?;
    let name = db.atomic_name(element)?;
    let molar_mass = db.molar_mass(element)?;
    let density = db.density(element)?;

    let record = ElementRecord {
        z,
        symbol,
        name,
        molar_mass,
        density,
    };
    let pairs = [
        ("atomic number", z.to_string()),
        ("symbol", symbol.to_string()),
        ("name", name.to_string()),
        ("molar mass (g/mol)", num(molar_mass)),
        ("density (g/cm^3)", num(density)),
    ];
    write_record(out, format, &pairs, &record)
}

fn cmd_edges(db: &XrayDb, out: &mut impl Write, format: Format, element: &str) -> Result<()> {
    let edges = db.xray_edges(element)?;
    let mut sorted: Vec<_> = edges.into_iter().collect();
    sorted.sort_by(|a, b| b.1.energy.total_cmp(&a.1.energy));

    let mut table = Table::new(&["edge", "energy_eV", "fluor_yield", "jump_ratio"]);
    for (label, edge) in sorted {
        table.row(vec![
            label,
            num(edge.energy),
            num(edge.fluorescence_yield),
            num(edge.jump_ratio),
        ]);
    }
    table.write(out, format)
}

fn cmd_edge(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    element: &str,
    edge: &str,
) -> Result<()> {
    let e = db.xray_edge(element, edge)?;

    #[derive(Serialize)]
    struct Rec {
        edge: String,
        energy_ev: f64,
        fluorescence_yield: f64,
        jump_ratio: f64,
    }
    let record = Rec {
        edge: edge.to_string(),
        energy_ev: e.energy,
        fluorescence_yield: e.fluorescence_yield,
        jump_ratio: e.jump_ratio,
    };
    let pairs = [
        ("edge", edge.to_string()),
        ("energy (eV)", num(e.energy)),
        ("fluorescence yield", num(e.fluorescence_yield)),
        ("jump ratio", num(e.jump_ratio)),
    ];
    write_record(out, format, &pairs, &record)
}

fn cmd_lines(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    element: &str,
    level: Option<&str>,
    excitation: Option<f64>,
    limit: usize,
) -> Result<()> {
    let lines = db.xray_lines_by_intensity(element, level, excitation)?;

    let mut table = Table::new(&["line", "energy_eV", "intensity", "initial", "final"]);
    for (label, line) in lines.into_iter().take(limit) {
        table.row(vec![
            label,
            num(line.energy),
            num(line.intensity),
            line.initial_level,
            line.final_level,
        ]);
    }
    table.write(out, format)
}

/// Resolve a target that may be an element, a material name, or a bare formula.
fn resolve_target(db: &XrayDb, target: &str, density: Option<f64>) -> Result<(String, f64)> {
    if let Some(material) = db.find_material(target) {
        return Ok((
            material.formula.to_string(),
            density.unwrap_or(material.density),
        ));
    }
    if let Ok(symbol) = db.symbol(target) {
        let d = match density {
            Some(d) => d,
            None => db.density(target)?,
        };
        return Ok((symbol.to_string(), d));
    }
    match density {
        Some(d) => Ok((target.to_string(), d)),
        None => bail!(
            "'{target}' is not a known element or material, so --density is required \
             to interpret it as a formula"
        ),
    }
}

fn cmd_mu(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    target: &str,
    energy_spec: &str,
    density: Option<f64>,
    kind: Kind,
) -> Result<()> {
    let energies = energy::parse_spec(energy_spec)?;
    let (formula, resolved_density) = resolve_target(db, target, density)?;

    // Mass attenuation (cm^2/g) is the density-independent quantity; linear
    // attenuation (1/cm) is what you get after multiplying by density.
    let mu_linear = db.material_mu(&formula, resolved_density, &energies, kind.into())?;

    let mut table = Table::new(&["energy_eV", "mu_cm2_per_g", "mu_per_cm", "atten_len_cm"]);
    for (e, mu) in energies.iter().zip(mu_linear.iter()) {
        let mass_mu = mu / resolved_density;
        let atten = if *mu > 0.0 { 1.0 / mu } else { f64::INFINITY };
        table.row(vec![num(*e), num(mass_mu), num(*mu), num(atten)]);
    }
    table.write(out, format)
}

fn cmd_f1f2(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    element: &str,
    energy_spec: &str,
) -> Result<()> {
    let energies = energy::parse_spec(energy_spec)?;
    let f1 = db.f1_chantler(element, &energies)?;
    let f2 = db.f2_chantler(element, &energies)?;
    let z = f64::from(db.atomic_number(element)?);

    let mut table = Table::new(&["energy_eV", "f1_prime", "f1_total", "f2"]);
    for i in 0..energies.len() {
        table.row(vec![
            num(energies[i]),
            num(f1[i]),
            num(z + f1[i]),
            num(f2[i]),
        ]);
    }
    table.write(out, format)
}

fn cmd_mu_chantler(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    element: &str,
    energy_spec: &str,
    kind: CKind,
) -> Result<()> {
    let energies = energy::parse_spec(energy_spec)?;
    let mu = db.mu_chantler(element, &energies, kind.into())?;

    let mut table = Table::new(&["energy_eV", "mu_cm2_per_g"]);
    for (e, m) in energies.iter().zip(mu.iter()) {
        table.row(vec![num(*e), num(*m)]);
    }
    table.write(out, format)
}

fn cmd_f0(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    ion: Option<&str>,
    q_spec: &str,
    list: bool,
) -> Result<()> {
    if list {
        let ions = db.f0_ions(ion)?;
        let mut table = Table::new(&["ion"]);
        for i in ions {
            table.row(vec![i.to_string()]);
        }
        return table.write(out, format);
    }

    let ion = ion.context("an ion name is required unless --list is given")?;
    let qs = energy::parse_spec(q_spec)?;
    let f0 = db.f0(ion, &qs)?;

    let mut table = Table::new(&["q_inv_angstrom", "f0"]);
    for (q, v) in qs.iter().zip(f0.iter()) {
        table.row(vec![num(*q), num(*v)]);
    }
    table.write(out, format)
}

fn cmd_delta_beta(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    formula: &str,
    energy_spec: &str,
    density: Option<f64>,
) -> Result<()> {
    let energies = energy::parse_spec(energy_spec)?;
    let (resolved_formula, resolved_density) = resolve_target(db, formula, density)?;

    let mut table = Table::new(&[
        "energy_eV",
        "delta",
        "beta",
        "atten_len_cm",
        "crit_angle_mrad",
    ]);
    for &e in &energies {
        let n = db.xray_delta_beta(&resolved_formula, resolved_density, e)?;
        // Critical angle for total external reflection: theta_c = sqrt(2*delta).
        let theta_c = if n.delta > 0.0 {
            (2.0 * n.delta).sqrt() * 1000.0
        } else {
            0.0
        };
        table.row(vec![
            num(e),
            num(n.delta),
            num(n.beta),
            num(n.attenuation_length_cm),
            num(theta_c),
        ]);
    }
    table.write(out, format)
}

fn cmd_material(db: &XrayDb, out: &mut impl Write, format: Format, name: &str) -> Result<()> {
    let material = db
        .find_material(name)
        .with_context(|| format!("'{name}' is not in the materials table"))?;

    #[derive(Serialize)]
    struct Rec<'a> {
        name: &'a str,
        formula: &'a str,
        density: f64,
    }
    let record = Rec {
        name: material.name,
        formula: material.formula,
        density: material.density,
    };
    let pairs = [
        ("name", material.name.to_string()),
        ("formula", material.formula.to_string()),
        ("density (g/cm^3)", num(material.density)),
    ];
    write_record(out, format, &pairs, &record)
}

fn cmd_materials(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    filter: Option<&str>,
) -> Result<()> {
    let needle = filter.map(str::to_lowercase);
    let mut table = Table::new(&["name", "formula", "density_g_per_cm3"]);
    for m in db.materials() {
        if let Some(ref needle) = needle
            && !m.name.to_lowercase().contains(needle)
            && !m.formula.to_lowercase().contains(needle)
        {
            continue;
        }
        table.row(vec![
            m.name.to_string(),
            m.formula.to_string(),
            num(m.density),
        ]);
    }
    table.write(out, format)
}

fn cmd_guess_edge(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    energy: f64,
    edges: &str,
    tolerance: Option<f64>,
) -> Result<()> {
    let labels: Vec<&str> = edges
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect();
    if labels.is_empty() {
        bail!("--edges must list at least one edge label");
    }

    let guess = db
        .guess_edge(energy, Some(&labels), tolerance)
        .with_context(|| match tolerance {
            Some(t) => format!("no {edges} edge within {t} eV of {energy} eV"),
            None => format!("no edge found near {energy} eV"),
        })?;

    #[derive(Serialize)]
    struct Rec {
        element: String,
        edge: String,
        energy_ev: f64,
        difference_ev: f64,
    }
    let record = Rec {
        element: guess.element.clone(),
        edge: guess.edge.clone(),
        energy_ev: guess.energy,
        difference_ev: guess.difference,
    };
    let pairs = [
        ("element", guess.element.clone()),
        ("edge", guess.edge.clone()),
        ("edge energy (eV)", num(guess.energy)),
        ("difference (eV)", num(guess.difference)),
    ];
    write_record(out, format, &pairs, &record)
}

fn cmd_core_width(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    element: &str,
    edge: Option<&str>,
) -> Result<()> {
    let widths = db.core_width(element, edge)?;
    let mut sorted: Vec<_> = widths.into_iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(&b.0));

    let mut table = Table::new(&["edge", "width_eV"]);
    for (label, width) in sorted {
        table.row(vec![label, num(width)]);
    }
    table.write(out, format)
}

#[allow(clippy::too_many_arguments)]
fn cmd_ionchamber(
    db: &XrayDb,
    out: &mut impl Write,
    format: Format,
    gases: &[String],
    volts: f64,
    length: f64,
    energy: f64,
    sensitivity: f64,
    with_compton: bool,
    both_carriers: bool,
) -> Result<()> {
    let mut parsed: Vec<(String, f64)> = Vec::with_capacity(gases.len());
    for spec in gases {
        match spec.split_once(':') {
            Some((name, fraction)) => {
                let f: f64 = fraction
                    .trim()
                    .parse()
                    .with_context(|| format!("invalid gas fraction in '{spec}'"))?;
                parsed.push((name.trim().to_string(), f));
            }
            None => parsed.push((spec.trim().to_string(), 1.0)),
        }
    }
    let refs: Vec<(&str, f64)> = parsed.iter().map(|(n, f)| (n.as_str(), *f)).collect();

    let fluxes = db.ionchamber_fluxes(
        &refs,
        volts,
        length,
        energy,
        sensitivity,
        with_compton,
        both_carriers,
    )?;

    #[derive(Serialize)]
    struct Rec {
        incident: f64,
        transmitted: f64,
        photo: f64,
        incoherent: f64,
        coherent: f64,
    }
    let record = Rec {
        incident: fluxes.incident,
        transmitted: fluxes.transmitted,
        photo: fluxes.photo,
        incoherent: fluxes.incoherent,
        coherent: fluxes.coherent,
    };
    let pairs = [
        ("incident flux (photons/s)", num(fluxes.incident)),
        ("transmitted (photons/s)", num(fluxes.transmitted)),
        ("photoabsorbed (photons/s)", num(fluxes.photo)),
        ("incoherent (photons/s)", num(fluxes.incoherent)),
        ("coherent (photons/s)", num(fluxes.coherent)),
    ];
    write_record(out, format, &pairs, &record)
}

fn cmd_compton(db: &XrayDb, out: &mut impl Write, format: Format, energy: f64) -> Result<()> {
    let c = db.compton_energies(energy);

    #[derive(Serialize)]
    struct Rec {
        incident: f64,
        xray_90deg: f64,
        xray_mean: f64,
        electron_mean: f64,
    }
    let record = Rec {
        incident: c.incident,
        xray_90deg: c.xray_90deg,
        xray_mean: c.xray_mean,
        electron_mean: c.electron_mean,
    };
    let pairs = [
        ("incident (eV)", num(c.incident)),
        ("scattered at 90 deg (eV)", num(c.xray_90deg)),
        ("mean scattered photon (eV)", num(c.xray_mean)),
        ("mean recoil electron (eV)", num(c.electron_mean)),
    ];
    write_record(out, format, &pairs, &record)
}

fn cmd_info(db: &XrayDb, out: &mut impl Write, format: Format) -> Result<()> {
    let raw = db.raw();

    #[derive(Serialize)]
    struct Rec<'a> {
        xraydb_version: Option<&'a str>,
        elements: usize,
        xray_levels: usize,
        xray_transitions: usize,
        chantler_elements: usize,
        waasmaier_ions: usize,
        materials: usize,
    }
    let record = Rec {
        xraydb_version: db.data_version(),
        elements: raw.elements.len(),
        xray_levels: raw.xray_levels.len(),
        xray_transitions: raw.xray_transitions.len(),
        chantler_elements: raw.chantler.len(),
        waasmaier_ions: raw.waasmaier.len(),
        materials: db.materials().len(),
    };
    let pairs = [
        (
            "XrayDB data version",
            db.data_version().unwrap_or("unknown").to_string(),
        ),
        ("elements", raw.elements.len().to_string()),
        ("absorption edges", raw.xray_levels.len().to_string()),
        ("emission lines", raw.xray_transitions.len().to_string()),
        (
            "elements with Chantler data",
            raw.chantler.len().to_string(),
        ),
        ("Waasmaier-Kirfel ions", raw.waasmaier.len().to_string()),
        ("materials", db.materials().len().to_string()),
    ];
    write_record(out, format, &pairs, &record)
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn format_shorthands_win_over_the_default() {
        let cli = Cli::parse_from(["xraydb", "--json", "element", "Fe"]);
        assert_eq!(cli.resolved_format(), Format::Json);

        let cli = Cli::parse_from(["xraydb", "--csv", "element", "Fe"]);
        assert_eq!(cli.resolved_format(), Format::Csv);

        let cli = Cli::parse_from(["xraydb", "element", "Fe"]);
        assert_eq!(cli.resolved_format(), Format::Text);
    }

    #[test]
    fn conflicting_format_flags_are_rejected() {
        assert!(Cli::try_parse_from(["xraydb", "--json", "--csv", "element", "Fe"]).is_err());
    }

    #[test]
    fn resolve_target_prefers_materials_then_elements() {
        let db = XrayDb::new();

        // A material name brings its own density.
        let (formula, density) = resolve_target(&db, "kapton", None).unwrap();
        assert_eq!(formula, "C22H10N2O5");
        assert_eq!(density, 1.42);

        // The materials table is consulted first, so "Fe" resolves through its
        // "iron" entry (7.88) rather than the elemental table (7.86). Both values
        // are legitimate; the materials table is the curated one.
        let (formula, density) = resolve_target(&db, "Fe", None).unwrap();
        assert_eq!(formula, "Fe");
        assert!((density - 7.88).abs() < 0.01, "got {density}");

        // An element with no materials-table entry falls back to elemental density.
        assert!(db.find_material("Mn").is_none());
        let (formula, density) = resolve_target(&db, "Mn", None).unwrap();
        assert_eq!(formula, "Mn");
        assert_eq!(density, db.density("Mn").unwrap());

        // An explicit density overrides.
        let (_, density) = resolve_target(&db, "kapton", Some(2.0)).unwrap();
        assert_eq!(density, 2.0);

        // A bare formula needs a density.
        assert!(resolve_target(&db, "Fe2O3", None).is_err());
        assert!(resolve_target(&db, "Fe2O3", Some(5.24)).is_ok());
    }
}
