mod chantler;
mod elam;
mod parsers;
mod waasmaier;

use std::io::Write;
use std::path::Path;

use xraydb_data::XrayDatabase;

fn main() {
    let data_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("XrayDB")
        .join("data_sources");

    if !data_dir.exists() {
        eprintln!(
            "Error: XrayDB data_sources directory not found at {:?}",
            data_dir
        );
        eprintln!("Clone the upstream repo: git clone https://github.com/xraypy/XrayDB.git XrayDB");
        std::process::exit(1);
    }

    println!("Parsing raw data files from {:?}...", data_dir);

    let version = parsers::parse_version(&data_dir.join("Version.dat"));
    println!("  Version: {} entries", version.len());

    let elements = parsers::parse_elements(&data_dir.join("elemental_data.txt"));
    println!("  Elements: {} entries", elements.len());

    let ionization_potentials =
        parsers::parse_ionization_potentials(&data_dir.join("ion_chamber_potentials.txt"));
    println!(
        "  Ionization potentials: {} entries",
        ionization_potentials.len()
    );

    let compton = parsers::parse_compton_energies(&data_dir.join("Compton_energies.txt"));
    println!("  Compton energies: {} points", compton.incident.len());

    let (kk, ko, corelevel) = parsers::parse_corehole_data(
        &data_dir.join("keskirahkonen_krause.dat"),
        &data_dir.join("krause_oliver1979.dat"),
    );
    println!(
        "  Core widths: KK={}, KO={}, merged={}",
        kk.len(),
        ko.len(),
        corelevel.len()
    );

    let waasmaier = waasmaier::parse_waasmaier(&data_dir.join("waasmaeir_kirfel.dat"));
    println!("  Waasmaier: {} ions", waasmaier.len());

    let (xray_levels, xray_transitions, coster_kronig, photoabsorption, scattering) =
        elam::parse_elam(&data_dir.join("elam.dat"));
    println!(
        "  Elam: {} levels, {} transitions, {} CK, {} photo, {} scatter",
        xray_levels.len(),
        xray_transitions.len(),
        coster_kronig.len(),
        photoabsorption.len(),
        scattering.len(),
    );

    let chantler = chantler::parse_chantler(&data_dir.join("chantler").join("fine"));
    println!("  Chantler: {} elements", chantler.len());

    let db = XrayDatabase {
        version,
        elements,
        xray_levels,
        xray_transitions,
        coster_kronig,
        photoabsorption,
        scattering,
        chantler,
        waasmaier,
        compton_energies: compton,
        keski_rahkonen_krause: kk,
        krause_oliver: ko,
        corelevel_widths: corelevel,
        ionization_potentials,
    };

    println!("\nSerializing with postcard...");
    let serialized = postcard::to_allocvec(&db).expect("postcard serialization failed");
    println!(
        "  Serialized size: {} bytes ({:.2} MB)",
        serialized.len(),
        serialized.len() as f64 / 1_048_576.0
    );

    println!("Compressing with zstd (level 19)...");
    let compressed = zstd::encode_all(&serialized[..], 19).expect("zstd compression failed");
    println!(
        "  Compressed size: {} bytes ({:.2} MB)",
        compressed.len(),
        compressed.len() as f64 / 1_048_576.0
    );
    println!(
        "  Compression ratio: {:.1}x",
        serialized.len() as f64 / compressed.len() as f64
    );

    let out_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("xraydb-lib")
        .join("data")
        .join("xraydb.bin.zst");

    std::fs::create_dir_all(out_path.parent().unwrap()).unwrap();
    let mut f = std::fs::File::create(&out_path).expect("failed to create output file");
    f.write_all(&compressed)
        .expect("failed to write compressed data");

    println!("\nWrote {:?}", out_path);

    // Verify round-trip
    println!("Verifying round-trip deserialization...");
    let decompressed = zstd::decode_all(&compressed[..]).expect("zstd decompression failed");
    assert_eq!(decompressed.len(), serialized.len());
    let _db2: XrayDatabase =
        postcard::from_bytes(&decompressed).expect("postcard deserialization failed");
    println!("  Round-trip OK!");
}
