//! Benchmarks covering the paths optimised in 0.2.0, so the gains stay guarded.
//!
//! Reference numbers on an M-series Mac (0.1.2 → 0.2.0):
//!
//! | Benchmark | 0.1.2 | 0.2.0 |
//! |---|---|---|
//! | `chantler/f2_scalar` | 19.7 µs | ~120 ns |
//! | `materials/xray_delta_beta_sio2` | 142 µs | ~1.1 µs |
//! | `transitions/xray_lines_u_excitation` | 84 µs | ~2.5 µs |
//! | `transitions/xray_lines_fe` | 13.5 µs | ~1.2 µs |

use criterion::{Criterion, black_box, criterion_group, criterion_main};
use xraydb::{ChantlerKind, CrossSectionKind, XrayDb};

fn energies(n: usize) -> Vec<f64> {
    (0..n).map(|i| 1000.0 + i as f64 * 150.0).collect()
}

fn bench_element_lookup(c: &mut Criterion) {
    let db = XrayDb::new();
    let ids = ["Fe", "iron", "26", "Si", "silicon", "14"];

    let mut group = c.benchmark_group("elements");
    group.bench_function("atomic_number_mixed_identifiers", |b| {
        b.iter(|| {
            for id in ids {
                black_box(db.atomic_number(black_box(id)).unwrap());
            }
        });
    });
    group.bench_function("symbol_mixed_identifiers", |b| {
        b.iter(|| {
            for id in ids {
                black_box(db.symbol(black_box(id)).unwrap());
            }
        });
    });
    group.bench_function("molar_mass_mixed_identifiers", |b| {
        b.iter(|| {
            for id in ids {
                black_box(db.molar_mass(black_box(id)).unwrap());
            }
        });
    });
    group.finish();
}

fn bench_elam(c: &mut Criterion) {
    let db = XrayDb::new();
    let batch = energies(200);

    let mut group = c.benchmark_group("elam");
    group.bench_function("mu_elam_scalar_total", |b| {
        b.iter(|| {
            black_box(
                db.mu_elam_at(
                    black_box("Fe"),
                    black_box(10_000.0),
                    CrossSectionKind::Total,
                )
                .unwrap(),
            );
        });
    });
    group.bench_function("mu_elam_batch200_total", |b| {
        b.iter(|| {
            black_box(
                db.mu_elam(black_box("Fe"), black_box(&batch), CrossSectionKind::Total)
                    .unwrap(),
            );
        });
    });
    group.bench_function("mu_elam_batch200_photo", |b| {
        b.iter(|| {
            black_box(
                db.mu_elam(black_box("Fe"), black_box(&batch), CrossSectionKind::Photo)
                    .unwrap(),
            );
        });
    });
    group.finish();
}

/// The headline optimisation: log-space Chantler tables are precomputed at load
/// time instead of being rebuilt on every call.
fn bench_chantler(c: &mut Criterion) {
    let db = XrayDb::new();
    let batch = energies(200);
    let big = energies(1000);

    let mut group = c.benchmark_group("chantler");
    group.bench_function("f2_scalar", |b| {
        b.iter(|| {
            black_box(
                db.f2_chantler_at(black_box("Fe"), black_box(10_000.0))
                    .unwrap(),
            )
        });
    });
    group.bench_function("f1_scalar", |b| {
        b.iter(|| {
            black_box(
                db.f1_chantler_at(black_box("Fe"), black_box(10_000.0))
                    .unwrap(),
            )
        });
    });
    group.bench_function("f2_batch200", |b| {
        b.iter(|| black_box(db.f2_chantler(black_box("Fe"), black_box(&batch)).unwrap()));
    });
    group.bench_function("f2_batch1000", |b| {
        b.iter(|| black_box(db.f2_chantler(black_box("Fe"), black_box(&big)).unwrap()));
    });
    group.bench_function("mu_chantler_batch200_total", |b| {
        b.iter(|| {
            black_box(
                db.mu_chantler(black_box("Fe"), black_box(&batch), ChantlerKind::Total)
                    .unwrap(),
            );
        });
    });
    group.finish();
}

/// Previously full-table linear scans; now index lookups.
fn bench_indexed_lookups(c: &mut Criterion) {
    let db = XrayDb::new();

    let mut group = c.benchmark_group("transitions");
    group.bench_function("xray_edge_fe_k", |b| {
        b.iter(|| black_box(db.xray_edge(black_box("Fe"), black_box("K")).unwrap()));
    });
    group.bench_function("xray_edges_fe", |b| {
        b.iter(|| black_box(db.xray_edges(black_box("Fe")).unwrap()));
    });
    group.bench_function("xray_lines_fe", |b| {
        b.iter(|| black_box(db.xray_lines(black_box("Fe"), None, None).unwrap()));
    });
    // Worst case in 0.1.2: a nested 1807 x 1430 scan.
    group.bench_function("xray_lines_u_excitation", |b| {
        b.iter(|| {
            black_box(
                db.xray_lines(black_box("U"), None, black_box(Some(20_000.0)))
                    .unwrap(),
            );
        });
    });
    group.bench_function("guess_edge", |b| {
        b.iter(|| black_box(db.guess_edge(black_box(7112.0), None, None)));
    });
    group.bench_function("core_width_fe_k", |b| {
        b.iter(|| black_box(db.core_width_at(black_box("Fe"), black_box("K")).unwrap()));
    });
    group.bench_function("ck_probability_u", |b| {
        b.iter(|| {
            black_box(
                db.ck_probability(black_box("U"), black_box("L1"), black_box("L3"), true)
                    .unwrap(),
            );
        });
    });
    group.bench_function("ionization_potential_argon", |b| {
        b.iter(|| black_box(db.ionization_potential(black_box("argon")).unwrap()));
    });
    group.finish();
}

fn bench_materials(c: &mut Criterion) {
    let db = XrayDb::new();
    let batch = energies(200);

    let mut group = c.benchmark_group("materials");
    group.bench_function("material_mu_at_sio2", |b| {
        b.iter(|| {
            black_box(
                db.material_mu_at(black_box("SiO2"), 2.2, 10_000.0, CrossSectionKind::Total)
                    .unwrap(),
            );
        });
    });
    group.bench_function("material_mu_batch200_kapton", |b| {
        b.iter(|| {
            black_box(
                db.material_mu(
                    black_box("C22H10N2O5"),
                    1.42,
                    black_box(&batch),
                    CrossSectionKind::Total,
                )
                .unwrap(),
            );
        });
    });
    group.bench_function("xray_delta_beta_sio2", |b| {
        b.iter(|| {
            black_box(
                db.xray_delta_beta(black_box("SiO2"), 2.2, 10_000.0)
                    .unwrap(),
            )
        });
    });
    group.bench_function("parse_formula_nested", |b| {
        b.iter(|| black_box(db.parse_formula(black_box("Mn(SO4)2(H2O)7")).unwrap()));
    });
    group.bench_function("ionchamber_fluxes_nitrogen", |b| {
        b.iter(|| {
            black_box(
                db.ionchamber_fluxes(&[("nitrogen", 1.0)], 1.0, 30.0, 10_000.0, 1e-6, true, true)
                    .unwrap(),
            );
        });
    });
    group.finish();
}

criterion_group!(
    benches,
    bench_element_lookup,
    bench_elam,
    bench_chantler,
    bench_indexed_lookups,
    bench_materials
);
criterion_main!(benches);
