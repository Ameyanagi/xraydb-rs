use criterion::{Criterion, black_box, criterion_group, criterion_main};
use xraydb::{ChantlerKind, CrossSectionKind, XrayDb};

fn bench_lookup_apis(c: &mut Criterion) {
    let db = XrayDb::new();
    let ids = ["Fe", "iron", "26", "Si", "silicon", "14"];

    c.bench_function("atomic_number_mixed_identifiers", |b| {
        b.iter(|| {
            for id in ids {
                black_box(db.atomic_number(black_box(id)).unwrap());
            }
        });
    });

    c.bench_function("symbol_mixed_identifiers", |b| {
        b.iter(|| {
            for id in ids {
                black_box(db.symbol(black_box(id)).unwrap());
            }
        });
    });

    c.bench_function("molar_mass_mixed_identifiers", |b| {
        b.iter(|| {
            for id in ids {
                black_box(db.molar_mass(black_box(id)).unwrap());
            }
        });
    });
}

fn bench_mu_elam_vector(c: &mut Criterion) {
    let db = XrayDb::new();
    let energies: Vec<f64> = (0..200).map(|i| 1000.0 + i as f64 * 150.0).collect();

    c.bench_function("mu_elam_fe_vector_total", |b| {
        b.iter(|| {
            black_box(
                db.mu_elam(
                    black_box("Fe"),
                    black_box(&energies),
                    black_box(CrossSectionKind::Total),
                )
                .unwrap(),
            );
        });
    });
}

fn bench_f2_chantler_vector(c: &mut Criterion) {
    let db = XrayDb::new();
    let energies: Vec<f64> = (0..200).map(|i| 1000.0 + i as f64 * 200.0).collect();

    c.bench_function("f2_chantler_fe_vector", |b| {
        b.iter(|| {
            black_box(
                db.f2_chantler(black_box("Fe"), black_box(&energies))
                    .unwrap(),
            );
        });
    });

    c.bench_function("mu_chantler_fe_vector_total", |b| {
        b.iter(|| {
            black_box(
                db.mu_chantler(
                    black_box("Fe"),
                    black_box(&energies),
                    black_box(ChantlerKind::Total),
                )
                .unwrap(),
            );
        });
    });
}

criterion_group!(
    benches,
    bench_lookup_apis,
    bench_mu_elam_vector,
    bench_f2_chantler_vector
);
criterion_main!(benches);
