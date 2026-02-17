use xraydb::XrayDb;

#[cfg(feature = "optics")]
mod optics {
    use super::*;
    use xraydb::Polarization;

    #[test]
    fn test_darwin_width_si_111() {
        let db = XrayDb::new();
        let dw = db
            .darwin_width(10000.0, "Si", (1, 1, 1), None, Polarization::S, false, false, 1)
            .unwrap();
        let dw = dw.expect("Bragg condition should be satisfied at 10 keV");

        // Bragg angle for Si(111) at 10 keV should be around 11.4 degrees
        let theta_deg = dw.theta.to_degrees();
        assert!(
            theta_deg > 10.0 && theta_deg < 13.0,
            "theta = {theta_deg} deg"
        );

        // Energy width should be a few eV
        assert!(dw.energy_width > 0.1, "energy_width = {}", dw.energy_width);
        assert!(
            dw.energy_width < 100.0,
            "energy_width = {}",
            dw.energy_width
        );

        // Theta width should be positive
        assert!(dw.theta_width > 0.0);
        assert!(dw.theta_fwhm > 0.0);

        // Intensity should have values near 1.0 (total reflection)
        let max_intensity = dw.intensity.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            max_intensity > 0.8,
            "max intensity = {max_intensity}"
        );

        // Rocking curve should have meaningful FWHM
        assert!(dw.rocking_energy_fwhm > 0.0);
    }

    #[test]
    fn test_darwin_width_si_220() {
        let db = XrayDb::new();
        let dw = db
            .darwin_width(10000.0, "Si", (2, 2, 0), None, Polarization::S, false, false, 1)
            .unwrap();
        let dw = dw.expect("Bragg condition should be satisfied");

        // Si(220) at 10 keV: larger angle than (111)
        let dw_111 = db
            .darwin_width(10000.0, "Si", (1, 1, 1), None, Polarization::S, false, false, 1)
            .unwrap()
            .unwrap();
        assert!(dw.theta > dw_111.theta);

        // Energy width should be narrower for higher-order reflection
        assert!(dw.energy_width < dw_111.energy_width);
    }

    #[test]
    fn test_darwin_width_bragg_impossible() {
        let db = XrayDb::new();
        // Very low energy should make Bragg condition impossible
        let result = db
            .darwin_width(100.0, "Si", (1, 1, 1), None, Polarization::S, false, false, 1)
            .unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_darwin_width_p_polarization() {
        let db = XrayDb::new();
        let dw_s = db
            .darwin_width(10000.0, "Si", (1, 1, 1), None, Polarization::S, false, false, 1)
            .unwrap()
            .unwrap();
        let dw_p = db
            .darwin_width(10000.0, "Si", (1, 1, 1), None, Polarization::P, false, false, 1)
            .unwrap()
            .unwrap();

        // P-polarization width should be narrower than S
        assert!(dw_p.energy_width < dw_s.energy_width);
    }

    #[test]
    fn test_darwin_width_ge() {
        let db = XrayDb::new();
        let dw = db
            .darwin_width(10000.0, "Ge", (1, 1, 1), None, Polarization::S, false, false, 1)
            .unwrap();
        assert!(dw.is_some());
    }

    #[test]
    fn test_mirror_reflectivity_si() {
        let db = XrayDb::new();
        // Silicon mirror at 10 keV
        let theta: Vec<f64> = (1..100)
            .map(|i| i as f64 * 0.1e-3) // 0.1 to 10 mrad
            .collect();

        let refl = db
            .mirror_reflectivity("Si", &theta, 10000.0, 2.33, 0.0, Polarization::S)
            .unwrap();

        assert_eq!(refl.len(), theta.len());

        // At very small angles, reflectivity should be near 1
        assert!(refl[0] > 0.9, "R at 0.1 mrad = {}", refl[0]);

        // At large angles, reflectivity should be small
        let last = refl.len() - 1;
        assert!(refl[last] < 0.1, "R at 10 mrad = {}", refl[last]);

        // Should be monotonically decreasing overall
        // (may not be perfectly monotonic due to absorption, but generally true)
    }

    #[test]
    fn test_mirror_reflectivity_pt() {
        let db = XrayDb::new();
        let theta: Vec<f64> = (1..50).map(|i| i as f64 * 0.5e-3).collect();

        let refl = db
            .mirror_reflectivity("Pt", &theta, 10000.0, 21.45, 0.0, Polarization::S)
            .unwrap();

        // Pt has higher critical angle than Si
        let refl_si = db
            .mirror_reflectivity("Si", &theta, 10000.0, 2.33, 0.0, Polarization::S)
            .unwrap();

        // At intermediate angles, Pt should reflect more than Si
        let mid = theta.len() / 2;
        assert!(
            refl[mid] > refl_si[mid],
            "Pt R = {}, Si R = {} at theta = {}",
            refl[mid],
            refl_si[mid],
            theta[mid]
        );
    }

    #[test]
    fn test_mirror_reflectivity_roughness() {
        let db = XrayDb::new();
        let theta: Vec<f64> = (1..50).map(|i| i as f64 * 0.2e-3).collect();

        let smooth = db
            .mirror_reflectivity("Si", &theta, 10000.0, 2.33, 0.0, Polarization::S)
            .unwrap();
        let rough = db
            .mirror_reflectivity("Si", &theta, 10000.0, 2.33, 5.0, Polarization::S)
            .unwrap();

        // Roughness should reduce reflectivity
        for i in 0..smooth.len() {
            assert!(
                rough[i] <= smooth[i] + 1e-10,
                "rough > smooth at i={}: {} vs {}",
                i,
                rough[i],
                smooth[i]
            );
        }
    }

    #[test]
    fn test_multilayer_si_w() {
        let db = XrayDb::new();
        let theta: Vec<f64> = (1..100).map(|i| i as f64 * 0.1e-3).collect();

        // Simple W/Si bilayer, repeated 20 times on Si substrate
        let refl = db
            .multilayer_reflectivity(
                &["W", "Si"],
                &[20.0, 20.0], // 20 Å each
                "Si",
                &theta,
                10000.0,
                20,
                &[19.25, 2.33],
                2.33,
                0.0,
                0.0,
                Polarization::S,
            )
            .unwrap();

        assert_eq!(refl.len(), theta.len());

        // Should have a Bragg peak somewhere
        let max_refl = refl.iter().cloned().fold(0.0_f64, f64::max);
        assert!(
            max_refl > 0.01,
            "max multilayer reflectivity = {max_refl}"
        );
    }

    #[test]
    fn test_coated_reflectivity_rh() {
        let db = XrayDb::new();
        let theta: Vec<f64> = (1..50).map(|i| i as f64 * 0.2e-3).collect();

        let refl = db
            .coated_reflectivity(
                "Rh",
                500.0, // 500 Å coating
                "Si",
                &theta,
                10000.0,
                12.41, // Rh density
                0.0,
                2.33, // Si density
                0.0,
                None,
                Polarization::S,
            )
            .unwrap();

        assert_eq!(refl.len(), theta.len());

        // At small angles, should have good reflectivity
        assert!(refl[0] > 0.5, "R at small angle = {}", refl[0]);
    }
}

#[test]
fn test_find_material() {
    let db = XrayDb::new();

    let (formula, density) = db.find_material("water").unwrap();
    assert_eq!(formula, "H2O");
    assert!((density - 1.0).abs() < 1e-6);

    let (formula, density) = db.find_material("silicon").unwrap();
    assert_eq!(formula, "Si");
    assert!((density - 2.329).abs() < 0.01);

    let (formula, density) = db.find_material("kapton").unwrap();
    assert_eq!(formula, "C22H10N2O5");
    assert!((density - 1.42).abs() < 0.01);

    // Case insensitive
    assert!(db.find_material("Water").is_some());
    assert!(db.find_material("SILICON").is_some());

    // Unknown material
    assert!(db.find_material("unobtainium").is_none());
}

#[test]
fn test_material_mu_named() {
    let db = XrayDb::new();
    let mu = db
        .material_mu_named("water", &[10000.0], xraydb::CrossSectionKind::Total, None)
        .unwrap();
    assert!(mu[0] > 0.0);
    assert!(mu[0] < 100.0);
}

#[test]
fn test_ionchamber_nitrogen() {
    let db = XrayDb::new();
    let fluxes = db
        .ionchamber_fluxes(
            &[("nitrogen", 1.0)],
            1.0,   // 1 V
            100.0, // 100 cm
            10000.0,
            1e-6,  // 1 µA/V
            true,
            true,
        )
        .unwrap();

    // Incident flux should be positive
    assert!(fluxes.incident > 0.0, "incident = {}", fluxes.incident);
    // Transmitted should be less than incident
    assert!(fluxes.transmitted < fluxes.incident);
    // Photo absorption flux should be positive
    assert!(fluxes.photo > 0.0);
}

#[test]
fn test_ionchamber_argon() {
    let db = XrayDb::new();
    let fluxes = db
        .ionchamber_fluxes(
            &[("argon", 1.0)],
            1.0,
            30.0, // 30 cm
            10000.0,
            1e-6,
            true,
            true,
        )
        .unwrap();

    assert!(fluxes.incident > 0.0);

    // Argon absorbs more than nitrogen at same conditions
    let fluxes_n = db
        .ionchamber_fluxes(
            &[("nitrogen", 1.0)],
            1.0,
            30.0,
            10000.0,
            1e-6,
            true,
            true,
        )
        .unwrap();

    // Higher absorption means lower flux for same voltage
    // (more photons needed to produce same signal in N2 which absorbs less)
    assert!(
        fluxes.incident < fluxes_n.incident,
        "Ar incident = {}, N2 incident = {}",
        fluxes.incident,
        fluxes_n.incident
    );
}
