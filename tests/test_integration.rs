use feynman_engine::diagrams::generate::{generate_diagrams, Observable};
use feynman_engine::diagrams::classify::classify_diagrams;
use feynman_engine::diagrams::rules::apply_feynman_rules;
use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::{HubbardModel, ThermalParams};
use feynman_engine::numerical::matsubara::compute_pp_susceptibility;
use feynman_engine::resummation::bethe_salpeter::{find_tc, solve_tmatrix};

/// End-to-end: Wick contractions -> diagram generation -> classification -> Feynman rules
#[test]
fn test_wick_to_diagrams_to_classification() {
    let lat = SquareLattice::new(4, 4);
    let model = HubbardModel::new(lat, 1.0, -2.0);

    // Generate 2nd order self-energy diagrams
    let diagrams = generate_diagrams(&model, 2, Observable::SelfEnergy);
    assert!(!diagrams.is_empty(), "Should generate at least one diagram");

    // Classify by topology
    let classified = classify_diagrams(&diagrams);
    assert!(
        !classified.is_empty(),
        "Classified diagrams must be non-empty"
    );
    assert!(
        classified.len() <= diagrams.len(),
        "Classification should group diagrams: {} unique vs {} total",
        classified.len(),
        diagrams.len()
    );

    // Apply Feynman rules to each unique topology
    for (diag, weight) in &classified {
        assert_ne!(*weight, 0, "Each topology must have nonzero weight");

        let expr = apply_feynman_rules(diag, &model);
        assert!(
            expr.n_loops > 0,
            "2nd order self-energy diagrams must have at least 1 loop, got 0"
        );
        assert_ne!(expr.prefactor, 0.0, "Prefactor must be nonzero");
        assert!(
            expr.description.contains("Order 2"),
            "Description should mention order 2"
        );

        // Expression should contain Green's functions
        let s = format!("{}", expr.expr);
        assert!(
            s.contains("G"),
            "Expression should contain Green's function symbols"
        );
    }
}

/// End-to-end: Ladder resummation to find Tc for attractive Hubbard model
#[test]
fn test_ladder_resummation_find_tc() {
    let lat = SquareLattice::new(8, 8);
    let model = HubbardModel::new(lat, 1.0, -2.0);

    // Find Tc with bisection
    let tc = find_tc(&model, (0.01, 2.0), 256, 1e-3);

    assert!(tc > 0.0, "Tc must be positive, got {}", tc);
    assert!(tc < 4.0, "Tc must be less than 4t, got {}", tc);

    // Verify T-matrix diverges near Tc: evaluate just above Tc
    let beta_near_tc = 1.0 / (tc * 1.01);
    let params_near = ThermalParams::new(beta_near_tc, 0.0, 256);
    let t_matrix = solve_tmatrix(&model, &params_near, &[0.0, 0.0], 0);

    // Near Tc the denominator 1 - U*chi0 -> 0, so |T| should be large
    assert!(
        t_matrix.norm() > 10.0,
        "T-matrix should be large near Tc, got |T| = {}",
        t_matrix.norm()
    );
}

/// Repulsive Hubbard (U > 0) should have no particle-particle instability
#[test]
fn test_repulsive_hubbard_no_pp_instability() {
    let lat = SquareLattice::new(8, 8);
    let model = HubbardModel::new(lat, 1.0, 2.0); // U = +2 (repulsive)

    // Check at low temperature (high beta)
    let params = ThermalParams::new(20.0, 0.0, 256);
    let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);

    // For repulsive U and negative chi0.re, 1 - U*chi0 should stay positive
    // (no Thouless criterion satisfied => no superconducting instability)
    let stoner = (1.0 - model.u * chi0).re;
    assert!(
        stoner > 0.0,
        "1 - U*chi0 must be positive for repulsive U, got {} (chi0 = {})",
        stoner,
        chi0
    );
}
