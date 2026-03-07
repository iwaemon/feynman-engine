//! Example: Generate and analyze 1st and 2nd order self-energy diagrams
//! for the attractive Hubbard model on a 4x4 square lattice.

use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::HubbardModel;
use feynman_engine::diagrams::generate::{generate_diagrams, Observable};
use feynman_engine::diagrams::classify::classify_diagrams;
use feynman_engine::diagrams::rules::apply_feynman_rules;
use feynman_engine::visualization::dot::to_dot_all;
use feynman_engine::visualization::json::to_json_all;

fn main() {
    // 1. Create a 4x4 lattice with attractive Hubbard model (t=1, U=-2)
    let lattice = SquareLattice::new(4, 4);
    let model = HubbardModel::new(lattice, 1.0, -2.0);
    println!("=== Attractive Hubbard Model ===");
    println!("Lattice: 4x4 square, t = {}, U = {}", model.t, model.u);
    println!();

    // 2. Generate 1st and 2nd order self-energy diagrams
    let diagrams_1st = generate_diagrams(&model, 1, Observable::SelfEnergy);
    let diagrams_2nd = generate_diagrams(&model, 2, Observable::SelfEnergy);
    println!("1st order: {} Wick contractions", diagrams_1st.len());
    println!("2nd order: {} Wick contractions", diagrams_2nd.len());
    println!();

    // 3. Classify diagrams by topology
    let classified_1st = classify_diagrams(&diagrams_1st);
    let classified_2nd = classify_diagrams(&diagrams_2nd);
    println!("1st order: {} unique topologies", classified_1st.len());
    println!("2nd order: {} unique topologies", classified_2nd.len());
    println!();

    // 4. Apply Feynman rules and print results
    println!("=== 1st Order Self-Energy Diagrams ===");
    for (i, (diag, weight)) in classified_1st.iter().enumerate() {
        let expr = apply_feynman_rules(diag, &model);
        println!("Topology {} (weight = {}):", i + 1, weight);
        println!("  {}", expr.description);
        println!("  Prefactor: {:.4}", expr.prefactor);
        println!("  Expression: {}", expr.expr);
        println!();
    }

    println!("=== 2nd Order Self-Energy Diagrams ===");
    for (i, (diag, weight)) in classified_2nd.iter().enumerate() {
        let expr = apply_feynman_rules(diag, &model);
        println!("Topology {} (weight = {}):", i + 1, weight);
        println!("  {}", expr.description);
        println!("  Prefactor: {:.4}", expr.prefactor);
        println!("  Expression: {}", expr.expr);
        println!();
    }

    // 5. Output DOT visualization to stdout
    let mut all_classified = classified_1st;
    all_classified.extend(classified_2nd);
    let dot = to_dot_all(&all_classified);
    println!("=== DOT Visualization ===");
    println!("{}", dot);

    // 6. Output JSON for D3.js visualization
    let json = to_json_all(&all_classified);
    println!("=== JSON Data ===");
    println!("{}", json);
}
