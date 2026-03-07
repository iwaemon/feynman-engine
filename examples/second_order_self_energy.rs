//! Example: Generate and analyze 1st and 2nd order self-energy diagrams
//! for the attractive Hubbard model on a 4x4 square lattice.
//!
//! By default, prints a text summary to stdout. Use flags to export:
//!   --dot <file>   Write Graphviz DOT to file (requires `dot` from Graphviz to render)
//!   --json <file>  Write JSON (for D3.js visualization) to file
//!   --dot -        Write DOT to stdout (suppresses text summary)
//!   --json -       Write JSON to stdout (suppresses text summary)

use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::HubbardModel;
use feynman_engine::diagrams::generate::{generate_diagrams, Observable};
use feynman_engine::diagrams::classify::classify_diagrams;
use feynman_engine::diagrams::rules::apply_feynman_rules;
use feynman_engine::visualization::dot::to_dot_all;
use feynman_engine::visualization::json::to_json_all;
use std::fs;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut dot_output: Option<String> = None;
    let mut json_output: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--dot" => {
                i += 1;
                dot_output = Some(args.get(i).cloned().unwrap_or_else(|| "-".to_string()));
            }
            "--json" => {
                i += 1;
                json_output = Some(args.get(i).cloned().unwrap_or_else(|| "-".to_string()));
            }
            _ => {
                eprintln!("Unknown option: {}", args[i]);
                eprintln!("Usage: second_order_self_energy [--dot <file|->] [--json <file|->]");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    // 1. Create a 4x4 lattice with attractive Hubbard model (t=1, U=-2)
    let lattice = SquareLattice::new(4, 4);
    let model = HubbardModel::new(lattice, 1.0, -2.0);

    // 2. Generate 1st and 2nd order self-energy diagrams
    let diagrams_1st = generate_diagrams(&model, 1, Observable::SelfEnergy);
    let diagrams_2nd = generate_diagrams(&model, 2, Observable::SelfEnergy);

    // 3. Classify diagrams by topology
    let classified_1st = classify_diagrams(&diagrams_1st);
    let classified_2nd = classify_diagrams(&diagrams_2nd);

    // 4. Combine all classified diagrams
    let mut all_classified = classified_1st.clone();
    all_classified.extend(classified_2nd.clone());

    // 5. Handle DOT output
    if let Some(ref path) = dot_output {
        let dot = to_dot_all(&all_classified);
        if path == "-" {
            println!("{}", dot);
        } else {
            fs::write(path, &dot).unwrap_or_else(|e| {
                eprintln!("Error writing DOT to {}: {}", path, e);
                std::process::exit(1);
            });
            eprintln!("DOT written to {}", path);
        }
    }

    // 6. Handle JSON output
    if let Some(ref path) = json_output {
        let json = to_json_all(&all_classified);
        if path == "-" {
            println!("{}", json);
        } else {
            fs::write(path, &json).unwrap_or_else(|e| {
                eprintln!("Error writing JSON to {}: {}", path, e);
                std::process::exit(1);
            });
            eprintln!("JSON written to {}", path);
        }
    }

    // 7. Print text summary (unless stdout is used for DOT/JSON)
    let stdout_used = matches!(dot_output.as_deref(), Some("-"))
        || matches!(json_output.as_deref(), Some("-"));

    if !stdout_used {
        println!("=== Attractive Hubbard Model ===");
        println!("Lattice: 4x4 square, t = {}, U = {}", model.t, model.u);
        println!();
        println!("1st order: {} Wick contractions", diagrams_1st.len());
        println!("2nd order: {} Wick contractions", diagrams_2nd.len());
        println!();
        println!("1st order: {} unique topologies", classified_1st.len());
        println!("2nd order: {} unique topologies", classified_2nd.len());
        println!();

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

        if dot_output.is_none() && json_output.is_none() {
            eprintln!("Hint: use --dot diagrams.dot to export Graphviz, --json data.json for D3.js");
        }
    }
}
