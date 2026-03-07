//! Example: Particle-particle ladder resummation for the attractive Hubbard model.
//! Finds Tc via the Thouless criterion and computes the pp susceptibility
//! at several temperatures.

use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::{HubbardModel, ThermalParams};
use feynman_engine::resummation::bethe_salpeter::find_tc;
use feynman_engine::numerical::matsubara::compute_pp_susceptibility;

fn main() {
    // 1. Create an 8x8 lattice with attractive Hubbard model (t=1, U=-2)
    let lattice = SquareLattice::new(8, 8);
    let model = HubbardModel::new(lattice, 1.0, -2.0);
    println!("=== Ladder Resummation: Attractive Hubbard Model ===");
    println!("Lattice: 8x8 square, t = {}, U = {}", model.t, model.u);
    println!();

    // 2. Find Tc via Thouless criterion with n_matsubara=256
    println!("Searching for Tc via Thouless criterion (n_matsubara = 256)...");
    let tc = find_tc(&model, (0.01, 2.0), 256, 1e-4);
    println!("Tc = {:.6}", tc);
    println!();

    // 3. Compute chi_0 at several temperatures and print the Thouless criterion
    println!("=== Thouless Criterion: 1 - U * chi_0(q=0, nu=0) ===");
    println!("{:<12} {:<20} {:<20}", "T", "chi_0(q=0)", "1 - U*chi_0");
    println!("{}", "-".repeat(52));

    let temperatures = [
        tc * 0.8,
        tc * 0.9,
        tc * 0.95,
        tc * 1.0,
        tc * 1.05,
        tc * 1.1,
        tc * 1.5,
        tc * 2.0,
    ];

    for temp in &temperatures {
        let beta = 1.0 / temp;
        let params = ThermalParams::new(beta, 0.0, 256);
        let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);
        let thouless = 1.0 - model.u * chi0.re;
        println!("{:<12.6} {:<20.8} {:<20.8}", temp, chi0.re, thouless);
    }

    println!();
    println!("At T = Tc, the Thouless criterion 1 - U*chi_0 should vanish,");
    println!("signaling the onset of superconducting instability.");
}
