use num_complex::Complex64;
use std::f64::consts::PI;
use crate::models::hubbard::{HubbardModel, ThermalParams};

/// Fermionic Matsubara frequency ωₙ = (2n+1)π/β
pub fn matsubara_freq(n: i32, beta: f64) -> f64 {
    (2 * n + 1) as f64 * PI / beta
}

/// Free Matsubara Green's function G₀(k, iωₙ) = 1/(iωₙ - εₖ + μ)
pub fn evaluate_g0(model: &HubbardModel, params: &ThermalParams,
                   k: &[f64; 2], n: i32) -> Complex64 {
    let omega = matsubara_freq(n, params.beta);
    let eps = model.dispersion(k);
    1.0 / Complex64::new(-eps + params.mu, omega)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;

    #[test]
    fn test_evaluate_g0() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let params = ThermalParams::new(10.0, 0.0, 128);

        let k = [0.0, 0.0]; // ε = -4
        let n = 0;
        let g = evaluate_g0(&model, &params, &k, n);

        let omega = PI / params.beta;
        let eps = model.dispersion(&k);
        let expected = 1.0 / Complex64::new(-eps + params.mu, omega);
        assert!((g - expected).norm() < 1e-12);
    }

    #[test]
    fn test_g0_sum_rule() {
        // (1/β) Σₙ G₀(k, iωₙ) → n_F(εₖ-μ)
        // At k=(π/2, π/2), ε≈0, μ=0: n_F(0) = 0.5
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, 0.0);
        let params = ThermalParams::new(20.0, 0.0, 2048);

        let k = [PI / 2.0, PI / 2.0];
        let eps = model.dispersion(&k);
        assert!(eps.abs() < 1e-10);

        let mut sum = 0.0;
        let n_max = params.n_matsubara as i32;
        for n in -n_max..n_max {
            let g = evaluate_g0(&model, &params, &k, n);
            sum += g.re / params.beta;
        }
        // (1/β) Σₙ Re G₀(k, iωₙ) = n_F(εₖ-μ) - 1/2
        // so n_F = 1/2 + (1/β) Σₙ Re G₀
        let nk = 0.5 + sum;
        assert!((nk - 0.5).abs() < 0.01, "nk = {}, expected 0.5", nk);
    }
}
