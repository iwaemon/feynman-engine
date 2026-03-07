use num_complex::Complex64;
use std::f64::consts::PI;
use crate::models::hubbard::{HubbardModel, ThermalParams};
use crate::numerical::greens_function::evaluate_g0;

/// Bosonic Matsubara frequency νₘ = 2mπ/β
pub fn bosonic_matsubara_freq(m: i32, beta: f64) -> f64 {
    2.0 * m as f64 * PI / beta
}

/// Particle-particle bubble:
/// χ₀(q, iνₘ) = -(1/Nβ) Σₖ Σₙ G₀(k, iωₙ) G₀(q-k, iνₘ - iωₙ)
pub fn compute_pp_susceptibility(
    model: &HubbardModel,
    params: &ThermalParams,
    q: &[f64; 2],
    m: i32,
) -> Complex64 {
    let kpoints = model.lattice.k_grid();
    let n_sites = kpoints.len() as f64;
    let n_max = params.n_matsubara as i32;

    let mut chi0 = Complex64::new(0.0, 0.0);

    for k in &kpoints {
        let qmk = [q[0] - k[0], q[1] - k[1]];
        for n in -n_max..n_max {
            let g1 = evaluate_g0(model, params, k, n);
            let g2 = evaluate_g0(model, params, &qmk, m - 1 - n);
            chi0 += g1 * g2;
        }
    }

    -chi0 / (n_sites * params.beta)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;

    #[test]
    fn test_pp_susceptibility_real_at_q0() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let params = ThermalParams::new(5.0, 0.0, 512);

        let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);
        assert!(chi0.im.abs() < 1e-6, "imaginary part = {}", chi0.im);
        assert!(chi0.re < 0.0, "real part = {}", chi0.re);
    }

    #[test]
    fn test_pp_susceptibility_increases_with_beta() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);

        let params_high_t = ThermalParams::new(2.0, 0.0, 256);
        let params_low_t = ThermalParams::new(10.0, 0.0, 1024);

        let chi_high = compute_pp_susceptibility(&model, &params_high_t, &[0.0, 0.0], 0).re.abs();
        let chi_low = compute_pp_susceptibility(&model, &params_low_t, &[0.0, 0.0], 0).re.abs();

        assert!(chi_low > chi_high, "chi_low={}, chi_high={}", chi_low, chi_high);
    }
}
