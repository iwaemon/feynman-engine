use num_complex::Complex64;
use crate::models::hubbard::{HubbardModel, ThermalParams};
use crate::numerical::matsubara::compute_pp_susceptibility;

/// T-matrix from pp-ladder resummation: T(q,iνₘ) = U / (1 - U·χ₀(q,iνₘ))
pub fn solve_tmatrix(
    model: &HubbardModel,
    params: &ThermalParams,
    q: &[f64; 2],
    m: i32,
) -> Complex64 {
    let chi0 = compute_pp_susceptibility(model, params, q, m);
    model.u / (1.0 - model.u * chi0)
}

/// Find superconducting Tc via Thouless criterion: 1 - U·χ₀(q=0, 0; T) = 0
/// Uses bisection method.
pub fn find_tc(
    model: &HubbardModel,
    t_range: (f64, f64),
    n_matsubara: usize,
    tol: f64,
) -> f64 {
    let (mut t_lo, mut t_hi) = t_range;

    let criterion = |temp: f64| -> f64 {
        let params = ThermalParams::new(1.0 / temp, 0.0, n_matsubara);
        let chi0 = compute_pp_susceptibility(model, &params, &[0.0, 0.0], 0);
        (1.0 - model.u * chi0).re
    };

    let mut f_lo = criterion(t_lo);
    let f_hi = criterion(t_hi);

    assert!(
        f_lo * f_hi < 0.0,
        "Thouless criterion not bracketed: f({})={}, f({})={}",
        t_lo, f_lo, t_hi, f_hi
    );

    while (t_hi - t_lo) > tol {
        let t_mid = (t_lo + t_hi) / 2.0;
        let f_mid = criterion(t_mid);
        if f_mid * f_lo < 0.0 {
            t_hi = t_mid;
        } else {
            t_lo = t_mid;
            f_lo = f_mid;
        }
    }

    (t_lo + t_hi) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;

    #[test]
    fn test_tmatrix() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let params = ThermalParams::new(5.0, 0.0, 512);

        let t_val = solve_tmatrix(&model, &params, &[0.0, 0.0], 0);
        let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);
        let expected = model.u / (1.0 - model.u * chi0);

        assert!((t_val - expected).norm() < 1e-10);
    }

    #[test]
    fn test_find_tc() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);

        let tc = find_tc(&model, (0.01, 2.0), 512, 1e-4);
        assert!(tc > 0.0, "tc = {}", tc);
        assert!(tc < 4.0, "tc = {}", tc);

        // Verify bracket around Tc
        let params_above = ThermalParams::new(1.0 / (tc * 1.05), 0.0, 512);
        let params_below = ThermalParams::new(1.0 / (tc * 0.95), 0.0, 512);
        let chi_above = compute_pp_susceptibility(&model, &params_above, &[0.0, 0.0], 0);
        let chi_below = compute_pp_susceptibility(&model, &params_below, &[0.0, 0.0], 0);

        assert!((1.0 - model.u * chi_above).re > 0.0);
        assert!((1.0 - model.u * chi_below).re < 0.0);
    }

    #[test]
    fn test_repulsive_no_instability() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, 2.0); // repulsive
        let params = ThermalParams::new(20.0, 0.0, 512);

        let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);
        assert!((1.0 - model.u * chi0).re > 0.0);
    }
}
