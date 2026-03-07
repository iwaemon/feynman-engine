use crate::models::lattice::SquareLattice;

#[derive(Debug, Clone)]
pub struct HubbardModel {
    pub lattice: SquareLattice,
    pub t: f64,
    pub u: f64,
}

impl HubbardModel {
    pub fn new(lattice: SquareLattice, t: f64, u: f64) -> Self {
        Self { lattice, t, u }
    }

    pub fn dispersion(&self, k: &[f64; 2]) -> f64 {
        self.lattice.dispersion(k, self.t)
    }
}

#[derive(Debug, Clone)]
pub struct ThermalParams {
    pub beta: f64,
    pub mu: f64,
    pub n_matsubara: usize,
}

impl ThermalParams {
    pub fn new(beta: f64, mu: f64, n_matsubara: usize) -> Self {
        Self { beta, mu, n_matsubara }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;

    #[test]
    fn test_hubbard_model() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        assert_eq!(model.t, 1.0);
        assert_eq!(model.u, -2.0);
    }

    #[test]
    fn test_thermal_params() {
        let params = ThermalParams::new(10.0, 0.0, 128);
        assert_eq!(params.beta, 10.0);
        assert_eq!(params.mu, 0.0);
        assert_eq!(params.n_matsubara, 128);
    }
}
