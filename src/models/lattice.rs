use std::f64::consts::PI;

#[derive(Debug, Clone)]
pub struct SquareLattice {
    pub lx: usize,
    pub ly: usize,
}

impl SquareLattice {
    pub fn new(lx: usize, ly: usize) -> Self {
        Self { lx, ly }
    }

    pub fn num_sites(&self) -> usize {
        self.lx * self.ly
    }

    pub fn k_grid(&self) -> Vec<[f64; 2]> {
        let mut points = Vec::with_capacity(self.num_sites());
        for ix in 0..self.lx {
            for iy in 0..self.ly {
                let kx = 2.0 * PI * (ix as f64) / (self.lx as f64);
                let ky = 2.0 * PI * (iy as f64) / (self.ly as f64);
                points.push([kx, ky]);
            }
        }
        points
    }

    pub fn dispersion(&self, k: &[f64; 2], t: f64) -> f64 {
        -2.0 * t * (k[0].cos() + k[1].cos())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_square_lattice() {
        let lat = SquareLattice::new(4, 4);
        assert_eq!(lat.num_sites(), 16);
        assert_eq!(lat.k_grid().len(), 16);
    }

    #[test]
    fn test_dispersion() {
        let lat = SquareLattice::new(4, 4);
        assert!((lat.dispersion(&[0.0, 0.0], 1.0) - (-4.0)).abs() < 1e-12);
        assert!((lat.dispersion(&[PI, PI], 1.0) - 4.0).abs() < 1e-12);
    }
}
