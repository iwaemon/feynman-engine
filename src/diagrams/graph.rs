use crate::algebra::operators::Spin;

#[derive(Debug, Clone, serde::Serialize)]
pub struct Vertex {
    pub id: usize,
    pub site: String,
    pub time: String,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct Propagator {
    pub from: usize,
    pub to: usize,
    pub spin: Spin,
    pub external: bool,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct FeynmanDiagram {
    pub order: usize,
    pub vertices: Vec<Vertex>,
    pub propagators: Vec<Propagator>,
    pub sign: i32,
    pub symmetry_factor: f64,
}

impl FeynmanDiagram {
    /// Count independent loops = internal_edges - vertices + 1 (connected)
    pub fn count_loops(&self) -> usize {
        let n_internal = self.propagators.iter().filter(|p| !p.external).count();
        let n_verts = self.vertices.len();
        if n_verts == 0 {
            return 0;
        }
        n_internal - n_verts + 1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_feynman_diagram() {
        let d = FeynmanDiagram {
            order: 1,
            vertices: vec![
                Vertex { id: 1, site: "i".into(), time: "τ1".into() },
                Vertex { id: 2, site: "j".into(), time: "τ2".into() },
            ],
            propagators: vec![
                Propagator { from: 1, to: 2, spin: Spin::Up, external: false },
                Propagator { from: 2, to: 1, spin: Spin::Down, external: false },
            ],
            sign: -1,
            symmetry_factor: 1.0,
        };
        assert_eq!(d.order, 1);
        assert_eq!(d.vertices.len(), 2);
        assert_eq!(d.propagators.len(), 2);
    }

    #[test]
    fn test_count_loops() {
        // 2 vertices, 2 internal propagators → 1 loop
        let d = FeynmanDiagram {
            order: 1,
            vertices: vec![
                Vertex { id: 1, site: "i".into(), time: "τ1".into() },
                Vertex { id: 2, site: "j".into(), time: "τ2".into() },
            ],
            propagators: vec![
                Propagator { from: 1, to: 2, spin: Spin::Up, external: false },
                Propagator { from: 2, to: 1, spin: Spin::Down, external: false },
            ],
            sign: 1,
            symmetry_factor: 1.0,
        };
        assert_eq!(d.count_loops(), 1);
    }

    #[test]
    fn test_count_loops_with_external() {
        // 1 vertex, 2 external + 1 internal (self-loop) → 1 loop
        let d = FeynmanDiagram {
            order: 1,
            vertices: vec![
                Vertex { id: 1, site: "i".into(), time: "τ1".into() },
            ],
            propagators: vec![
                Propagator { from: 0, to: 1, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 2, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 1, spin: Spin::Down, external: false },
            ],
            sign: 1,
            symmetry_factor: 1.0,
        };
        assert_eq!(d.count_loops(), 1);
    }
}
