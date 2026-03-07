use crate::diagrams::graph::FeynmanDiagram;

/// Serializable wrapper for classified diagrams
#[derive(serde::Serialize)]
pub struct ClassifiedDiagramSet {
    pub diagrams: Vec<ClassifiedDiagram>,
}

#[derive(serde::Serialize)]
pub struct ClassifiedDiagram {
    pub topology_id: usize,
    pub weight: i32,
    pub diagram: FeynmanDiagram,
}

/// Convert classified diagrams to JSON string
pub fn to_json_all(classified: &[(FeynmanDiagram, i32)]) -> String {
    let set = ClassifiedDiagramSet {
        diagrams: classified
            .iter()
            .enumerate()
            .map(|(i, (diag, weight))| ClassifiedDiagram {
                topology_id: i + 1,
                weight: *weight,
                diagram: diag.clone(),
            })
            .collect(),
    };
    serde_json::to_string_pretty(&set).expect("Failed to serialize diagrams to JSON")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::operators::Spin;
    use crate::diagrams::graph::{Vertex, Propagator};

    #[test]
    fn test_to_json_all_basic() {
        let d = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex {
                id: 1,
                site: "v1".into(),
                time: "τ1".into(),
            }],
            propagators: vec![
                Propagator { from: 0, to: 1, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 2, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 1, spin: Spin::Down, external: false },
            ],
            sign: -1,
            symmetry_factor: 1.0,
        };

        let classified = vec![(d, -1)];
        let json = to_json_all(&classified);

        assert!(json.contains("\"topology_id\": 1"));
        assert!(json.contains("\"weight\": -1"));
        assert!(json.contains("\"spin\": \"Up\""));
        assert!(json.contains("\"spin\": \"Down\""));
        assert!(json.contains("\"external\": true"));
        assert!(json.contains("\"external\": false"));
    }

    #[test]
    fn test_to_json_all_multiple() {
        let d1 = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex { id: 1, site: "v1".into(), time: "τ1".into() }],
            propagators: vec![
                Propagator { from: 1, to: 1, spin: Spin::Up, external: false },
            ],
            sign: 1,
            symmetry_factor: 1.0,
        };
        let d2 = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex { id: 1, site: "v1".into(), time: "τ1".into() }],
            propagators: vec![
                Propagator { from: 0, to: 1, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 2, spin: Spin::Up, external: true },
            ],
            sign: -1,
            symmetry_factor: 1.0,
        };

        let classified = vec![(d1, 1), (d2, -1)];
        let json = to_json_all(&classified);
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed["diagrams"].as_array().unwrap().len(), 2);
        assert_eq!(parsed["diagrams"][0]["topology_id"], 1);
        assert_eq!(parsed["diagrams"][1]["topology_id"], 2);
    }
}
