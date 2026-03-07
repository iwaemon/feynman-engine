use std::collections::BTreeMap;
use crate::algebra::operators::Spin;
use crate::diagrams::graph::FeynmanDiagram;

/// Classify diagrams by topology. Returns unique representatives with total weight.
pub fn classify_diagrams(diagrams: &[FeynmanDiagram]) -> Vec<(FeynmanDiagram, i32)> {
    let mut groups: Vec<(FeynmanDiagram, i32, Vec<(usize, usize, usize, usize, usize)>)> = Vec::new();

    for d in diagrams {
        let sig = diagram_signature(d);
        let mut found = false;
        for (_, weight, existing_sig) in groups.iter_mut() {
            if *existing_sig == sig {
                *weight += d.sign;
                found = true;
                break;
            }
        }
        if !found {
            groups.push((d.clone(), d.sign, sig));
        }
    }

    groups.into_iter()
        .filter(|(_, w, _)| *w != 0)
        .map(|(d, w, _)| (d, w))
        .collect()
}

/// Compute a canonical signature for diagram topology.
/// Per vertex: (in_up, out_up, in_down, out_down, self_loops), sorted.
fn diagram_signature(d: &FeynmanDiagram) -> Vec<(usize, usize, usize, usize, usize)> {
    let vert_ids: BTreeMap<usize, usize> = d.vertices.iter()
        .enumerate()
        .map(|(i, v)| (v.id, i))
        .collect();

    let n = d.vertices.len();
    if n == 0 {
        return vec![];
    }

    // Per vertex: (in_up, out_up, in_down, out_down, self_loops)
    let mut descriptors = vec![(0usize, 0usize, 0usize, 0usize, 0usize); n];

    for p in &d.propagators {
        if p.external {
            continue;
        }

        let from_idx = match vert_ids.get(&p.from) {
            Some(&i) => i,
            None => continue,
        };
        let to_idx = match vert_ids.get(&p.to) {
            Some(&i) => i,
            None => continue,
        };

        // Self-loop
        if from_idx == to_idx {
            descriptors[from_idx].4 += 1;
        }

        match p.spin {
            Spin::Up => {
                descriptors[from_idx].1 += 1; // out_up
                descriptors[to_idx].0 += 1;   // in_up
            }
            Spin::Down => {
                descriptors[from_idx].3 += 1; // out_down
                descriptors[to_idx].2 += 1;   // in_down
            }
        }
    }

    descriptors.sort();
    descriptors
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagrams::generate::*;
    use crate::models::lattice::SquareLattice;
    use crate::models::hubbard::HubbardModel;

    #[test]
    fn test_classify_1st_order() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let diagrams = generate_diagrams(&model, 1, Observable::SelfEnergy);
        let classified = classify_diagrams(&diagrams);

        // 1st order self-energy: Hartree (tadpole) + Fock (exchange) = 2 topologies
        assert_eq!(classified.len(), 2);

        // Each weight should be nonzero
        for (_, weight) in &classified {
            assert_ne!(*weight, 0);
        }
    }

    #[test]
    fn test_classify_2nd_order() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let diagrams = generate_diagrams(&model, 2, Observable::SelfEnergy);
        let classified = classify_diagrams(&diagrams);

        // Should have fewer unique topologies than total Wick terms
        assert!(classified.len() <= diagrams.len());
        assert!(!classified.is_empty());

        // Each weight should be nonzero
        for (_, weight) in &classified {
            assert_ne!(*weight, 0);
        }
    }
}
