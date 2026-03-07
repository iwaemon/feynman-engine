use std::collections::HashMap;
use crate::algebra::operators::*;
use crate::algebra::wick::wick_theorem;
use crate::diagrams::graph::*;
use crate::models::hubbard::HubbardModel;

#[derive(Debug, Clone, Copy)]
pub enum Observable {
    Vacuum,
    SelfEnergy,
}

/// Generate all Feynman diagrams at given perturbation order using Wick's theorem.
pub fn generate_diagrams(_model: &HubbardModel, order: usize, observable: Observable) -> Vec<FeynmanDiagram> {
    let mut ops = Vec::new();

    // External legs for self-energy
    if let Observable::SelfEnergy = observable {
        ops.push(FermionOperator::creator("ext_in", Spin::Up, "τ_in"));
        ops.push(FermionOperator::annihilator("ext_out", Spin::Up, "τ_out"));
    }

    // Interaction vertices: each contributes c†↑ c↑ c†↓ c↓
    for v in 1..=order {
        let site = format!("v{}", v);
        let time = format!("τ{}", v);
        ops.push(FermionOperator::creator(&site, Spin::Up, &time));
        ops.push(FermionOperator::annihilator(&site, Spin::Up, &time));
        ops.push(FermionOperator::creator(&site, Spin::Down, &time));
        ops.push(FermionOperator::annihilator(&site, Spin::Down, &time));
    }

    let terms = wick_theorem(&ops);
    let factorial: f64 = (1..=order as u64).product::<u64>() as f64;

    terms.iter().filter_map(|term| {
        wick_term_to_diagram(term, order, &ops, 1.0 / factorial)
    }).collect()
}

fn wick_term_to_diagram(
    term: &WickTerm,
    order: usize,
    _ops: &[FermionOperator],
    sym_factor: f64,
) -> Option<FeynmanDiagram> {
    let mut vertex_map: HashMap<(String, String), usize> = HashMap::new();
    let mut vertices = Vec::new();
    let mut vid = 1usize;

    for c in &term.contractions {
        for op in [&c.creator, &c.annihilator] {
            if !op.site.starts_with("ext") {
                let key = (op.site.clone(), op.time.clone());
                if !vertex_map.contains_key(&key) {
                    vertex_map.insert(key.clone(), vid);
                    vertices.push(Vertex { id: vid, site: op.site.clone(), time: op.time.clone() });
                    vid += 1;
                }
            }
        }
    }

    let ext_in_id = 0usize;
    let ext_out_id = vid;

    let mut propagators = Vec::new();
    for c in &term.contractions {
        let from_id = if c.creator.site.starts_with("ext") {
            ext_in_id
        } else {
            vertex_map[&(c.creator.site.clone(), c.creator.time.clone())]
        };

        let to_id = if c.annihilator.site.starts_with("ext") {
            ext_out_id
        } else {
            vertex_map[&(c.annihilator.site.clone(), c.annihilator.time.clone())]
        };

        let is_ext = c.creator.site.starts_with("ext") || c.annihilator.site.starts_with("ext");
        propagators.push(Propagator {
            from: from_id,
            to: to_id,
            spin: c.creator.spin,
            external: is_ext,
        });
    }

    Some(FeynmanDiagram {
        order,
        vertices,
        propagators,
        sign: term.sign,
        symmetry_factor: sym_factor,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;

    #[test]
    fn test_generate_1st_order_self_energy() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let diagrams = generate_diagrams(&model, 1, Observable::SelfEnergy);
        assert!(!diagrams.is_empty());
        assert!(diagrams.iter().all(|d| d.order == 1));
    }

    #[test]
    fn test_generate_2nd_order_self_energy() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let diagrams = generate_diagrams(&model, 2, Observable::SelfEnergy);
        assert!(!diagrams.is_empty());
        assert!(diagrams.iter().all(|d| d.order == 2));
    }

    #[test]
    fn test_generate_vacuum() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let diagrams = generate_diagrams(&model, 1, Observable::Vacuum);
        assert!(!diagrams.is_empty());
    }
}
