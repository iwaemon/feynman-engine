use crate::diagrams::graph::FeynmanDiagram;
use crate::models::hubbard::HubbardModel;
use crate::symbolic::expr::Expr;

/// Result of applying Feynman rules to a diagram
pub struct FeynmanExpression {
    pub n_loops: usize,
    pub prefactor: f64,
    pub expr: Expr,
    pub description: String,
}

/// Convert a Feynman diagram to its mathematical expression via Feynman rules.
pub fn apply_feynman_rules(d: &FeynmanDiagram, model: &HubbardModel) -> FeynmanExpression {
    let n_loops = d.count_loops();

    // Prefactor: (-U)^order × sign × symmetry_factor
    let prefactor = (-model.u).powi(d.order as i32) * (d.sign as f64) * d.symmetry_factor;

    // Build expression tree
    let internal_props: Vec<_> = d.propagators.iter()
        .filter(|p| !p.external)
        .enumerate()
        .collect();

    let mut factors: Vec<Expr> = Vec::new();
    factors.push(Expr::Scalar(prefactor));

    for (i, _prop) in &internal_props {
        let k_var = format!("k{}", i + 1);
        let omega_var = format!("ω{}", i + 1);
        factors.push(Expr::G0 { k: k_var, omega: omega_var });
    }

    let integrand = if factors.len() == 1 {
        factors.pop().unwrap()
    } else {
        Expr::Mul(factors)
    };

    // Wrap in loop summations
    let mut expr = integrand;
    for l in (1..=n_loops).rev() {
        expr = Expr::Sum {
            var: format!("k{},ω{}", l, l),
            body: Box::new(expr),
        };
    }

    let n_internal = internal_props.len();
    let desc = format!("Order {}: {} loop(s), {} internal propagators", d.order, n_loops, n_internal);

    FeynmanExpression {
        n_loops,
        prefactor,
        expr,
        description: desc,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagrams::generate::*;
    use crate::diagrams::classify::classify_diagrams;
    use crate::models::lattice::SquareLattice;

    #[test]
    fn test_feynman_rules_1st_order() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let diagrams = generate_diagrams(&model, 1, Observable::SelfEnergy);
        let classified = classify_diagrams(&diagrams);

        for (diag, _weight) in &classified {
            let expr = apply_feynman_rules(diag, &model);
            assert_eq!(expr.n_loops, diag.count_loops());
            assert_ne!(expr.prefactor, 0.0);
            let s = format!("{}", expr.expr);
            assert!(s.contains("G₀"));
        }
    }

    #[test]
    fn test_feynman_expression_description() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let diagrams = generate_diagrams(&model, 2, Observable::SelfEnergy);
        let classified = classify_diagrams(&diagrams);

        for (diag, _) in &classified {
            let expr = apply_feynman_rules(diag, &model);
            assert!(expr.description.contains("Order 2"));
            assert!(expr.description.contains("loop"));
        }
    }
}
