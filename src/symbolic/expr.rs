use std::fmt;

/// Symbolic expression tree for Feynman diagram expressions
#[derive(Debug, Clone)]
pub enum Expr {
    /// Free Green's function G₀(k, iωₙ)
    G0 { k: String, omega: String },
    /// Interaction vertex U
    U,
    /// Summation (1/Nβ) Σ_{var}
    Sum { var: String, body: Box<Expr> },
    /// Product of expressions
    Mul(Vec<Expr>),
    /// Sum of expressions
    Add(Vec<Expr>),
    /// Inverse 1/x
    Inv(Box<Expr>),
    /// Scalar value
    Scalar(f64),
    /// Negation
    Neg(Box<Expr>),
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Expr::G0 { k, omega } => write!(f, "G₀({}, i{})", k, omega),
            Expr::U => write!(f, "U"),
            Expr::Sum { var, body } => write!(f, "Σ_{} [{}]", var, body),
            Expr::Mul(terms) => {
                let parts: Vec<String> = terms.iter().map(|t| format!("{}", t)).collect();
                write!(f, "({})", parts.join(" × "))
            }
            Expr::Add(terms) => {
                let parts: Vec<String> = terms.iter().map(|t| format!("{}", t)).collect();
                write!(f, "({})", parts.join(" + "))
            }
            Expr::Inv(inner) => write!(f, "1/({})", inner),
            Expr::Scalar(v) => write!(f, "{}", v),
            Expr::Neg(inner) => write!(f, "-({})", inner),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_display() {
        let e = Expr::Mul(vec![
            Expr::Neg(Box::new(Expr::U)),
            Expr::G0 { k: "k".into(), omega: "ωₙ".into() },
        ]);
        let s = format!("{}", e);
        assert!(s.contains("U"));
        assert!(s.contains("G₀"));
    }

    #[test]
    fn test_sum_display() {
        let e = Expr::Sum {
            var: "k".into(),
            body: Box::new(Expr::G0 { k: "k".into(), omega: "ωₙ".into() }),
        };
        let s = format!("{}", e);
        assert!(s.contains("Σ_k"));
        assert!(s.contains("G₀"));
    }
}
