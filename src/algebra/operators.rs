#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FermionOperator {
    pub creation: bool,
    pub site: String,
    pub spin: Spin,
    pub time: String,
}

impl FermionOperator {
    pub fn creator(site: &str, spin: Spin, time: &str) -> Self {
        Self {
            creation: true,
            site: site.to_string(),
            spin,
            time: time.to_string(),
        }
    }

    pub fn annihilator(site: &str, spin: Spin, time: &str) -> Self {
        Self {
            creation: false,
            site: site.to_string(),
            spin,
            time: time.to_string(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Contraction {
    pub creator: FermionOperator,
    pub annihilator: FermionOperator,
}

impl Contraction {
    pub fn new(creator: FermionOperator, annihilator: FermionOperator) -> Self {
        Self {
            creator,
            annihilator,
        }
    }
}

#[derive(Debug, Clone)]
pub struct WickTerm {
    pub sign: i32,
    pub contractions: Vec<Contraction>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_operators() {
        let c_dag = FermionOperator::creator("i", Spin::Up, "τ1");
        let c = FermionOperator::annihilator("i", Spin::Up, "τ1");
        assert!(c_dag.creation);
        assert!(!c.creation);
        assert_eq!(c_dag.site, "i");
        assert_eq!(c_dag.spin, Spin::Up);
    }

    #[test]
    fn test_contraction() {
        let c_dag = FermionOperator::creator("i", Spin::Up, "τ1");
        let c = FermionOperator::annihilator("j", Spin::Up, "τ2");
        let cont = Contraction::new(c_dag.clone(), c.clone());
        assert_eq!(cont.creator.site, "i");
        assert_eq!(cont.annihilator.site, "j");
        assert_eq!(cont.creator.spin, cont.annihilator.spin);
    }

    #[test]
    fn test_wick_term() {
        let c1 = FermionOperator::creator("i", Spin::Up, "τ1");
        let a1 = FermionOperator::annihilator("j", Spin::Up, "τ2");
        let c2 = FermionOperator::creator("i", Spin::Down, "τ1");
        let a2 = FermionOperator::annihilator("j", Spin::Down, "τ2");
        let term = WickTerm {
            sign: 1,
            contractions: vec![Contraction::new(c1, a1), Contraction::new(c2, a2)],
        };
        assert_eq!(term.sign, 1);
        assert_eq!(term.contractions.len(), 2);
    }
}
