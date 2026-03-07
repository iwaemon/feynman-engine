use crate::algebra::operators::*;
use itertools::Itertools;

/// Apply Wick's theorem to a time-ordered product of fermion operators.
pub fn wick_theorem(ops: &[FermionOperator]) -> Vec<WickTerm> {
    let n = ops.len();
    if n % 2 != 0 {
        return vec![];
    }

    let creators: Vec<(usize, &FermionOperator)> =
        ops.iter().enumerate().filter(|(_, op)| op.creation).collect();
    let annihilators: Vec<(usize, &FermionOperator)> =
        ops.iter().enumerate().filter(|(_, op)| !op.creation).collect();

    if creators.len() != annihilators.len() {
        return vec![];
    }

    // Group by spin
    let spins = [Spin::Up, Spin::Down];
    let mut spin_creators: Vec<Vec<(usize, &FermionOperator)>> = Vec::new();
    let mut spin_annihilators: Vec<Vec<(usize, &FermionOperator)>> = Vec::new();

    for &spin in &spins {
        let cs: Vec<_> = creators.iter().filter(|(_, op)| op.spin == spin).cloned().collect();
        let as_: Vec<_> = annihilators.iter().filter(|(_, op)| op.spin == spin).cloned().collect();
        if cs.len() != as_.len() {
            return vec![];
        }
        spin_creators.push(cs);
        spin_annihilators.push(as_);
    }

    // Generate all matchings per spin
    let mut spin_matchings: Vec<Vec<Vec<(usize, usize)>>> = Vec::new();
    for (cs, as_) in spin_creators.iter().zip(spin_annihilators.iter()) {
        if cs.is_empty() {
            spin_matchings.push(vec![vec![]]);
            continue;
        }
        let perms: Vec<Vec<(usize, usize)>> = (0..as_.len())
            .permutations(as_.len())
            .map(|perm| {
                perm.iter()
                    .enumerate()
                    .map(|(k, &p)| (cs[k].0, as_[p].0))
                    .collect()
            })
            .collect();
        spin_matchings.push(perms);
    }

    // Cartesian product
    let mut terms = Vec::new();
    let all_combos = cartesian_product(&spin_matchings);

    for combo in all_combos {
        let pairs: Vec<(usize, usize)> = combo.into_iter().flatten().collect();
        let contractions: Vec<Contraction> = pairs
            .iter()
            .map(|&(ci, ai)| Contraction::new(ops[ci].clone(), ops[ai].clone()))
            .collect();
        let sign = permutation_sign(&pairs);
        terms.push(WickTerm { sign, contractions });
    }

    terms
}

fn cartesian_product(lists: &[Vec<Vec<(usize, usize)>>]) -> Vec<Vec<Vec<(usize, usize)>>> {
    let mut result: Vec<Vec<Vec<(usize, usize)>>> = vec![vec![]];
    for list in lists {
        let mut new_result = Vec::new();
        for prev in &result {
            for item in list {
                let mut combined = prev.clone();
                combined.push(item.clone());
                new_result.push(combined);
            }
        }
        result = new_result;
    }
    result
}

/// Compute fermion sign from permutation that brings ops into contracted pairs.
fn permutation_sign(pairs: &[(usize, usize)]) -> i32 {
    let target: Vec<usize> = pairs.iter().flat_map(|&(c, a)| vec![c, a]).collect();
    let inversions: usize = target
        .iter()
        .enumerate()
        .map(|(i, &a)| target[i + 1..].iter().filter(|&&b| a > b).count())
        .sum();
    if inversions % 2 == 0 { 1 } else { -1 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wick_2_operators() {
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("j", Spin::Up, "τ2"),
        ];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].contractions[0].creator.site, "i");
        assert_eq!(terms[0].contractions[0].annihilator.site, "j");
    }

    #[test]
    fn test_wick_4_same_spin() {
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("j", Spin::Up, "τ2"),
            FermionOperator::creator("k", Spin::Up, "τ3"),
            FermionOperator::annihilator("l", Spin::Up, "τ4"),
        ];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 2);
        let mut signs: Vec<i32> = terms.iter().map(|t| t.sign).collect();
        signs.sort();
        assert_eq!(signs, vec![-1, 1]);
    }

    #[test]
    fn test_wick_4_different_spins() {
        // Hubbard vertex: spin conservation → only 1 term
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("i", Spin::Up, "τ1"),
            FermionOperator::creator("i", Spin::Down, "τ1"),
            FermionOperator::annihilator("i", Spin::Down, "τ1"),
        ];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0].contractions.len(), 2);
    }

    #[test]
    fn test_wick_8_operators() {
        // Two Hubbard vertices → 2!×2! = 4 terms
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("i", Spin::Up, "τ1"),
            FermionOperator::creator("i", Spin::Down, "τ1"),
            FermionOperator::annihilator("i", Spin::Down, "τ1"),
            FermionOperator::creator("j", Spin::Up, "τ2"),
            FermionOperator::annihilator("j", Spin::Up, "τ2"),
            FermionOperator::creator("j", Spin::Down, "τ2"),
            FermionOperator::annihilator("j", Spin::Down, "τ2"),
        ];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 4);
        assert!(terms.iter().all(|t| t.contractions.len() == 4));
    }

    #[test]
    fn test_wick_odd() {
        let ops = vec![FermionOperator::creator("i", Spin::Up, "τ1")];
        assert_eq!(wick_theorem(&ops).len(), 0);
    }

    #[test]
    fn test_wick_spin_mismatch() {
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("j", Spin::Down, "τ2"),
        ];
        assert_eq!(wick_theorem(&ops).len(), 0);
    }
}
