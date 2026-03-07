# FeynmanEngine (Rust) Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a Rust perturbation engine that auto-generates Feynman diagrams via Wick's theorem and resums ladder diagrams for the attractive Hubbard model to find Tc.

**Architecture:** Bottom-up from operator algebra (Wick's theorem) through diagram generation/classification, to Bethe-Salpeter resummation. Expression tree (AST) for symbolic expressions. Matsubara formalism for finite temperature. rayon for parallel numerical evaluation.

**Tech Stack:** Rust, num-complex, rayon, rustfft, plotters, itertools

---

## Task 1: Project Scaffolding

**Files:**
- Create: `Cargo.toml`
- Create: `src/lib.rs`
- Create: `src/main.rs`

**Step 1: Initialize Rust project**

```bash
cd /Users/shumpei/work/feynman
cargo init --name feynman-engine
```

**Step 2: Set up Cargo.toml with dependencies**

Write `Cargo.toml`:

```toml
[package]
name = "feynman-engine"
version = "0.1.0"
edition = "2021"

[dependencies]
num-complex = "0.4"
rayon = "1.10"
itertools = "0.13"

[dev-dependencies]
approx = "0.5"
```

Note: rustfft and plotters are added in later tasks when needed.

**Step 3: Create module structure**

Write `src/lib.rs`:

```rust
pub mod algebra;
pub mod models;
pub mod diagrams;
pub mod symbolic;
pub mod resummation;
pub mod numerical;
pub mod visualization;

pub mod prelude {
    pub use crate::algebra::operators::*;
    pub use crate::algebra::wick::*;
    pub use crate::models::lattice::*;
    pub use crate::models::hubbard::*;
    pub use crate::diagrams::graph::*;
    pub use crate::diagrams::generate::*;
    pub use crate::diagrams::classify::*;
    pub use crate::diagrams::rules::*;
    pub use crate::symbolic::expr::*;
    pub use crate::resummation::channels::*;
    pub use crate::resummation::bethe_salpeter::*;
    pub use crate::numerical::greens_function::*;
    pub use crate::numerical::matsubara::*;
    pub use crate::visualization::dot::*;
}
```

Create all module directories and mod.rs files:

```
src/algebra/mod.rs          → pub mod operators; pub mod wick;
src/models/mod.rs           → pub mod lattice; pub mod hubbard;
src/diagrams/mod.rs         → pub mod graph; pub mod generate; pub mod classify; pub mod rules;
src/symbolic/mod.rs         → pub mod expr; pub mod simplify; pub mod display;
src/resummation/mod.rs      → pub mod channels; pub mod bethe_salpeter;
src/numerical/mod.rs        → pub mod greens_function; pub mod matsubara;
src/visualization/mod.rs    → pub mod dot; pub mod plot;
```

Each leaf module file starts empty or with a placeholder comment.

Write `src/main.rs`:

```rust
use feynman_engine::prelude::*;

fn main() {
    println!("FeynmanEngine v0.1.0");
}
```

**Step 4: Verify it compiles**

```bash
cargo build
cargo test
```

Expected: Compiles and tests pass (no tests yet, just compilation check).

**Step 5: Commit**

```bash
git add Cargo.toml Cargo.lock src/
git commit -m "feat: scaffold feynman-engine Rust project with module structure"
```

---

## Task 2: Fermion Operator Types

**Files:**
- Create: `src/algebra/operators.rs`
- Modify: `src/algebra/mod.rs`

**Step 1: Write failing tests**

At bottom of `src/algebra/operators.rs`:

```rust
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
            contractions: vec![
                Contraction::new(c1, a1),
                Contraction::new(c2, a2),
            ],
        };
        assert_eq!(term.sign, 1);
        assert_eq!(term.contractions.len(), 2);
    }
}
```

**Step 2: Run test to verify it fails**

```bash
cargo test
```

Expected: FAIL — types not defined.

**Step 3: Write implementation (above tests in same file)**

```rust
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
        Self { creation: true, site: site.to_string(), spin, time: time.to_string() }
    }

    pub fn annihilator(site: &str, spin: Spin, time: &str) -> Self {
        Self { creation: false, site: site.to_string(), spin, time: time.to_string() }
    }
}

#[derive(Debug, Clone)]
pub struct Contraction {
    pub creator: FermionOperator,
    pub annihilator: FermionOperator,
}

impl Contraction {
    pub fn new(creator: FermionOperator, annihilator: FermionOperator) -> Self {
        Self { creator, annihilator }
    }
}

#[derive(Debug, Clone)]
pub struct WickTerm {
    pub sign: i32,
    pub contractions: Vec<Contraction>,
}
```

**Step 4: Run tests**

```bash
cargo test algebra::operators
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/algebra/
git commit -m "feat: add FermionOperator, Contraction, WickTerm types"
```

---

## Task 3: Wick's Theorem

**Files:**
- Create: `src/algebra/wick.rs`
- Modify: `src/algebra/mod.rs`

**Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::operators::*;

    #[test]
    fn test_wick_2_operators() {
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("j", Spin::Up, "τ2"),
        ];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 1);
        // Sign: bringing c† and c together from T-ordered product
        assert_eq!(terms[0].contractions[0].creator.site, "i");
        assert_eq!(terms[0].contractions[0].annihilator.site, "j");
    }

    #[test]
    fn test_wick_4_operators_same_spin() {
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("j", Spin::Up, "τ2"),
            FermionOperator::creator("k", Spin::Up, "τ3"),
            FermionOperator::annihilator("l", Spin::Up, "τ4"),
        ];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 2);
        let signs: Vec<i32> = {
            let mut s: Vec<i32> = terms.iter().map(|t| t.sign).collect();
            s.sort();
            s
        };
        assert_eq!(signs, vec![-1, 1]);
    }

    #[test]
    fn test_wick_4_operators_different_spins() {
        // Hubbard vertex: c†↑ c↑ c†↓ c↓ — spin conservation → 1 term
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
    fn test_wick_8_operators_2nd_order() {
        // Two Hubbard vertices: 8 operators
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
        // 2!×2! = 4 terms
        assert_eq!(terms.len(), 4);
        assert!(terms.iter().all(|t| t.contractions.len() == 4));
    }

    #[test]
    fn test_wick_odd_operators() {
        let ops = vec![FermionOperator::creator("i", Spin::Up, "τ1")];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 0);
    }

    #[test]
    fn test_wick_spin_mismatch() {
        let ops = vec![
            FermionOperator::creator("i", Spin::Up, "τ1"),
            FermionOperator::annihilator("j", Spin::Down, "τ2"),
        ];
        let terms = wick_theorem(&ops);
        assert_eq!(terms.len(), 0);
    }
}
```

**Step 2: Run test to verify it fails**

```bash
cargo test algebra::wick
```

Expected: FAIL.

**Step 3: Write implementation**

```rust
use crate::algebra::operators::*;
use itertools::Itertools;

/// Apply Wick's theorem to a time-ordered product of fermion operators.
/// Returns all fully-contracted terms with correct fermion signs.
/// Enforces spin conservation (G₀ is spin-diagonal).
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

    // Generate all matchings per spin, then take Cartesian product
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

    // Cartesian product of spin matchings
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

/// Compute fermion sign from the permutation that brings operators into contracted pairs.
fn permutation_sign(pairs: &[(usize, usize)]) -> i32 {
    let target: Vec<usize> = pairs.iter().flat_map(|&(c, a)| vec![c, a]).collect();
    let inversions = target.iter().enumerate().map(|(i, &a)| {
        target[i + 1..].iter().filter(|&&b| a > b).count()
    }).sum::<usize>();
    if inversions % 2 == 0 { 1 } else { -1 }
}
```

**Step 4: Run tests**

```bash
cargo test algebra::wick
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/algebra/
git commit -m "feat: implement Wick's theorem with spin conservation"
```

---

## Task 4: Lattice and Hubbard Model

**Files:**
- Create: `src/models/lattice.rs`
- Create: `src/models/hubbard.rs`

**Step 1: Write failing tests**

In `src/models/lattice.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_square_lattice() {
        let lat = SquareLattice::new(4, 4);
        assert_eq!(lat.num_sites(), 16);
        assert_eq!(lat.k_grid().len(), 16);
    }

    #[test]
    fn test_dispersion() {
        let lat = SquareLattice::new(4, 4);
        // Γ point: ε(0,0) = -4t
        assert!((lat.dispersion(&[0.0, 0.0], 1.0) - (-4.0)).abs() < 1e-12);
        // M point: ε(π,π) = 4t
        assert!((lat.dispersion(&[PI, PI], 1.0) - 4.0).abs() < 1e-12);
    }
}
```

In `src/models/hubbard.rs`:

```rust
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
```

**Step 2: Run test to verify it fails**

```bash
cargo test models
```

**Step 3: Write implementation**

`src/models/lattice.rs`:

```rust
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
```

`src/models/hubbard.rs`:

```rust
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
```

**Step 4: Run tests**

```bash
cargo test models
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/models/
git commit -m "feat: add SquareLattice, HubbardModel, ThermalParams"
```

---

## Task 5: Free Green's Function

**Files:**
- Create: `src/numerical/greens_function.rs`

**Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;
    use crate::models::hubbard::{HubbardModel, ThermalParams};
    use std::f64::consts::PI;

    #[test]
    fn test_evaluate_g0() {
        let lat = SquareLattice::new(4, 4);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let params = ThermalParams::new(10.0, 0.0, 128);

        let k = [0.0, 0.0]; // ε = -4
        let n = 0; // ω₀ = π/β
        let g = evaluate_g0(&model, &params, &k, n);

        let omega = PI / params.beta;
        let eps = model.dispersion(&k);
        let expected = 1.0 / Complex64::new(-eps + params.mu, omega);
        assert!((g - expected).norm() < 1e-12);
    }

    #[test]
    fn test_g0_sum_rule() {
        // (1/β) Σₙ G₀(k, iωₙ) → n_F(εₖ-μ)
        // For εₖ=0, μ=0: n_F(0) = 0.5
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, 0.0);
        let params = ThermalParams::new(20.0, 0.0, 2048);

        let k = [PI / 2.0, PI / 2.0]; // ε ≈ 0
        let eps = model.dispersion(&k);
        assert!(eps.abs() < 1e-10);

        let mut nk = 0.0;
        let n_max = params.n_matsubara as i32;
        for n in -n_max..n_max {
            let g = evaluate_g0(&model, &params, &k, n);
            nk += g.re / params.beta;
        }
        assert!((nk - 0.5).abs() < 0.01);
    }
}
```

**Step 2: Run test, verify fails. Step 3: Implement.**

```rust
use num_complex::Complex64;
use std::f64::consts::PI;
use crate::models::hubbard::{HubbardModel, ThermalParams};

/// Fermionic Matsubara frequency ωₙ = (2n+1)π/β
pub fn matsubara_freq(n: i32, beta: f64) -> f64 {
    (2 * n + 1) as f64 * PI / beta
}

/// Free Matsubara Green's function G₀(k, iωₙ) = 1/(iωₙ - εₖ + μ)
pub fn evaluate_g0(model: &HubbardModel, params: &ThermalParams,
                   k: &[f64; 2], n: i32) -> Complex64 {
    let omega = matsubara_freq(n, params.beta);
    let eps = model.dispersion(k);
    1.0 / Complex64::new(-eps + params.mu, omega)
}
```

**Step 4: Run tests. Step 5: Commit.**

```bash
cargo test numerical::greens_function
git add src/numerical/
git commit -m "feat: add free Matsubara Green's function G₀(k, iωₙ)"
```

---

## Task 6: PP Susceptibility

**Files:**
- Create: `src/numerical/matsubara.rs`

**Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;
    use crate::models::hubbard::{HubbardModel, ThermalParams};

    #[test]
    fn test_pp_susceptibility_real_at_q0() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let params = ThermalParams::new(5.0, 0.0, 512);

        let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);
        assert!(chi0.im.abs() < 1e-6); // should be real at q=0, ν=0
        assert!(chi0.re < 0.0);         // negative for pp bubble
    }

    #[test]
    fn test_pp_susceptibility_increases_with_beta() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);

        let params_high_t = ThermalParams::new(2.0, 0.0, 256);
        let params_low_t = ThermalParams::new(10.0, 0.0, 1024);

        let chi_high = compute_pp_susceptibility(&model, &params_high_t, &[0.0, 0.0], 0).re.abs();
        let chi_low = compute_pp_susceptibility(&model, &params_low_t, &[0.0, 0.0], 0).re.abs();

        assert!(chi_low > chi_high);
    }
}
```

**Step 2: Run test, verify fails. Step 3: Implement.**

```rust
use num_complex::Complex64;
use std::f64::consts::PI;
use crate::models::hubbard::{HubbardModel, ThermalParams};
use crate::numerical::greens_function::evaluate_g0;

/// Bosonic Matsubara frequency νₘ = 2mπ/β
pub fn bosonic_matsubara_freq(m: i32, beta: f64) -> f64 {
    2.0 * m as f64 * PI / beta
}

/// Particle-particle bubble:
/// χ₀(q, iνₘ) = -(1/Nβ) Σₖ Σₙ G₀(k, iωₙ) G₀(q-k, iνₘ - iωₙ)
pub fn compute_pp_susceptibility(
    model: &HubbardModel,
    params: &ThermalParams,
    q: &[f64; 2],
    m: i32,
) -> Complex64 {
    let kpoints = model.lattice.k_grid();
    let n_sites = kpoints.len() as f64;
    let n_max = params.n_matsubara as i32;

    let mut chi0 = Complex64::new(0.0, 0.0);

    for k in &kpoints {
        let qmk = [q[0] - k[0], q[1] - k[1]];
        for n in -n_max..n_max {
            let g1 = evaluate_g0(model, params, k, n);
            // iνₘ - iωₙ: bosonic freq - fermionic freq = fermionic freq
            // index for second G₀: m - 1 - n
            let g2 = evaluate_g0(model, params, &qmk, m - 1 - n);
            chi0 += g1 * g2;
        }
    }

    -chi0 / (n_sites * params.beta)
}
```

Note: The frequency index `m - 1 - n` comes from: iνₘ - iωₙ = 2mπ/β - (2n+1)π/β = (2(m-1-n)+1)π/β = iω_{m-1-n}. The physics tests (real at q=0, correct sign, temperature dependence) will validate this.

**Step 4: Run tests. Step 5: Commit.**

```bash
cargo test numerical::matsubara
git add src/numerical/
git commit -m "feat: add particle-particle susceptibility χ₀(q, iνₘ)"
```

---

## Task 7: Bethe-Salpeter Equation and Tc Finder

**Files:**
- Create: `src/resummation/channels.rs`
- Create: `src/resummation/bethe_salpeter.rs`

**Step 1: Write failing tests**

In `src/resummation/bethe_salpeter.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;
    use crate::models::hubbard::{HubbardModel, ThermalParams};
    use crate::numerical::matsubara::compute_pp_susceptibility;

    #[test]
    fn test_tmatrix() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);
        let params = ThermalParams::new(5.0, 0.0, 512);

        let t_val = solve_tmatrix(&model, &params, &[0.0, 0.0], 0);
        let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);
        let expected = model.u / (1.0 - model.u * chi0);

        assert!((t_val - expected).norm() < 1e-10);
    }

    #[test]
    fn test_find_tc() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, -2.0);

        let tc = find_tc(&model, (0.01, 2.0), 512, 1e-4);

        assert!(tc > 0.0);
        assert!(tc < 4.0);

        // Verify bracket
        let params_above = ThermalParams::new(1.0 / (tc * 1.05), 0.0, 512);
        let params_below = ThermalParams::new(1.0 / (tc * 0.95), 0.0, 512);
        let chi_above = compute_pp_susceptibility(&model, &params_above, &[0.0, 0.0], 0);
        let chi_below = compute_pp_susceptibility(&model, &params_below, &[0.0, 0.0], 0);

        assert!((1.0 - model.u * chi_above).re > 0.0);
        assert!((1.0 - model.u * chi_below).re < 0.0);
    }

    #[test]
    fn test_repulsive_no_pp_instability() {
        let lat = SquareLattice::new(8, 8);
        let model = HubbardModel::new(lat, 1.0, 2.0); // repulsive
        let params = ThermalParams::new(20.0, 0.0, 512);

        let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);
        assert!((1.0 - model.u * chi0).re > 0.0);
    }
}
```

**Step 2: Run test, verify fails. Step 3: Implement.**

`src/resummation/channels.rs`:

```rust
// Channel definitions for resummation
// Currently pp-ladder only; ph-bubble for future extension
```

`src/resummation/bethe_salpeter.rs`:

```rust
use num_complex::Complex64;
use crate::models::hubbard::{HubbardModel, ThermalParams};
use crate::numerical::matsubara::compute_pp_susceptibility;

/// T-matrix from pp-ladder resummation: T(q,iνₘ) = U / (1 - U·χ₀(q,iνₘ))
pub fn solve_tmatrix(
    model: &HubbardModel,
    params: &ThermalParams,
    q: &[f64; 2],
    m: i32,
) -> Complex64 {
    let chi0 = compute_pp_susceptibility(model, params, q, m);
    model.u / (1.0 - model.u * chi0)
}

/// Find superconducting Tc via Thouless criterion: 1 - U·χ₀(q=0, 0; T) = 0.
/// Uses bisection method.
pub fn find_tc(
    model: &HubbardModel,
    t_range: (f64, f64),
    n_matsubara: usize,
    tol: f64,
) -> f64 {
    let (mut t_lo, mut t_hi) = t_range;

    let criterion = |temp: f64| -> f64 {
        let params = ThermalParams::new(1.0 / temp, 0.0, n_matsubara);
        let chi0 = compute_pp_susceptibility(model, &params, &[0.0, 0.0], 0);
        (1.0 - model.u * chi0).re
    };

    let mut f_lo = criterion(t_lo);
    let f_hi = criterion(t_hi);

    assert!(
        f_lo * f_hi < 0.0,
        "Thouless criterion not bracketed: f({})={}, f({})={}",
        t_lo, f_lo, t_hi, f_hi
    );

    while (t_hi - t_lo) > tol {
        let t_mid = (t_lo + t_hi) / 2.0;
        let f_mid = criterion(t_mid);
        if f_mid * f_lo < 0.0 {
            t_hi = t_mid;
        } else {
            t_lo = t_mid;
            f_lo = f_mid;
        }
    }

    (t_lo + t_hi) / 2.0
}
```

**Step 4: Run tests. Step 5: Commit.**

```bash
cargo test resummation
git add src/resummation/
git commit -m "feat: add PP ladder resummation and Thouless Tc finder"
```

---

## Task 8: Diagram Graph Types

**Files:**
- Create: `src/diagrams/graph.rs`

**Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::operators::Spin;

    #[test]
    fn test_feynman_diagram() {
        let v1 = Vertex { id: 1, site: "i".into(), time: "τ1".into() };
        let v2 = Vertex { id: 2, site: "j".into(), time: "τ2".into() };

        let p1 = Propagator { from: 1, to: 2, spin: Spin::Up, external: false };
        let p2 = Propagator { from: 2, to: 1, spin: Spin::Down, external: false };

        let d = FeynmanDiagram {
            order: 1,
            vertices: vec![v1, v2],
            propagators: vec![p1, p2],
            sign: -1,
            symmetry_factor: 1.0,
        };

        assert_eq!(d.order, 1);
        assert_eq!(d.vertices.len(), 2);
        assert_eq!(d.propagators.len(), 2);
    }

    #[test]
    fn test_count_loops() {
        // 2 vertices, 2 internal props → 1 loop
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
}
```

**Step 3: Implement**

```rust
use crate::algebra::operators::Spin;

#[derive(Debug, Clone)]
pub struct Vertex {
    pub id: usize,
    pub site: String,
    pub time: String,
}

#[derive(Debug, Clone)]
pub struct Propagator {
    pub from: usize,  // vertex id
    pub to: usize,    // vertex id
    pub spin: Spin,
    pub external: bool,
}

#[derive(Debug, Clone)]
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
```

**Step 4: Run tests. Step 5: Commit.**

```bash
cargo test diagrams::graph
git add src/diagrams/
git commit -m "feat: add Vertex, Propagator, FeynmanDiagram types with loop counting"
```

---

## Task 9: Diagram Generation

**Files:**
- Create: `src/diagrams/generate.rs`

**Step 1: Write failing tests**

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::lattice::SquareLattice;
    use crate::models::hubbard::HubbardModel;

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
```

**Step 3: Implement**

```rust
use crate::algebra::operators::*;
use crate::algebra::wick::wick_theorem;
use crate::diagrams::graph::*;
use crate::models::hubbard::HubbardModel;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
pub enum Observable {
    Vacuum,
    SelfEnergy,
}

pub fn generate_diagrams(model: &HubbardModel, order: usize, observable: Observable) -> Vec<FeynmanDiagram> {
    let mut ops = Vec::new();
    let mut ext_indices = Vec::new();

    // External legs
    if let Observable::SelfEnergy = observable {
        ops.push(FermionOperator::creator("ext_in", Spin::Up, "τ_in"));
        ext_indices.push(ops.len() - 1);
        ops.push(FermionOperator::annihilator("ext_out", Spin::Up, "τ_out"));
        ext_indices.push(ops.len() - 1);
    }

    // Interaction vertices
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
        wick_term_to_diagram(term, order, &ext_indices, &ops, 1.0 / factorial)
    }).collect()
}

fn wick_term_to_diagram(
    term: &WickTerm,
    order: usize,
    ext_indices: &[usize],
    ops: &[FermionOperator],
    sym_factor: f64,
) -> Option<FeynmanDiagram> {
    let mut vertex_map: HashMap<(String, String), usize> = HashMap::new();
    let mut vertices = Vec::new();
    let mut vid = 1usize; // 0 reserved for ext

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

    let ext_in_id = 0usize;   // pseudo-vertex for external in
    let ext_out_id = vid;      // pseudo-vertex for external out

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
```

**Step 4: Run tests. Step 5: Commit.**

```bash
cargo test diagrams::generate
git add src/diagrams/
git commit -m "feat: generate Feynman diagrams from Wick contractions"
```

---

## Task 10: Diagram Classification

**Files:**
- Create: `src/diagrams/classify.rs`

**Step 1: Write tests, Step 3: Implement**

Uses sorted vertex-descriptor signatures for topology matching. Groups diagrams and sums signs.

```rust
pub fn classify_diagrams(diagrams: &[FeynmanDiagram]) -> Vec<(FeynmanDiagram, i32)> {
    // For each diagram, compute canonical signature
    // Group by signature, sum signs
    // Return unique representatives with total weight
}
```

Key function: `diagram_signature()` computes per-vertex (in_degree_up, out_degree_up, in_degree_down, out_degree_down, self_loops, ext_connections), then sorts.

**Commit:**

```bash
git commit -m "feat: classify diagrams by topology using graph signatures"
```

---

## Task 11: Feynman Rules (Diagram → Expression)

**Files:**
- Create: `src/symbolic/expr.rs`
- Create: `src/diagrams/rules.rs`

**Step 1: Tests, Step 3: Implement**

`src/symbolic/expr.rs` — the expression tree:

```rust
#[derive(Debug, Clone)]
pub enum Expr {
    G0 { k: String, omega: String },
    U,
    Sum { var: String, body: Box<Expr> },
    Mul(Vec<Expr>),
    Add(Vec<Expr>),
    Inv(Box<Expr>),
    Scalar(f64),
    Neg(Box<Expr>),
}

impl std::fmt::Display for Expr { /* human-readable output */ }
```

`src/diagrams/rules.rs`:

```rust
pub struct FeynmanExpression {
    pub n_loops: usize,
    pub prefactor: f64,
    pub expr: Expr,
    pub description: String,
}

pub fn apply_feynman_rules(d: &FeynmanDiagram, model: &HubbardModel) -> FeynmanExpression {
    // Build expression tree from diagram
}
```

**Commit:**

```bash
git commit -m "feat: add expression tree and Feynman rules"
```

---

## Task 12: Diagram Visualization (DOT output)

**Files:**
- Create: `src/visualization/dot.rs`

**Step 1: Tests, Step 3: Implement**

```rust
pub fn to_dot(diagram: &FeynmanDiagram) -> String {
    // Generate Graphviz DOT format
    // Vertices as nodes, propagators as edges
    // Spin-up: blue, spin-down: red
    // External legs: dashed
}

pub fn to_dot_all(classified: &[(FeynmanDiagram, i32)]) -> String {
    // Multiple diagrams in subgraphs
}
```

**Commit:**

```bash
git commit -m "feat: add Graphviz DOT visualization for Feynman diagrams"
```

---

## Task 13: Example Binaries

**Files:**
- Create: `examples/ladder_resummation.rs`
- Create: `examples/second_order_self_energy.rs`

End-to-end examples that use the library. The ladder example finds Tc and outputs data. The self-energy example generates, classifies, and prints diagrams.

**Commit:**

```bash
git commit -m "feat: add ladder resummation and self-energy example scripts"
```

---

## Task 14: Integration Tests

**Files:**
- Create: `tests/test_integration.rs`

End-to-end tests: Wick→Diagrams→Classification pipeline, Ladder→Tc pipeline, repulsive Hubbard has no pp instability.

**Commit:**

```bash
git commit -m "feat: add end-to-end integration tests"
```

---

## Summary

| Task | Description | Dependencies |
|------|-------------|-------------|
| 1 | Project scaffolding | — |
| 2 | FermionOperator types | 1 |
| 3 | Wick's theorem | 2 |
| 4 | Lattice + Hubbard model | 1 |
| 5 | Free Green's function | 4 |
| 6 | PP susceptibility | 5 |
| 7 | Bethe-Salpeter + Tc | 6 |
| 8 | Diagram graph types | 2 |
| 9 | Diagram generation | 3, 8 |
| 10 | Diagram classification | 9 |
| 11 | Feynman rules + Expr tree | 10 |
| 12 | DOT visualization | 8 |
| 13 | Example scripts | all |
| 14 | Integration tests | all |

**Parallel tracks:** Tasks 2→3→8→9→10→11 (algebra/diagrams) and Tasks 4→5→6→7 (numerics) can proceed in parallel after Task 1. Task 12 can start after Task 8.
