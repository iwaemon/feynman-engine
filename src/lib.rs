//! # feynman-engine
//!
//! A perturbation theory engine for the attractive Hubbard model at finite temperature.
//!
//! This crate automates the perturbative expansion: starting from fermion operators, it
//! applies Wick's theorem to generate all contractions with proper fermion signs and spin
//! conservation, constructs Feynman diagrams, classifies them by topology, and converts
//! them into symbolic expressions via Feynman rules. On the numerical side, it evaluates
//! the particle-particle susceptibility in the Matsubara formalism and performs ladder
//! resummation to locate the superconducting critical temperature Tc.
//!
//! ## Data flow
//!
//! **Symbolic pipeline** (diagram generation):
//! ```text
//! FermionOperator → [Wick's theorem] → WickTerm (contractions + sign)
//!     → [generate_diagrams] → FeynmanDiagram (graph representation)
//!     → [classify_diagrams] → unique topologies with weights
//!     → [apply_feynman_rules] → FeynmanExpression (symbolic Expr tree)
//! ```
//!
//! **Numerical pipeline** (resummation):
//! ```text
//! HubbardModel + ThermalParams → evaluate_g0(k, iωₙ)
//!     → compute_pp_susceptibility(q, iνₘ)
//!     → find_tc (Thouless criterion, bisection)
//! ```
//!
//! ## Quick start
//!
//! ```rust
//! use feynman_engine::prelude::*;
//!
//! // Generate 1st-order self-energy diagrams
//! let lattice = SquareLattice::new(4, 4);
//! let model = HubbardModel::new(lattice, 1.0, -2.0);
//! let diagrams = generate_diagrams(&model, 1, Observable::SelfEnergy);
//! let classified = classify_diagrams(&diagrams);
//! for (diag, weight) in &classified {
//!     let expr = apply_feynman_rules(diag, &model);
//!     println!("weight={}, expr={}", weight, expr.expr);
//! }
//! ```

pub mod algebra;
pub mod models;
pub mod diagrams;
pub mod symbolic;
pub mod resummation;
pub mod numerical;
pub mod visualization;

/// Convenience re-exports of the most commonly used types and functions.
pub mod prelude {
    pub use crate::algebra::operators::{FermionOperator, Spin, Contraction, WickTerm};
    pub use crate::algebra::wick::wick_theorem;
    pub use crate::models::lattice::SquareLattice;
    pub use crate::models::hubbard::{HubbardModel, ThermalParams};
    pub use crate::diagrams::graph::{FeynmanDiagram, Vertex, Propagator};
    pub use crate::diagrams::generate::{generate_diagrams, Observable};
    pub use crate::diagrams::classify::classify_diagrams;
    pub use crate::diagrams::rules::{apply_feynman_rules, FeynmanExpression};
    pub use crate::symbolic::expr::Expr;
    pub use crate::numerical::greens_function::evaluate_g0;
    pub use crate::numerical::matsubara::compute_pp_susceptibility;
    pub use crate::resummation::bethe_salpeter::{find_tc, solve_tmatrix};
    pub use crate::visualization::dot::{to_dot, to_dot_all};
    pub use crate::visualization::json::to_json_all;
}
