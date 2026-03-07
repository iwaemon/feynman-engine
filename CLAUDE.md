# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
cargo build                              # Build library + binary
cargo test                               # Run all tests (unit + integration, ~0.5s)
cargo test algebra::wick                  # Run tests in a specific module
cargo test --test test_integration        # Run only integration tests
cargo run --example ladder_resummation    # Run an example (may take seconds for Matsubara sums)
cargo run --example second_order_self_energy
```

## What This Project Is

A Rust perturbation theory engine for the attractive Hubbard model. It generates Feynman diagrams automatically from Wick's theorem, classifies them by topology, converts them to symbolic expressions, and performs ladder resummation to find the superconducting critical temperature Tc via the Thouless criterion.

Finite temperature (Matsubara formalism) is the default. Zero temperature is the β→∞ limit.

## Architecture

The data flows bottom-up through these layers:

```
FermionOperator → [Wick's theorem] → WickTerm (contractions + sign)
    → [generate_diagrams] → FeynmanDiagram (graph representation)
    → [classify_diagrams] → unique topologies with weights
    → [apply_feynman_rules] → FeynmanExpression (symbolic Expr tree)
```

Separately, the numerical pipeline:

```
HubbardModel + ThermalParams → evaluate_g0(k, iωₙ)
    → compute_pp_susceptibility(q, iνₘ)  [χ₀ = pp bubble]
    → solve_tmatrix / find_tc            [T = U/(1-Uχ₀), Thouless criterion]
```

### Key modules

- **`algebra/`** — `FermionOperator`, `Spin`, `Contraction`, `WickTerm`, `wick_theorem()`. Spin conservation filters matchings (↑ pairs with ↑ only). Fermion sign computed via inversion counting.
- **`models/`** — `SquareLattice` (k-grid, dispersion ε(k)=-2t(cos kx+cos ky)), `HubbardModel` (t, U), `ThermalParams` (β, μ, n_matsubara).
- **`diagrams/`** — `FeynmanDiagram` as a graph of `Vertex`/`Propagator`. `generate_diagrams()` takes `Observable` enum (Vacuum, SelfEnergy). `classify_diagrams()` uses sorted vertex-descriptor signatures for topology matching. `apply_feynman_rules()` produces `FeynmanExpression` with symbolic `Expr` tree.
- **`symbolic/`** — `Expr` enum (G0, U, Sum, Mul, Add, Inv, Scalar, Neg) with `Display` for human-readable math.
- **`resummation/`** — `PPLadder` channel. `solve_tmatrix()` computes T=U/(1-Uχ₀). `find_tc()` uses bisection on Thouless criterion 1-Uχ₀(T)=0.
- **`numerical/`** — `evaluate_g0()` for G₀(k,iωₙ)=1/(iωₙ-εₖ+μ). `compute_pp_susceptibility()` sums over k and Matsubara frequencies. Frequency index for iνₘ-iωₙ maps to fermionic index `m-1-n`.
- **`visualization/`** — `to_dot()`/`to_dot_all()` export Graphviz DOT. Spin-up=blue, spin-down=red, external legs=dashed.

## Physics Conventions

- Interaction vertex: c†↑ c↑ c†↓ c↓ at same site and time
- n-th order perturbation: (-U)^n / n! × Wick contractions of 4n operators
- PP bubble: χ₀(q,iνₘ) = -(1/Nβ) Σₖ Σₙ G₀(k,iωₙ) G₀(q-k,iνₘ-iωₙ)
- Attractive Hubbard (U<0): pp ladder diverges at Tc (Cooper instability)
- Repulsive Hubbard (U>0): no pp instability

## Design Documents

Detailed design and implementation plans are in `docs/plans/` (Japanese).
