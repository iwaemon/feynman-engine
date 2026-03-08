# feynman-engine

[日本語ガイド](docs/usage-ja.md)

A perturbation theory engine for the Hubbard model in Rust.

## Overview

`feynman-engine` automates the perturbative expansion of the attractive Hubbard model at finite temperature. Starting from fermion operators, it applies Wick's theorem to generate all contractions with proper fermion signs and spin conservation, constructs Feynman diagrams, classifies them by topology, and converts them into symbolic expressions via Feynman rules. On the numerical side, it evaluates the particle-particle susceptibility in the Matsubara formalism and performs ladder resummation (Bethe-Salpeter equation) to locate the superconducting critical temperature Tc through the Thouless criterion.

## Features

- Automatic Wick contraction with spin conservation and fermion sign
- Feynman diagram generation and topological classification
- Symbolic expression tree with Feynman rules
- Particle-particle ladder resummation (Bethe-Salpeter)
- Thouless criterion Tc finder (bisection)
- Graphviz DOT visualization (spin-colored propagators)

## Installation

### Prerequisites

- [Rust](https://www.rust-lang.org/tools/install) toolchain (edition 2021, MSRV 1.70)
- [Graphviz](https://graphviz.org/) (optional, for rendering DOT diagrams to PNG/SVG)
  - macOS: `brew install graphviz`
  - Ubuntu: `apt install graphviz`

```bash
git clone <repo-url>
cd feynman
cargo build
cargo test
```

## Quick Start

### Ladder resummation

Find the superconducting Tc for the attractive Hubbard model on a square lattice:

```bash
cargo run --example ladder_resummation
```

Output (trimmed):

```
=== Ladder Resummation: Attractive Hubbard Model ===
Lattice: 8x8 square, t = 1, U = -2

Searching for Tc via Thouless criterion (n_matsubara = 256)...
Tc = 0.223253

=== Thouless Criterion: 1 - U * chi_0(q=0, nu=0) ===
T            chi_0(q=0)           1 - U*chi_0
----------------------------------------------------
0.178603     -0.56792543          -0.13585086
0.200928     -0.53070584          -0.06141169
0.212091     -0.51465041          -0.02930081
0.223253     -0.49996228          0.00007544
...
```

At T = Tc the Thouless criterion 1 - U*chi_0 vanishes, signaling the onset of superconducting instability.

### Self-energy diagrams

Generate and classify self-energy diagrams up to second order:

```bash
cargo run --example second_order_self_energy
```

Output (trimmed):

```
=== Attractive Hubbard Model ===
Lattice: 4x4 square, t = 1, U = -2

1st order: 2 Wick contractions
2nd order: 12 Wick contractions

1st order: 2 unique topologies
2nd order: 7 unique topologies

=== 1st Order Self-Energy Diagrams ===
Topology 1 (weight = 1):
  Order 1: 2 loop(s), 2 internal propagators
  Prefactor: 2.0000
  Expression: Σ_k1,ω1 [Σ_k2,ω2 [(2 × G₀(k1, iω1) × G₀(k2, iω2))]]

Topology 2 (weight = -1):
  Order 1: 1 loop(s), 1 internal propagators
  Prefactor: -2.0000
  Expression: Σ_k1,ω1 [(-2 × G₀(k1, iω1))]

... 7 more topologies for 2nd order
```

#### Exporting diagrams

Export Graphviz DOT and/or JSON (for D3.js visualization) to files:

```bash
# Export DOT and render to PNG
cargo run --example second_order_self_energy -- --dot diagrams.dot
dot -Tpng diagrams.dot -o diagrams.png

# Export JSON for D3.js
cargo run --example second_order_self_energy -- --json diagrams.json

# Export both at once
cargo run --example second_order_self_energy -- --dot diagrams.dot --json diagrams.json

# Write DOT to stdout (suppresses text summary)
cargo run --example second_order_self_energy -- --dot -
```

Spin-up propagators are blue, spin-down are red, and external legs are dashed in the DOT output.

The JSON output contains an array of classified diagrams with `topology_id`, `weight`, `vertices`, `propagators` (with `from`, `to`, `spin`, `external`), `sign`, and `symmetry_factor` fields.

## Architecture

The data flows bottom-up through two pipelines.

**Symbolic pipeline** (diagram generation):

```
FermionOperator → [Wick's theorem] → WickTerm (contractions + sign)
    → [generate_diagrams] → FeynmanDiagram (graph representation)
    → [classify_diagrams] → unique topologies with weights
    → [apply_feynman_rules] → FeynmanExpression (symbolic Expr tree)
```

**Numerical pipeline** (resummation):

```
HubbardModel + ThermalParams → evaluate_g0(k, iωₙ)
    → compute_pp_susceptibility(q, iνₘ)  [χ₀ = pp bubble]
    → solve_tmatrix / find_tc            [T = U/(1-Uχ₀), Thouless criterion]
```

### Module map

| Module | Description |
|---|---|
| `algebra/` | Fermion operators, spin, contractions, Wick's theorem |
| `models/` | Square lattice (k-grid, dispersion), Hubbard model parameters |
| `diagrams/` | Feynman diagram graph, generation, topological classification, Feynman rules |
| `symbolic/` | Expression tree (`G0`, `U`, `Sum`, `Mul`, `Add`, ...) with `Display` |
| `resummation/` | PP ladder channel, T-matrix, Thouless criterion Tc finder |
| `numerical/` | Green's function G0(k, iwn), Matsubara sums, PP susceptibility |
| `visualization/` | Graphviz DOT export (spin-up = blue, spin-down = red) |

## Physics Background

### The Hubbard model

The Hubbard model describes electrons hopping on a lattice with an on-site interaction. On a 2D square lattice the dispersion is epsilon(k) = -2t(cos kx + cos ky), and the interaction vertex couples four fermion operators at the same site: c_up^dagger c_up c_down^dagger c_down with coupling U. For U < 0 (attractive), the model supports s-wave superconductivity; for U > 0 (repulsive), the particle-particle channel does not diverge.

### Matsubara formalism

At finite temperature T = 1/beta, time-ordered perturbation theory is formulated in imaginary time. Fermionic Matsubara frequencies are iwn = i(2n+1)pi/beta. The bare Green's function is G0(k, iwn) = 1/(iwn - epsilon_k + mu). All internal momentum and frequency sums run over the Brillouin zone and the discrete Matsubara set.

### Particle-particle susceptibility and ladder resummation

The particle-particle bubble (bare susceptibility) is chi_0(q, i_nu_m) = -1/(N*beta) sum_k sum_n G0(k, iwn) G0(q-k, i_nu_m - iwn). Ladder resummation sums repeated particle-particle scattering to all orders via the T-matrix: T(q, i_nu_m) = U / (1 - U * chi_0(q, i_nu_m)). This is equivalent to solving the Bethe-Salpeter equation in the particle-particle channel.

### Thouless criterion

The superconducting instability occurs when the T-matrix diverges, i.e., when 1 - U * chi_0(q=0, i_nu_0=0) = 0. This is the Thouless criterion. The critical temperature Tc is found by bisection: scanning temperature until the criterion is satisfied. For the attractive Hubbard model (U < 0), chi_0 is negative, so U*chi_0 = |U|*|chi_0| > 0 and the criterion becomes 1 - |U|*|chi_0| = 0. The divergence signals the Cooper instability and onset of pair condensation.

## Documentation

This project uses [mdBook](https://rust-lang.github.io/mdBook/) for detailed documentation covering theory, usage guides, and visualization.

### Building the docs

```bash
cargo install mdbook    # install mdBook (once)

mdbook build book/      # build English docs -> target/book/
mdbook build book-ja/   # build Japanese docs -> target/book-ja/
```

### Previewing locally

```bash
mdbook serve book/      # English: http://localhost:3000
mdbook serve book-ja/   # Japanese: http://localhost:3000
```

### Online documentation

The documentation is published on GitHub Pages:

- [English](https://iwaemon.github.io/feynman-engine/)
- [日本語](https://iwaemon.github.io/feynman-engine/ja/)

### Documentation structure

| Section | Content |
|---|---|
| **Theory** | Hubbard model, perturbation theory, Wick's theorem, Feynman diagrams, ladder resummation |
| **Getting Started** | Installation, quick start with examples |
| **Guide** | Module-by-module usage with code examples |
| **Visualization** | Graphviz DOT export, D3.js interactive diagrams |
