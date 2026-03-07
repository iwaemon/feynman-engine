# Introduction

## What feynman-engine does

`feynman-engine` is a Rust library that automates perturbative many-body calculations for the attractive Hubbard model at finite temperature. Starting from bare fermion operators, it applies Wick's theorem to generate all contractions with correct fermion signs and spin conservation, constructs Feynman diagrams, classifies them by topology, and converts them into symbolic mathematical expressions via standard Feynman rules.

On the numerical side, the library evaluates the particle-particle susceptibility in the Matsubara formalism and performs ladder resummation (solving the Bethe-Salpeter equation) to locate the superconducting critical temperature Tc through the Thouless criterion.

## Two pipelines

The engine is organized into two complementary pipelines.

**Symbolic pipeline** (diagram generation):

```text
FermionOperator -> [Wick's theorem] -> WickTerm (contractions + sign)
    -> [generate_diagrams] -> FeynmanDiagram (graph representation)
    -> [classify_diagrams] -> unique topologies with weights
    -> [apply_feynman_rules] -> FeynmanExpression (symbolic Expr tree)
```

**Numerical pipeline** (resummation):

```text
HubbardModel + ThermalParams -> evaluate_g0(k, iwn)
    -> compute_pp_susceptibility(q, inu_m)  [chi_0 = pp bubble]
    -> solve_tmatrix / find_tc              [T = U/(1-U*chi_0), Thouless criterion]
```

## Module map

| Module | Description |
|---|---|
| `algebra/` | Fermion operators, spin, contractions, Wick's theorem |
| `models/` | Square lattice (k-grid, dispersion), Hubbard model parameters |
| `diagrams/` | Feynman diagram graph, generation, topological classification, Feynman rules |
| `symbolic/` | Expression tree (`G0`, `U`, `Sum`, `Mul`, `Add`, ...) with `Display` |
| `resummation/` | PP ladder channel, T-matrix, Thouless criterion Tc finder |
| `numerical/` | Green's function G0(k, iwn), Matsubara sums, PP susceptibility |
| `visualization/` | Graphviz DOT and D3.js JSON export |

## Who this book is for

This book is intended for two audiences: **developers** who want to use `feynman-engine` as a library for automated diagram generation and numerical many-body calculations, and **physics students and researchers** who want to understand how perturbation theory, Wick's theorem, and ladder resummation are implemented in practice. No prior Rust experience is assumed for the physics chapters, and no prior many-body theory is assumed for the code walkthrough chapters.
