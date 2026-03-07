# Quick Start

This chapter walks through the two built-in examples to give you a feel for
what feynman-engine can do.

## Ladder Resummation -- Finding \\( T_c \\)

The `ladder_resummation` example sets up an 8x8 square lattice with hopping
\\( t = 1 \\) and attractive interaction \\( U = -2 \\), then uses bisection on the
Thouless criterion to locate the superconducting critical temperature.

```bash
cargo run --example ladder_resummation
```

Expected output:

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

The table shows the particle-particle susceptibility \\( \chi_0(q{=}0, i\nu_0) \\)
and the Thouless criterion \\( 1 - U \chi_0 \\) at several temperatures around
\\( T_c \\). At \\( T = T_c \approx 0.2233 \\), the criterion crosses zero,
signaling the Cooper instability and onset of s-wave superconductivity.

**What this demonstrates:**

- Building a `HubbardModel` on a `SquareLattice`
- Computing the particle-particle susceptibility via Matsubara sums
- Finding \\( T_c \\) with the built-in bisection solver (`find_tc`)

## Self-Energy Diagrams

The `second_order_self_energy` example generates all first- and second-order
self-energy diagrams from Wick's theorem, classifies them by topology, and
prints the corresponding symbolic expressions from Feynman rules.

```bash
cargo run --example second_order_self_energy
```

Expected output (trimmed):

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

=== 2nd Order Self-Energy Diagrams ===
...
```

**What this demonstrates:**

- Automatic Wick contraction with spin conservation and fermion sign
- Feynman diagram generation and topological classification
- Symbolic expression tree construction via Feynman rules

## Exporting Diagrams

Both DOT (Graphviz) and JSON (for D3.js visualization) export are available
via command-line flags on the `second_order_self_energy` example:

```bash
# Export Graphviz DOT and render to PNG
cargo run --example second_order_self_energy -- --dot diagrams.dot
dot -Tpng diagrams.dot -o diagrams.png

# Export JSON for D3.js visualization
cargo run --example second_order_self_energy -- --json diagrams.json

# Export both at once
cargo run --example second_order_self_energy -- --dot diagrams.dot --json diagrams.json

# Write DOT to stdout (suppresses the text summary)
cargo run --example second_order_self_energy -- --dot -
```

In the DOT output, spin-up propagators are colored blue, spin-down are red,
and external legs are drawn with dashed lines.
