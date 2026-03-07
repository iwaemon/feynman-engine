# Wick's Theorem

This chapter shows how to use the `algebra` module to create fermion operators
and apply Wick's theorem to decompose time-ordered products into all possible
contractions.

## Creating Fermion Operators

A `FermionOperator` represents a single creation or annihilation operator with
a site label, spin, and imaginary-time label. Use the convenience constructors
`creator()` and `annihilator()`:

```rust
use feynman_engine::prelude::*;

// Creation operator: c†_{i,↑}(τ₁)
let c_dag = FermionOperator::creator("i", Spin::Up, "τ1");

// Annihilation operator: c_{j,↓}(τ₂)
let c = FermionOperator::annihilator("j", Spin::Down, "τ2");

assert!(c_dag.creation);   // true for c†
assert!(!c.creation);       // false for c
```

Each operator carries four fields:

| Field      | Type     | Description                              |
|------------|----------|------------------------------------------|
| `creation` | `bool`   | `true` for \\( c^\dagger \\), `false` for \\( c \\) |
| `site`     | `String` | Site label (e.g. `"i"`, `"v1"`)          |
| `spin`     | `Spin`   | `Spin::Up` or `Spin::Down`               |
| `time`     | `String` | Imaginary-time label (e.g. `"τ1"`)       |

## Applying Wick's Theorem

Pass a slice of `FermionOperator` to `wick_theorem()` to obtain all valid
contractions. The function returns a `Vec<WickTerm>`, where each term contains:

- **`contractions`**: a `Vec<Contraction>`, each pairing a creator with an
  annihilator (representing a free propagator \\( G_0 \\)).
- **`sign`**: `+1` or `-1`, the fermion sign from the permutation needed to
  bring operators into contracted pairs.

```rust
use feynman_engine::prelude::*;

// A single Hubbard interaction vertex: c†↑ c↑ c†↓ c↓
let ops = vec![
    FermionOperator::creator("i", Spin::Up, "τ1"),
    FermionOperator::annihilator("i", Spin::Up, "τ1"),
    FermionOperator::creator("i", Spin::Down, "τ1"),
    FermionOperator::annihilator("i", Spin::Down, "τ1"),
];

let terms = wick_theorem(&ops);

// With spin conservation, only 1 contraction survives
assert_eq!(terms.len(), 1);
assert_eq!(terms[0].contractions.len(), 2);
println!("Sign: {:+}", terms[0].sign);
for c in &terms[0].contractions {
    println!(
        "  G₀: c†({}, {:?}, {}) contracted with c({}, {:?}, {})",
        c.creator.site, c.creator.spin, c.creator.time,
        c.annihilator.site, c.annihilator.spin, c.annihilator.time,
    );
}
```

## Spin Conservation

A key physical constraint enforced by `wick_theorem()` is **spin conservation**:
a creation operator with \\( \text{Spin::Up} \\) can only contract with an
annihilation operator that is also \\( \text{Spin::Up} \\), and likewise for
\\( \text{Spin::Down} \\). This dramatically reduces the number of terms.

If the number of up-spin creators does not match the number of up-spin
annihilators (and similarly for down-spin), the function returns an empty vector:

```rust
use feynman_engine::prelude::*;

// Spin mismatch: one ↑ creator, one ↓ annihilator
let ops = vec![
    FermionOperator::creator("i", Spin::Up, "τ1"),
    FermionOperator::annihilator("j", Spin::Down, "τ2"),
];
let terms = wick_theorem(&ops);
assert!(terms.is_empty()); // no valid contractions
```

## Scaling: Two Interaction Vertices

For two Hubbard vertices (8 operators: 2 up-creators, 2 up-annihilators,
2 down-creators, 2 down-annihilators), there are \\( 2! \times 2! = 4 \\) Wick
contractions:

```rust
use feynman_engine::prelude::*;

let ops = vec![
    // Vertex 1
    FermionOperator::creator("v1", Spin::Up, "τ1"),
    FermionOperator::annihilator("v1", Spin::Up, "τ1"),
    FermionOperator::creator("v1", Spin::Down, "τ1"),
    FermionOperator::annihilator("v1", Spin::Down, "τ1"),
    // Vertex 2
    FermionOperator::creator("v2", Spin::Up, "τ2"),
    FermionOperator::annihilator("v2", Spin::Up, "τ2"),
    FermionOperator::creator("v2", Spin::Down, "τ2"),
    FermionOperator::annihilator("v2", Spin::Down, "τ2"),
];

let terms = wick_theorem(&ops);
assert_eq!(terms.len(), 4);

// Each term has 4 contractions (one per creator-annihilator pair)
for (i, term) in terms.iter().enumerate() {
    println!("Term {}: sign = {:+}, {} contractions",
             i, term.sign, term.contractions.len());
}
```

## Odd Operator Count

An odd number of operators always yields zero terms, as required by physics:

```rust
use feynman_engine::prelude::*;

let ops = vec![FermionOperator::creator("i", Spin::Up, "τ1")];
assert!(wick_theorem(&ops).is_empty());
```

## What Comes Next

The Wick contractions produced here are the raw input for the diagram generation
pipeline. The next chapter shows how `generate_diagrams()` converts these
contractions into `FeynmanDiagram` graph objects and classifies them by topology.
