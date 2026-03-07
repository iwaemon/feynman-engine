# Diagram Generation

This chapter covers creating a Hubbard model, generating Feynman diagrams from
Wick contractions, and classifying them into unique topologies.

## Setting Up the Model

Start by creating a square lattice and a Hubbard model. The lattice defines the
momentum grid, and the model adds the hopping parameter \\( t \\) and on-site
interaction \\( U \\):

```rust
use feynman_engine::prelude::*;

// 4x4 square lattice with t=1, U=-2 (attractive)
let lattice = SquareLattice::new(4, 4);
let model = HubbardModel::new(lattice, 1.0, -2.0);

println!("Sites: {}", model.lattice.num_sites()); // 16
println!("t = {}, U = {}", model.t, model.u);     // t = 1, U = -2
```

`SquareLattice::new(lx, ly)` creates an \\( L_x \times L_y \\) lattice with
tight-binding dispersion \\( \varepsilon(\mathbf{k}) = -2t(\cos k_x + \cos k_y) \\).

## Generating Diagrams

`generate_diagrams()` takes the model, perturbation order, and an `Observable`
enum specifying which correlation function to expand:

```rust
use feynman_engine::prelude::*;

let lattice = SquareLattice::new(4, 4);
let model = HubbardModel::new(lattice, 1.0, -2.0);

// Generate all 1st-order self-energy diagrams
let diagrams_1st = generate_diagrams(&model, 1, Observable::SelfEnergy);
println!("1st order: {} Wick contractions", diagrams_1st.len());

// Generate all 2nd-order self-energy diagrams
let diagrams_2nd = generate_diagrams(&model, 2, Observable::SelfEnergy);
println!("2nd order: {} Wick contractions", diagrams_2nd.len());
```

The `Observable` enum has two variants:

| Variant                    | Description                                                  |
|----------------------------|--------------------------------------------------------------|
| `Observable::SelfEnergy`   | Adds external legs (one \\( c^\dagger_\uparrow \\) in, one \\( c_\uparrow \\) out) |
| `Observable::Vacuum`       | No external legs; generates vacuum (bubble) diagrams         |

Internally, `generate_diagrams()` constructs the operator string for the given
order (4 operators per interaction vertex, plus 2 for external legs if
self-energy), calls `wick_theorem()`, and converts each `WickTerm` into a
`FeynmanDiagram`.

## The FeynmanDiagram Structure

Each diagram is a graph with vertices and propagators:

```rust,ignore
pub struct FeynmanDiagram {
    pub order: usize,              // perturbation order
    pub vertices: Vec<Vertex>,     // interaction vertices
    pub propagators: Vec<Propagator>, // fermion lines (Green's functions)
    pub sign: i32,                 // fermion sign (+1 or -1)
    pub symmetry_factor: f64,      // 1/n! prefactor
}
```

- **`Vertex`** has an `id`, `site` label, and `time` label.
- **`Propagator`** connects two vertex IDs with a `spin` and an `external` flag.
  External propagators represent the incoming/outgoing legs of the self-energy.
- **`count_loops()`** computes the number of independent momentum loops:
  \\( L = E_\text{internal} - V + 1 \\).

## Classifying by Topology

Many Wick contractions yield the same Feynman diagram topology (they differ only
in labeling). `classify_diagrams()` groups them and sums the fermion signs:

```rust
use feynman_engine::prelude::*;

let lattice = SquareLattice::new(4, 4);
let model = HubbardModel::new(lattice, 1.0, -2.0);

let diagrams = generate_diagrams(&model, 1, Observable::SelfEnergy);
let classified = classify_diagrams(&diagrams);

println!("{} unique topologies from {} contractions",
         classified.len(), diagrams.len());

for (i, (diag, weight)) in classified.iter().enumerate() {
    println!("Topology {}: weight = {}, loops = {}, propagators = {}",
             i + 1, weight, diag.count_loops(),
             diag.propagators.len());
}
```

The return type is `Vec<(FeynmanDiagram, i32)>`:

- The `FeynmanDiagram` is one representative of the topology class.
- The `i32` weight is the sum of fermion signs across all Wick terms that share
  this topology. Topologies whose signs cancel to zero are automatically removed.

For first-order self-energy, you get exactly **2 topologies**:
1. The **Hartree** (tadpole) diagram -- a self-loop on the vertex.
2. The **Fock** (exchange) diagram -- propagator goes to the vertex and back.

## Full Example

This mirrors the `second_order_self_energy` example:

```rust
use feynman_engine::prelude::*;

fn main() {
    let lattice = SquareLattice::new(4, 4);
    let model = HubbardModel::new(lattice, 1.0, -2.0);

    // Generate and classify 1st and 2nd order
    let diag_1 = generate_diagrams(&model, 1, Observable::SelfEnergy);
    let diag_2 = generate_diagrams(&model, 2, Observable::SelfEnergy);
    let class_1 = classify_diagrams(&diag_1);
    let class_2 = classify_diagrams(&diag_2);

    println!("1st order: {} contractions -> {} topologies",
             diag_1.len(), class_1.len());
    println!("2nd order: {} contractions -> {} topologies",
             diag_2.len(), class_2.len());

    for (diag, weight) in &class_1 {
        println!("  weight={:+}, order={}, loops={}",
                 weight, diag.order, diag.count_loops());
    }
}
```

## What Comes Next

With classified diagrams in hand, the next step is to convert them into symbolic
mathematical expressions using Feynman rules. See the
[Symbolic Expressions](symbolic-expressions.md) chapter.
