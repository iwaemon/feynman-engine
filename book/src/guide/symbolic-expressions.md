# Symbolic Expressions

This chapter explains how to convert Feynman diagrams into symbolic mathematical
expressions using `apply_feynman_rules()`, and how to work with the resulting
`Expr` tree.

## Applying Feynman Rules

Given a `FeynmanDiagram` and a `HubbardModel`, `apply_feynman_rules()` produces
a `FeynmanExpression` that encodes the diagram's mathematical formula:

```rust
use feynman_engine::prelude::*;

let lattice = SquareLattice::new(4, 4);
let model = HubbardModel::new(lattice, 1.0, -2.0);

let diagrams = generate_diagrams(&model, 1, Observable::SelfEnergy);
let classified = classify_diagrams(&diagrams);

for (diag, weight) in &classified {
    let fe = apply_feynman_rules(diag, &model);
    println!("weight={}, {}", weight, fe.description);
    println!("  {}", fe.expr);
}
```

## The FeynmanExpression Structure

`apply_feynman_rules()` returns a `FeynmanExpression` with three fields:

```rust,ignore
pub struct FeynmanExpression {
    pub n_loops: usize,      // number of independent momentum loops
    pub prefactor: f64,       // (-U)^order * sign * symmetry_factor
    pub expr: Expr,           // symbolic expression tree
    pub description: String,  // human-readable summary
}
```

| Field         | Description                                                            |
|---------------|------------------------------------------------------------------------|
| `n_loops`     | Independent loop count \\( L = E_\text{internal} - V + 1 \\)          |
| `prefactor`   | Numerical prefactor \\( (-U)^n \times \text{sign} \times 1/n! \\)     |
| `expr`        | The full symbolic expression as an `Expr` tree                         |
| `description` | A string like `"Order 2: 2 loop(s), 4 internal propagators"`           |

The prefactor combines the interaction strength, fermion sign, and symmetry
factor into a single number. For example, at first order with \\( U = -2 \\):

\\[
  \text{prefactor} = (-U)^1 \times \text{sign} \times \frac{1}{1!} = 2 \times (\pm 1)
\\]

## The Expr Enum

The `Expr` enum is a recursive expression tree with these variants:

| Variant                  | Mathematical Meaning                                    | Display Example        |
|--------------------------|---------------------------------------------------------|------------------------|
| `Expr::G0 { k, omega }` | \\( G_0(\mathbf{k}, i\omega_n) \\)                     | `G₀(k1, iω1)`         |
| `Expr::U`               | Interaction \\( U \\)                                   | `U`                    |
| `Expr::Sum { var, body }`| \\( \frac{1}{N\beta} \sum_{\text{var}} [\text{body}] \\)| `Σ_{k1,ω1} [...]`     |
| `Expr::Mul(Vec<Expr>)`  | Product of expressions                                  | `(a × b × c)`         |
| `Expr::Add(Vec<Expr>)`  | Sum of expressions                                      | `(a + b + c)`          |
| `Expr::Inv(Box<Expr>)`  | Inverse \\( 1/x \\)                                     | `1/(...)`              |
| `Expr::Scalar(f64)`     | A numerical constant                                    | `2`                    |
| `Expr::Neg(Box<Expr>)`  | Negation \\( -x \\)                                     | `-(...)   `            |

The `Display` trait renders expressions in a human-readable format using Unicode
symbols (\\( G_0 \\) appears as `G₀`, summation as `Σ`).

## Interpreting the Output

The expression tree built by `apply_feynman_rules()` follows this structure:

1. **Innermost**: a product (`Mul`) of the scalar prefactor and all internal
   Green's functions \\( G_0(k_i, i\omega_i) \\).
2. **Wrapped in summations**: one `Sum` per independent loop, summing over the
   loop momentum and Matsubara frequency.

For a first-order self-energy diagram with 1 internal propagator and 1 loop:

\\[
  \Sigma^{(1)} = \sum_{k_1, \omega_1} \text{prefactor} \times G_0(k_1, i\omega_1)
\\]

For a second-order diagram with 3 internal propagators and 2 loops:

\\[
  \Sigma^{(2)} = \sum_{k_2, \omega_2} \sum_{k_1, \omega_1}
    \text{prefactor} \times G_0(k_1, i\omega_1) \times G_0(k_2, i\omega_2) \times G_0(k_3, i\omega_3)
\\]

## Complete Example: Printing All Topologies

```rust
use feynman_engine::prelude::*;

fn main() {
    let lattice = SquareLattice::new(4, 4);
    let model = HubbardModel::new(lattice, 1.0, -2.0);

    for order in 1..=2 {
        let diagrams = generate_diagrams(&model, order, Observable::SelfEnergy);
        let classified = classify_diagrams(&diagrams);

        println!("=== Order {} ({} topologies) ===", order, classified.len());
        for (i, (diag, weight)) in classified.iter().enumerate() {
            let fe = apply_feynman_rules(diag, &model);
            println!("Topology {} (weight = {}):", i + 1, weight);
            println!("  {}", fe.description);
            println!("  Prefactor: {:.4}", fe.prefactor);
            println!("  Expression: {}", fe.expr);
            println!();
        }
    }
}
```

Sample output (abbreviated):

```text
=== Order 1 (2 topologies) ===
Topology 1 (weight = 1):
  Order 1: 1 loop(s), 1 internal propagators
  Prefactor: 2.0000
  Expression: Σ_{k1,ω1} [(2 × G₀(k1, iω1))]

Topology 2 (weight = 1):
  Order 1: 1 loop(s), 1 internal propagators
  Prefactor: 2.0000
  Expression: Σ_{k1,ω1} [(2 × G₀(k1, iω1))]
```

## What Comes Next

The symbolic expressions show the structure of each diagram analytically. To
evaluate them numerically, you need the Matsubara Green's function and frequency
summation machinery covered in the [Numerical Evaluation](numerical.md) chapter.
