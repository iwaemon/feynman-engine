# Feynman Diagrams

This chapter explains how feynman-engine generates, classifies, and translates Feynman diagrams into symbolic mathematical expressions. The implementation lives in the `diagrams/` and `symbolic/` modules.

## From Wick Contractions to Diagrams

Each Wick contraction (a complete pairing of creation and annihilation operators) corresponds to a Feynman diagram. The function `generate_diagrams(model, order, observable)` in `src/diagrams/generate.rs` automates this process:

1. **Construct the operator string.** For the chosen `Observable`, the appropriate external legs are prepended. Then, for each of the `order` interaction vertices, four operators \\( c^\dagger_\uparrow c_\uparrow c^\dagger_\downarrow c_\downarrow \\) are appended, each labeled with a unique site and imaginary-time variable.

2. **Apply Wick's theorem.** The full operator string is passed to `wick_theorem()`, which returns all valid contractions with their fermion signs.

3. **Build graph representation.** Each Wick term is converted into a `FeynmanDiagram` by mapping vertices (interaction points) and propagators (contractions) onto a graph structure. The symmetry factor \\( 1/n! \\) is attached.

## Observable Types

The `Observable` enum determines which diagrams are generated:

- **`Vacuum`**: No external legs. All operators come from interaction vertices. These are the closed (vacuum) diagrams contributing to the thermodynamic potential.

- **`SelfEnergy`**: Two external legs (one incoming creator \\( c^\dagger_\uparrow \\), one outgoing annihilator \\( c_\uparrow \\)) are added to the operator string. The resulting diagrams are one-particle irreducible self-energy corrections \\( \Sigma(\mathbf{k}, i\omega_n) \\).

## Graph Representation

The diagram graph is defined in `src/diagrams/graph.rs`:

- **`Vertex`**: Represents an interaction point, with fields `id` (unique integer), `site` (lattice site label), and `time` (imaginary-time label).

- **`Propagator`**: Represents a bare Green's function \\( G_0 \\) line, with fields `from` and `to` (vertex IDs), `spin` (\\( \uparrow \\) or \\( \downarrow \\)), and `external` (whether the line connects to an external leg).

- **`FeynmanDiagram`**: The complete diagram, containing a list of vertices, a list of propagators, the perturbation `order`, the fermion `sign`, and the `symmetry_factor`.

The number of independent momentum loops is computed by the `count_loops()` method using the Euler formula for connected graphs:

\\[
L = E_{\text{internal}} - V + 1
\\]

where \\( E_{\text{internal}} \\) is the number of internal propagators and \\( V \\) is the number of vertices.

## Topological Classification

Many distinct Wick contractions produce diagrams with the same topology (the same graph structure up to relabeling of vertices). The function `classify_diagrams(diagrams)` in `src/diagrams/classify.rs` groups diagrams by topology and sums their weights.

The classification works by computing a **canonical signature** for each diagram. For every vertex, a descriptor tuple is computed:

\\[
(\text{in}\_\uparrow, \; \text{out}\_\uparrow, \; \text{in}\_\downarrow, \; \text{out}\_\downarrow, \; \text{self\\_loops})
\\]

counting the number of incoming and outgoing propagators of each spin, plus the number of self-loops. The list of all vertex descriptors is sorted to produce a permutation-invariant signature. Two diagrams with the same sorted signature are considered topologically equivalent.

For example, at first order in the self-energy, there are 2 Wick contractions that reduce to 2 unique topologies: the Hartree (tadpole) diagram and the Fock (exchange) diagram.

## Feynman Rules

The function `apply_feynman_rules(diagram, model)` in `src/diagrams/rules.rs` converts a `FeynmanDiagram` into a `FeynmanExpression` containing a symbolic expression tree. The rules are:

1. **Prefactor**: \\( (-U)^n \times \text{sign} \times \text{symmetry\\_factor} \\), where \\( n \\) is the perturbation order.

2. **Propagators**: Each internal propagator contributes a factor of \\( G_0(k_i, i\omega_i) \\).

3. **Loop summations**: Each independent loop contributes a summation \\( \sum_{k, \omega} \\) over internal momentum and Matsubara frequency.

The result is a `FeynmanExpression` with fields:
- `n_loops`: number of independent momentum/frequency loops,
- `prefactor`: the numerical prefactor,
- `expr`: a symbolic `Expr` tree,
- `description`: a human-readable summary.

## Symbolic Expression Tree

The `Expr` enum in the `symbolic/` module represents mathematical expressions as a tree:

| Variant | Meaning |
|---------|---------|
| `G0 { k, omega }` | Bare Green's function \\( G_0(k, i\omega) \\) |
| `U` | Interaction vertex |
| `Scalar(f64)` | Numerical constant |
| `Mul(Vec<Expr>)` | Product of expressions |
| `Add(Vec<Expr>)` | Sum of expressions |
| `Neg(Box<Expr>)` | Negation |
| `Inv(Box<Expr>)` | Inverse \\( 1/(\cdots) \\) |
| `Sum { var, body }` | Summation \\( \sum_{\text{var}} \text{body} \\) |

The `Display` trait is implemented to produce human-readable output, for example:

```text
Σ_k1,ω1 [Σ_k2,ω2 [(2 × G₀(k1, iω1) × G₀(k2, iω2))]]
```

This symbolic representation provides a bridge between the combinatorial diagram generation and the numerical evaluation pipeline.
