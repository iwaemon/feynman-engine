# Perturbation Theory & Wick's Theorem

This chapter describes the finite-temperature perturbation theory framework used by feynman-engine, implemented in the `algebra/` module.

## Matsubara Formalism

At finite temperature \\( T = 1/\beta \\), perturbation theory is formulated in imaginary time \\( \tau \in [0, \beta) \\). Fourier transforming to frequency space, the fermionic Matsubara frequencies are:

\\[
i\omega_n = i\frac{(2n+1)\pi}{\beta}, \quad n \in \mathbb{Z}
\\]

These odd multiples of \\( \pi/\beta \\) enforce the anti-periodicity \\( G(\tau + \beta) = -G(\tau) \\) required by Fermi-Dirac statistics. Bosonic quantities (such as the particle-particle susceptibility) use even frequencies \\( i\nu_m = i \cdot 2m\pi/\beta \\).

In the code, `matsubara_freq(n, beta)` (in `src/numerical/greens_function.rs`) computes \\( \omega_n = (2n+1)\pi/\beta \\).

## Bare Green's Function

The non-interacting (bare) Matsubara Green's function is:

\\[
G_0(\mathbf{k}, i\omega_n) = \frac{1}{i\omega_n - \varepsilon_\mathbf{k} + \mu}
\\]

where \\( \varepsilon_\mathbf{k} \\) is the lattice dispersion and \\( \mu \\) is the chemical potential. This propagator is the building block of all diagrammatic calculations.

The function `evaluate_g0(model, params, k, n)` (in `src/numerical/greens_function.rs`) evaluates this expression numerically as a complex number.

## Perturbative Expansion

The \\( n \\)-th order contribution to any observable involves the time-ordered product of \\( n \\) interaction vertices. Each Hubbard interaction vertex contributes four fermion operators \\( c^\dagger_\uparrow c_\uparrow c^\dagger_\downarrow c_\downarrow \\), so an \\( n \\)-th order term involves \\( 4n \\) operators (plus any external legs). The combinatorial prefactor is:

\\[
\frac{(-U)^n}{n!}
\\]

The factor \\( (-1)^n \\) comes from expanding the exponential \\( e^{-\int_0^\beta H_I(\tau) d\tau} \\), and the \\( 1/n! \\) accounts for the time ordering of \\( n \\) identical vertices.

## Wick's Theorem

Wick's theorem states that the expectation value of a time-ordered product of fermion operators in a non-interacting state equals the sum over all complete pairings (contractions) of the operators, with each pairing weighted by a fermion sign.

A **contraction** pairs a creation operator \\( c^\dagger \\) with an annihilation operator \\( c \\), representing a bare propagator \\( G_0 \\). For \\( 2m \\) operators, the number of complete pairings is \\( (2m-1)!! \\) in general, but spin conservation dramatically reduces this.

### Spin Conservation

The Hubbard interaction vertex couples \\( c^\dagger_\uparrow c_\uparrow c^\dagger_\downarrow c_\downarrow \\) at the same site. Since \\( G_0 \\) is diagonal in spin, Wick contractions must respect spin conservation:

- A spin-\\( \uparrow \\) creation operator can only contract with a spin-\\( \uparrow \\) annihilation operator.
- A spin-\\( \downarrow \\) creation operator can only contract with a spin-\\( \downarrow \\) annihilation operator.

This means the pairing problem factorizes by spin. For each spin species, the number of creators must equal the number of annihilators. The implementation groups operators by spin and generates permutations independently for each spin sector, then takes the Cartesian product.

### Fermion Sign

Each complete pairing carries a sign \\( (-1)^P \\), where \\( P \\) is the parity of the permutation that rearranges the operators into contracted pairs. The implementation computes this by counting inversions: given the sequence of original operator indices rearranged into \\( (c_1, a_1, c_2, a_2, \ldots) \\) pairs, the sign is \\( (-1)^{\text{inversions}} \\).

For example, with two Hubbard vertices (8 operators, 4 each of \\( \uparrow \\) and \\( \downarrow \\)), there are \\( 2! \times 2! = 4 \\) Wick contractions, each with 4 propagators and a definite sign.

## Implementation

The `algebra/` module provides the following types and functions:

- **`Spin`** (`operators.rs`): An enum with variants `Up` and `Down`.
- **`FermionOperator`** (`operators.rs`): Represents \\( c^\dagger_{i\sigma}(\tau) \\) or \\( c_{i\sigma}(\tau) \\), with fields `creation` (bool), `site`, `spin`, and `time`.
- **`Contraction`** (`operators.rs`): A pair of a creation and an annihilation operator, representing a single \\( G_0 \\) propagator.
- **`WickTerm`** (`operators.rs`): A complete set of contractions with an associated fermion `sign` (\\( \pm 1 \\)).
- **`wick_theorem(ops)`** (`wick.rs`): Takes a slice of `FermionOperator`s and returns all valid `WickTerm`s. It enforces spin conservation by grouping operators by spin, generates all permutations within each spin sector, and computes the fermion sign via inversion counting.

```rust
use feynman_engine::algebra::operators::*;
use feynman_engine::algebra::wick::wick_theorem;

// Single Hubbard vertex: c†↑ c↑ c†↓ c↓
let ops = vec![
    FermionOperator::creator("i", Spin::Up, "τ1"),
    FermionOperator::annihilator("i", Spin::Up, "τ1"),
    FermionOperator::creator("i", Spin::Down, "τ1"),
    FermionOperator::annihilator("i", Spin::Down, "τ1"),
];
let terms = wick_theorem(&ops);
// With spin conservation, only 1 contraction is possible
assert_eq!(terms.len(), 1);
```
