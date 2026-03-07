# The Hubbard Model

The Hubbard model is one of the simplest lattice models that captures the competition between kinetic energy (electron hopping) and electron-electron interactions. This chapter describes the model as implemented in the `models/` module of feynman-engine.

## Hamiltonian

The Hubbard Hamiltonian on a lattice is:

\\[
H = -t \sum_{\langle i,j \rangle, \sigma} c^\dagger_{i\sigma} c_{j\sigma} + U \sum_i n_{i\uparrow} n_{i\downarrow}
\\]

where:

- \\( c^\dagger_{i\sigma} \\) and \\( c_{i\sigma} \\) are creation and annihilation operators for an electron with spin \\( \sigma \\) at site \\( i \\),
- \\( n_{i\sigma} = c^\dagger_{i\sigma} c_{i\sigma} \\) is the number operator,
- \\( t \\) is the hopping amplitude between nearest-neighbor sites \\( \langle i, j \rangle \\),
- \\( U \\) is the on-site interaction strength.

## Square Lattice and Dispersion

feynman-engine works on a 2D square lattice. In momentum space, the tight-binding dispersion relation is:

\\[
\varepsilon(\mathbf{k}) = -2t(\cos k_x + \cos k_y)
\\]

This dispersion has a bandwidth of \\( 8t \\), ranging from \\( -4t \\) at \\( \mathbf{k} = (0, 0) \\) (the \\( \Gamma \\) point) to \\( +4t \\) at \\( \mathbf{k} = (\pi, \pi) \\) (the M point). At half-filling, the Fermi surface is a perfect square connecting the \\( (\pi, 0) \\) and \\( (0, \pi) \\) points (van Hove singularity), leading to a logarithmic divergence in the density of states.

The `SquareLattice` struct (in `src/models/lattice.rs`) stores the lattice dimensions `lx` and `ly`, provides a `k_grid()` method that generates the discrete Brillouin zone points \\( k_\alpha = 2\pi n_\alpha / L_\alpha \\), and a `dispersion(k, t)` method that evaluates \\( \varepsilon(\mathbf{k}) \\).

## Attractive vs. Repulsive Interaction

The physics depends crucially on the sign of \\( U \\):

- **Attractive (\\( U < 0 \\))**: Electrons experience an effective attraction at the same site. This drives s-wave superconductivity via Cooper pair formation. The particle-particle (pp) channel develops an instability at sufficiently low temperature, and the T-matrix diverges at \\( T_c \\). This is the primary focus of feynman-engine.

- **Repulsive (\\( U > 0 \\))**: The on-site repulsion disfavors double occupancy. The particle-particle ladder does not diverge (no Cooper instability in the pp channel). Magnetic instabilities (antiferromagnetism) may occur instead, but these are not currently implemented.

## Interaction Vertex

In the perturbative expansion, each order of perturbation theory introduces an interaction vertex. The Hubbard interaction couples four fermion operators at the same site and imaginary time:

\\[
U \, c^\dagger_{i\uparrow}(\tau) \, c_{i\uparrow}(\tau) \, c^\dagger_{i\downarrow}(\tau) \, c_{i\downarrow}(\tau)
\\]

This vertex structure enforces spin conservation: an up-spin creation operator must pair with an up-spin annihilation operator, and similarly for down-spin. This constraint dramatically reduces the number of allowed Wick contractions.

## Implementation

The `HubbardModel` struct (in `src/models/hubbard.rs`) holds a `SquareLattice`, the hopping parameter `t`, and the interaction strength `u`. It provides a `dispersion(k)` method that delegates to the lattice.

Finite-temperature calculations require the `ThermalParams` struct, which stores:

- `beta`: inverse temperature \\( \beta = 1/T \\),
- `mu`: chemical potential \\( \mu \\),
- `n_matsubara`: the number of Matsubara frequencies used in summations.

A typical setup looks like:

```rust
use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::{HubbardModel, ThermalParams};

let lattice = SquareLattice::new(8, 8);          // 8x8 square lattice
let model = HubbardModel::new(lattice, 1.0, -2.0); // t=1, U=-2 (attractive)
let params = ThermalParams::new(10.0, 0.0, 512);   // beta=10, mu=0, 512 frequencies
```
