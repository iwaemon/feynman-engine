# Ladder Resummation & Tc

This chapter describes the particle-particle ladder resummation and the Thouless criterion for finding the superconducting critical temperature \\( T_c \\). The implementation is in the `resummation/` and `numerical/` modules.

## Particle-Particle Susceptibility

The bare particle-particle (pp) bubble is the lowest-order contribution to the two-particle correlation function in the Cooper channel. It is defined as:

\\[
\chi_0(\mathbf{q}, i\nu_m) = -\frac{1}{N\beta} \sum_\mathbf{k} \sum_n G_0(\mathbf{k}, i\omega_n) \, G_0(\mathbf{q}-\mathbf{k}, i\nu_m - i\omega_n)
\\]

where:

- \\( G_0(\mathbf{k}, i\omega_n) = 1/(i\omega_n - \varepsilon_\mathbf{k} + \mu) \\) is the bare Green's function,
- \\( i\omega_n \\) are fermionic Matsubara frequencies,
- \\( i\nu_m = i \cdot 2m\pi/\beta \\) are bosonic Matsubara frequencies (the total pair frequency),
- \\( N \\) is the number of lattice sites.

The product \\( G_0(\mathbf{k}, i\omega_n) \, G_0(\mathbf{q}-\mathbf{k}, i\nu_m - i\omega_n) \\) represents a pair of electrons with total momentum \\( \mathbf{q} \\) and total frequency \\( i\nu_m \\).

The function `compute_pp_susceptibility(model, params, q, m)` in `src/numerical/matsubara.rs` evaluates this sum numerically. A key implementation detail is the frequency index mapping: the second Green's function uses the fermionic index \\( m - 1 - n \\) (not \\( m - n \\)), which correctly maps the bosonic frequency \\( i\nu_m - i\omega_n \\) to a fermionic Matsubara frequency.

At \\( \mathbf{q} = 0 \\) and \\( i\nu_0 = 0 \\), the susceptibility \\( \chi_0 \\) is purely real (by symmetry) and negative. Its magnitude \\( |\chi_0| \\) increases as the temperature decreases (i.e., as \\( \beta \\) increases).

## T-Matrix (Ladder Resummation)

The particle-particle ladder resummation sums repeated scattering of a Cooper pair to all orders in \\( U \\). The resulting T-matrix is:

\\[
T(\mathbf{q}, i\nu_m) = \frac{U}{1 - U \chi_0(\mathbf{q}, i\nu_m)}
\\]

This is the solution to the Bethe-Salpeter equation in the particle-particle channel with a local (momentum-independent) interaction \\( U \\). Diagrammatically, it sums the geometric series:

\\[
T = U + U \chi_0 U + U \chi_0 U \chi_0 U + \cdots = \frac{U}{1 - U \chi_0}
\\]

The function `solve_tmatrix(model, params, q, m)` in `src/resummation/bethe_salpeter.rs` evaluates this expression. The `PPLadder` struct in `src/resummation/channels.rs` represents the particle-particle channel.

## Thouless Criterion

The superconducting instability occurs when the T-matrix diverges, i.e., when its denominator vanishes:

\\[
1 - U \chi_0(\mathbf{q} = 0, i\nu_0 = 0) = 0
\\]

This is the **Thouless criterion**. It signals the onset of Cooper pair condensation at zero total pair momentum and zero pair frequency (the most favorable channel for s-wave superconductivity).

### Why the Attractive Case Works

For the attractive Hubbard model (\\( U < 0 \\)):

- The pp susceptibility satisfies \\( \chi_0 < 0 \\) (real and negative at \\( \mathbf{q} = 0 \\)).
- Therefore \\( U \chi_0 = |U| \cdot |\chi_0| > 0 \\).
- The Thouless criterion becomes \\( 1 - |U| \cdot |\chi_0| = 0 \\), or equivalently \\( |\chi_0| = 1/|U| \\).

As temperature decreases, \\( |\chi_0| \\) grows. At high temperature, \\( |\chi_0| < 1/|U| \\) and the system is in the normal state. At \\( T = T_c \\), \\( |\chi_0| \\) reaches \\( 1/|U| \\) and the instability occurs.

For the repulsive case (\\( U > 0 \\)):

- \\( U \chi_0 = U \cdot \chi_0 < 0 \\) (since \\( \chi_0 < 0 \\)).
- Therefore \\( 1 - U \chi_0 > 1 \\) for all temperatures.
- The T-matrix never diverges: there is no Cooper instability.

## Bisection Algorithm for Tc

The function `find_tc(model, t_range, n_matsubara, tol)` in `src/resummation/bethe_salpeter.rs` finds \\( T_c \\) using bisection on the Thouless criterion function:

\\[
f(T) = 1 - U \chi_0(\mathbf{q}=0, i\nu_0=0; T)
\\]

The algorithm:

1. **Bracket**: Provide a temperature range \\( (T_{\text{lo}}, T_{\text{hi}}) \\) such that \\( f(T_{\text{lo}}) \\) and \\( f(T_{\text{hi}}) \\) have opposite signs. That is, \\( T_{\text{lo}} < T_c < T_{\text{hi}} \\).

2. **Bisect**: Repeatedly evaluate \\( f \\) at the midpoint \\( T_{\text{mid}} = (T_{\text{lo}} + T_{\text{hi}})/2 \\). If \\( f(T_{\text{mid}}) \\) has the same sign as \\( f(T_{\text{lo}}) \\), update \\( T_{\text{lo}} \leftarrow T_{\text{mid}} \\); otherwise update \\( T_{\text{hi}} \leftarrow T_{\text{mid}} \\).

3. **Converge**: Stop when \\( T_{\text{hi}} - T_{\text{lo}} < \text{tol} \\). Return the midpoint as \\( T_c \\).

At each temperature evaluation, a new `ThermalParams` is constructed with \\( \beta = 1/T \\) and the pp susceptibility is recomputed.

### Example

For the attractive Hubbard model on an 8x8 square lattice with \\( t = 1 \\), \\( U = -2 \\), at half-filling (\\( \mu = 0 \\)):

```rust
use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::HubbardModel;
use feynman_engine::resummation::bethe_salpeter::find_tc;

let lattice = SquareLattice::new(8, 8);
let model = HubbardModel::new(lattice, 1.0, -2.0);
let tc = find_tc(&model, (0.01, 2.0), 512, 1e-4);
// tc ≈ 0.223
```

The critical temperature \\( T_c \approx 0.223t \\) marks the onset of s-wave superconductivity. Above \\( T_c \\), the system is a normal metal; below \\( T_c \\), Cooper pairs condense and the system becomes superconducting.

## Convergence Considerations

The accuracy of the \\( T_c \\) estimate depends on two numerical parameters:

- **Lattice size** (\\( L \times L \\)): Larger lattices give better momentum resolution. Finite-size effects are most significant near van Hove singularities.

- **Number of Matsubara frequencies** (`n_matsubara`): The frequency sum must be converged. Near \\( T_c \\), where \\( \beta \\) is large, more frequencies are needed. Typical values range from 256 to 1024.
