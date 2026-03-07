# Numerical Evaluation

This chapter covers the numerical tools for evaluating Green's functions and
susceptibilities in the Matsubara (imaginary-time) formalism.

## Thermal Parameters

All finite-temperature calculations require a `ThermalParams` instance that
specifies the inverse temperature, chemical potential, and Matsubara frequency
cutoff:

```rust
use feynman_engine::prelude::*;

// beta = 10 (T = 0.1), mu = 0, 128 positive Matsubara frequencies
let params = ThermalParams::new(10.0, 0.0, 128);
```

| Field          | Type    | Description                                        |
|----------------|---------|----------------------------------------------------|
| `beta`         | `f64`   | Inverse temperature \\( \beta = 1/T \\)            |
| `mu`           | `f64`   | Chemical potential \\( \mu \\)                      |
| `n_matsubara`  | `usize` | Number of positive Matsubara frequencies to include |

The Matsubara sum runs from \\( n = -N \\) to \\( n = N-1 \\), where
\\( N \\) = `n_matsubara`, giving \\( 2N \\) frequencies in total.

## Free Green's Function

`evaluate_g0()` computes the non-interacting Matsubara Green's function:

\\[
  G_0(\mathbf{k}, i\omega_n) = \frac{1}{i\omega_n - \varepsilon_\mathbf{k} + \mu}
\\]

where the fermionic Matsubara frequency is \\( \omega_n = (2n+1)\pi/\beta \\).

```rust
use feynman_engine::prelude::*;
use std::f64::consts::PI;

let lattice = SquareLattice::new(8, 8);
let model = HubbardModel::new(lattice, 1.0, -2.0);
let params = ThermalParams::new(10.0, 0.0, 128);

// Evaluate G₀ at k = (0, 0), n = 0 (lowest positive Matsubara frequency)
let k = [0.0, 0.0];
let g0 = evaluate_g0(&model, &params, &k, 0);

// At k=(0,0): ε = -2t(cos0 + cos0) = -4t = -4
// ω₀ = π/β = π/10
// G₀ = 1/(iω₀ - ε + μ) = 1/(iπ/10 + 4)
println!("G₀(k=(0,0), iω₀) = {:.6}", g0);
println!("  Real part: {:.6}", g0.re);
println!("  Imag part: {:.6}", g0.im);
```

The function signature:

```rust,ignore
pub fn evaluate_g0(
    model: &HubbardModel,
    params: &ThermalParams,
    k: &[f64; 2],    // momentum (kx, ky)
    n: i32,           // Matsubara index (can be negative)
) -> Complex64
```

The index `n` can be any integer; negative values give \\( \omega_n < 0 \\).

## Matsubara Frequency Conventions

Fermionic and bosonic Matsubara frequencies follow distinct indexing:

| Type      | Formula                                          | Spacing          |
|-----------|--------------------------------------------------|------------------|
| Fermionic | \\( i\omega_n = i(2n+1)\pi/\beta \\)            | Odd multiples    |
| Bosonic   | \\( i\nu_m = i \cdot 2m\pi/\beta \\)            | Even multiples   |

The particle-particle bubble involves a bosonic external frequency \\( i\nu_m \\)
and two internal fermionic frequencies. The key identity is:

\\[
  i\nu_m - i\omega_n = i(2m - 2n - 1)\pi/\beta = i\omega_{m-1-n}
\\]

This means the second Green's function \\( G_0(\mathbf{q}-\mathbf{k}, i\nu_m - i\omega_n) \\)
uses fermionic index \\( m - 1 - n \\). This mapping is implemented directly in
`compute_pp_susceptibility()`.

## Particle-Particle Susceptibility

The PP bubble (pair susceptibility) is defined as:

\\[
  \chi_0(\mathbf{q}, i\nu_m) = -\frac{1}{N\beta}
    \sum_\mathbf{k} \sum_n G_0(\mathbf{k}, i\omega_n) \,
    G_0(\mathbf{q}-\mathbf{k}, i\nu_m - i\omega_n)
\\]

```rust
use feynman_engine::prelude::*;

let lattice = SquareLattice::new(8, 8);
let model = HubbardModel::new(lattice, 1.0, -2.0);
let params = ThermalParams::new(5.0, 0.0, 512);

// PP susceptibility at q = (0,0) and bosonic frequency m = 0
let chi0 = compute_pp_susceptibility(&model, &params, &[0.0, 0.0], 0);

println!("chi_0(q=0, nu=0) = {:.8} + {:.8}i", chi0.re, chi0.im);
// At q=0, the imaginary part should vanish by symmetry
assert!(chi0.im.abs() < 1e-6);
```

The function signature:

```rust,ignore
pub fn compute_pp_susceptibility(
    model: &HubbardModel,
    params: &ThermalParams,
    q: &[f64; 2],    // center-of-mass momentum
    m: i32,           // bosonic Matsubara index
) -> Complex64
```

The sum runs over all \\( N \\) k-points on the lattice and over Matsubara
indices \\( n \in [-N_\text{mat}, N_\text{mat}) \\).

## Physical Properties of the PP Bubble

At \\( \mathbf{q} = 0 \\) and \\( \nu_m = 0 \\), the PP susceptibility:

- Is **purely real** (imaginary part vanishes by time-reversal symmetry).
- Is **negative** (for the attractive Hubbard model with half-filling).
- **Grows in magnitude** as temperature decreases (\\( \beta \\) increases),
  reflecting the buildup of Cooper pair correlations.

```rust
use feynman_engine::prelude::*;

let lattice = SquareLattice::new(8, 8);
let model = HubbardModel::new(lattice, 1.0, -2.0);

let params_high = ThermalParams::new(2.0, 0.0, 256);  // T = 0.5
let params_low = ThermalParams::new(10.0, 0.0, 1024);  // T = 0.1

let chi_high = compute_pp_susceptibility(&model, &params_high, &[0.0, 0.0], 0);
let chi_low = compute_pp_susceptibility(&model, &params_low, &[0.0, 0.0], 0);

println!("|chi_0| at T=0.5: {:.6}", chi_high.re.abs());
println!("|chi_0| at T=0.1: {:.6}", chi_low.re.abs());
assert!(chi_low.re.abs() > chi_high.re.abs());
```

## Convergence with n_matsubara

The Matsubara sum must be truncated at a finite \\( N_\text{mat} \\). The
required cutoff depends on \\( \beta \\): lower temperatures need more
frequencies because the Green's function decays as \\( 1/\omega_n \\). A good
rule of thumb:

- \\( N_\text{mat} \sim 10 \times \beta \\) for qualitative results
- \\( N_\text{mat} \sim 100 \times \beta \\) for quantitative convergence

For the examples in this book, `n_matsubara = 256` to `1024` is typical.

## What Comes Next

With the PP susceptibility in hand, the next chapter shows how to use it for
ladder resummation and finding the superconducting critical temperature
\\( T_c \\). See [Finding Tc](finding-tc.md).
