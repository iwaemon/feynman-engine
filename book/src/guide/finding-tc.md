# Finding Tc

This chapter shows how to use particle-particle ladder resummation to locate the
superconducting critical temperature \\( T_c \\) via the Thouless criterion.

## The T-Matrix

The `solve_tmatrix()` function computes the particle-particle T-matrix from the
ladder resummation of the Bethe-Salpeter equation:

\\[
  T(\mathbf{q}, i\nu_m) = \frac{U}{1 - U \, \chi_0(\mathbf{q}, i\nu_m)}
\\]

where \\( \chi_0 \\) is the PP bubble computed by `compute_pp_susceptibility()`.

```rust
use feynman_engine::prelude::*;

let lattice = SquareLattice::new(8, 8);
let model = HubbardModel::new(lattice, 1.0, -2.0);
let params = ThermalParams::new(5.0, 0.0, 512);

let t_val = solve_tmatrix(&model, &params, &[0.0, 0.0], 0);
println!("T-matrix at q=0, nu=0: {:.6}", t_val);
```

The function signature:

```rust,ignore
pub fn solve_tmatrix(
    model: &HubbardModel,
    params: &ThermalParams,
    q: &[f64; 2],    // center-of-mass momentum
    m: i32,           // bosonic Matsubara index
) -> Complex64
```

## The Thouless Criterion

The T-matrix diverges when its denominator vanishes:

\\[
  1 - U \, \chi_0(\mathbf{q}=0, i\nu_0=0;\, T) = 0
\\]

This is the **Thouless criterion**: the temperature at which this condition is
satisfied is the superconducting critical temperature \\( T_c \\). For the
attractive Hubbard model (\\( U < 0 \\)), the PP susceptibility
\\( \chi_0 < 0 \\), so \\( U \chi_0 > 0 \\), and the criterion can be
satisfied.

For repulsive \\( U > 0 \\), the product \\( U \chi_0 < 0 \\), so the
denominator \\( 1 - U\chi_0 > 1 \\) and no instability occurs in the PP channel.

## Finding Tc by Bisection

`find_tc()` uses bisection to locate the temperature where the Thouless
criterion crosses zero:

```rust
use feynman_engine::prelude::*;

let lattice = SquareLattice::new(8, 8);
let model = HubbardModel::new(lattice, 1.0, -2.0);

// Search for Tc in the range [0.01, 2.0] with tolerance 1e-4
let tc = find_tc(&model, (0.01, 2.0), 256, 1e-4);
println!("Tc = {:.6}", tc);
```

The function signature:

```rust,ignore
pub fn find_tc(
    model: &HubbardModel,
    t_range: (f64, f64),     // (T_low, T_high) temperature bracket
    n_matsubara: usize,       // Matsubara cutoff for each evaluation
    tol: f64,                 // convergence tolerance on T
) -> f64
```

The bracket `(T_low, T_high)` must satisfy:

- At \\( T_\text{low} \\): \\( 1 - U\chi_0 < 0 \\) (below \\( T_c \\), the
  criterion has crossed zero).
- At \\( T_\text{high} \\): \\( 1 - U\chi_0 > 0 \\) (above \\( T_c \\), no
  instability).

If the bracket does not straddle zero, the function will panic with a
descriptive error message.

## Full Workflow Example

This follows the `ladder_resummation` example:

```rust
use feynman_engine::prelude::*;

fn main() {
    // 1. Set up the model
    let lattice = SquareLattice::new(8, 8);
    let model = HubbardModel::new(lattice, 1.0, -2.0);
    println!("Lattice: 8x8 square, t = {}, U = {}", model.t, model.u);

    // 2. Find Tc
    let tc = find_tc(&model, (0.01, 2.0), 256, 1e-4);
    println!("Tc = {:.6}", tc);

    // 3. Scan the Thouless criterion across temperatures
    println!("{:<12} {:<20} {:<20}", "T", "chi_0(q=0)", "1 - U*chi_0");
    for factor in &[0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.5, 2.0] {
        let temp = tc * factor;
        let beta = 1.0 / temp;
        let params = ThermalParams::new(beta, 0.0, 256);
        let chi0 = compute_pp_susceptibility(
            &model, &params, &[0.0, 0.0], 0,
        );
        let thouless = 1.0 - model.u * chi0.re;
        println!("{:<12.6} {:<20.8} {:<20.8}", temp, chi0.re, thouless);
    }
}
```

Expected behavior:

- Below \\( T_c \\): the Thouless criterion \\( 1 - U\chi_0 \\) is **negative**
  (the T-matrix has already diverged).
- At \\( T_c \\): the criterion passes through **zero**.
- Above \\( T_c \\): the criterion is **positive** (normal state, no
  instability).

## Interpreting Results

The critical temperature \\( T_c \\) from ladder resummation corresponds to the
**BCS mean-field** \\( T_c \\) for the Hubbard model. Some key points:

- \\( T_c \\) increases with \\( |U| \\) (stronger attraction means stronger
  pairing).
- \\( T_c \\) depends on the lattice size through the density of states. Larger
  lattices give results closer to the thermodynamic limit.
- At half-filling (\\( \mu = 0 \\)), the van Hove singularity in the 2D square
  lattice enhances pairing.

## Tips for Convergence

### Choosing n_matsubara

The Matsubara cutoff affects both accuracy and runtime. Each evaluation of
\\( \chi_0 \\) sums over \\( 2 N_\text{mat} \\) frequencies for each of the
\\( N \\) k-points:

| `n_matsubara` | Use case                         | Typical runtime (8x8) |
|---------------|----------------------------------|-----------------------|
| 128           | Quick exploratory scan           | Fast                  |
| 256           | Reasonable accuracy for \\( T_c \\) | Moderate           |
| 512--1024     | Publication-quality results      | Slower                |

Lower temperatures require larger `n_matsubara` for convergence because more
frequencies carry significant spectral weight.

### Choosing the temperature bracket

- Start with a wide bracket like `(0.01, 2.0)`.
- If bisection panics, the bracket does not straddle \\( T_c \\). Try widening
  the range or checking that the model has an attractive interaction (\\( U < 0 \\)).
- For very weak coupling \\( |U| \ll t \\), \\( T_c \\) can be exponentially
  small, requiring a very low `T_low`.

### Tolerance

A tolerance of `1e-4` gives 4 significant digits in \\( T_c \\). The bisection
converges geometrically, so halving the tolerance costs only one extra iteration.
