# FeynmanEngine.jl Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build a Julia symbolic perturbation engine that auto-generates Feynman diagrams via Wick's theorem and resums ladder diagrams for the attractive Hubbard model to find Tc.

**Architecture:** Bottom-up from operator algebra (Wick's theorem) through diagram generation/classification, to Bethe-Salpeter resummation. Matsubara formalism for finite temperature. Symbolics.jl for symbolic expressions, Graphs.jl for diagram topology, Makie.jl for visualization.

**Tech Stack:** Julia, Symbolics.jl, Graphs.jl, CairoMakie.jl, FFTW.jl

---

## Task 1: Project Scaffolding

**Files:**
- Create: `Project.toml`
- Create: `src/FeynmanEngine.jl`
- Create: `test/runtests.jl`

**Step 1: Initialize Julia package**

```bash
cd /Users/shumpei/work/feynman
julia -e '
using Pkg
Pkg.generate("FeynmanEngine")
'
```

This creates `Project.toml` and `src/FeynmanEngine.jl`. Move generated files to project root:

```bash
mv FeynmanEngine/Project.toml .
mv FeynmanEngine/src/FeynmanEngine.jl src/FeynmanEngine.jl
rm -r FeynmanEngine
```

**Step 2: Add dependencies**

```bash
julia --project=. -e '
using Pkg
Pkg.add(["Symbolics", "Graphs", "CairoMakie", "FFTW", "Combinatorics"])
Pkg.add("Test")  # stdlib, but ensure available
'
```

**Step 3: Create directory structure**

```bash
mkdir -p src/algebra src/diagrams src/resummation src/models src/numerical src/visualization
mkdir -p test examples
```

**Step 4: Set up main module file**

Write `src/FeynmanEngine.jl`:

```julia
module FeynmanEngine

# Core algebra
include("algebra/operators.jl")
include("algebra/wick.jl")

# Models
include("models/lattice.jl")
include("models/hubbard.jl")

# Diagrams
include("diagrams/graph.jl")
include("diagrams/generate.jl")
include("diagrams/classify.jl")
include("diagrams/rules.jl")

# Resummation
include("resummation/channels.jl")
include("resummation/bethe_salpeter.jl")

# Numerical evaluation
include("numerical/greens_function.jl")
include("numerical/matsubara.jl")

# Visualization
include("visualization/diagram_plot.jl")

end # module
```

Create placeholder files (empty modules) for each include so the package loads:

```bash
for f in src/algebra/operators.jl src/algebra/wick.jl \
         src/models/lattice.jl src/models/hubbard.jl \
         src/diagrams/graph.jl src/diagrams/generate.jl \
         src/diagrams/classify.jl src/diagrams/rules.jl \
         src/resummation/channels.jl src/resummation/bethe_salpeter.jl \
         src/numerical/greens_function.jl src/numerical/matsubara.jl \
         src/visualization/diagram_plot.jl; do
    echo "# placeholder" > "$f"
done
```

**Step 5: Create test runner**

Write `test/runtests.jl`:

```julia
using Test
using FeynmanEngine

@testset "FeynmanEngine" begin
    @test true  # package loads
end
```

**Step 6: Run tests to verify package loads**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: All tests pass.

**Step 7: Commit**

```bash
git add -A
git commit -m "feat: scaffold FeynmanEngine.jl package with directory structure"
```

---

## Task 2: Fermion Operator Types

**Files:**
- Create: `src/algebra/operators.jl`
- Create: `test/algebra/test_operators.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/algebra/test_operators.jl`:

```julia
using Test
using FeynmanEngine

@testset "FermionOperator" begin
    # Basic construction
    c_dag = FermionOperator(true, :i, :↑, :τ₁)
    c = FermionOperator(false, :i, :↑, :τ₁)

    @test c_dag.creation == true
    @test c.creation == false
    @test c_dag.site == :i
    @test c_dag.spin == :↑
    @test c_dag.time == :τ₁

    # Convenience constructors
    c_dag2 = creator(:j, :↓, :τ₂)
    c2 = annihilator(:j, :↓, :τ₂)
    @test c_dag2.creation == true
    @test c2.creation == false
    @test c2.site == :j
    @test c2.spin == :↓
end

@testset "Contraction" begin
    c_dag = creator(:i, :↑, :τ₁)
    c = annihilator(:j, :↑, :τ₂)
    cont = Contraction(c_dag, c)

    @test cont.creator.site == :i
    @test cont.annihilator.site == :j
    # Spin must match for valid contraction
    @test cont.creator.spin == cont.annihilator.spin
end

@testset "WickTerm" begin
    c1 = creator(:i, :↑, :τ₁)
    a1 = annihilator(:j, :↑, :τ₂)
    c2 = creator(:i, :↓, :τ₁)
    a2 = annihilator(:j, :↓, :τ₂)

    term = WickTerm(1, [Contraction(c1, a1), Contraction(c2, a2)])
    @test term.sign == 1
    @test length(term.contractions) == 2
end
```

Add to `test/runtests.jl`:

```julia
using Test
using FeynmanEngine

@testset "FeynmanEngine" begin
    include("algebra/test_operators.jl")
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `FermionOperator` not defined.

**Step 3: Write implementation**

Write `src/algebra/operators.jl`:

```julia
export FermionOperator, Contraction, WickTerm, creator, annihilator

"""
Fermionic creation or annihilation operator c†_{site,spin}(τ) or c_{site,spin}(τ).
"""
struct FermionOperator
    creation::Bool
    site::Symbol
    spin::Symbol
    time::Symbol
end

"""Convenience constructor for creation operator c†."""
creator(site::Symbol, spin::Symbol, time::Symbol) = FermionOperator(true, site, spin, time)

"""Convenience constructor for annihilation operator c."""
annihilator(site::Symbol, spin::Symbol, time::Symbol) = FermionOperator(false, site, spin, time)

"""
Wick contraction between a creator and an annihilator.
Represents the free propagator G₀.
"""
struct Contraction
    creator::FermionOperator
    annihilator::FermionOperator
end

"""
A single term from Wick's theorem: a product of contractions with a sign.
"""
struct WickTerm
    sign::Int
    contractions::Vector{Contraction}
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/algebra/operators.jl test/algebra/test_operators.jl test/runtests.jl
git commit -m "feat: add FermionOperator, Contraction, WickTerm types"
```

---

## Task 3: Wick's Theorem

**Files:**
- Create: `src/algebra/wick.jl`
- Create: `test/algebra/test_wick.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/algebra/test_wick.jl`:

```julia
using Test
using FeynmanEngine

@testset "Wick's theorem - 2 operators" begin
    # ⟨T c†_{i↑}(τ₁) c_{j↑}(τ₂)⟩ = G₀(i,τ₁; j,τ₂) — one contraction
    ops = [creator(:i, :↑, :τ₁), annihilator(:j, :↑, :τ₂)]
    terms = wick_theorem(ops)

    @test length(terms) == 1
    @test terms[1].sign == -1  # one swap from normal ordering: c†c → -G₀
    @test terms[1].contractions[1].creator.site == :i
    @test terms[1].contractions[1].annihilator.site == :j
end

@testset "Wick's theorem - 4 operators, same spin" begin
    # ⟨T c†_{i↑} c_{j↑} c†_{k↑} c_{l↑}⟩ = 2 terms
    ops = [
        creator(:i, :↑, :τ₁), annihilator(:j, :↑, :τ₂),
        creator(:k, :↑, :τ₃), annihilator(:l, :↑, :τ₄)
    ]
    terms = wick_theorem(ops)

    @test length(terms) == 2
    # Check signs sum: should have one +1 and one -1
    signs = sort([t.sign for t in terms])
    @test signs == [-1, 1]
end

@testset "Wick's theorem - 4 operators, different spins" begin
    # Hubbard vertex: c†_{i↑} c_{i↑} c†_{i↓} c_{i↓}
    # Spin conservation: ↑ pairs with ↑, ↓ pairs with ↓ → only 1 term
    ops = [
        creator(:i, :↑, :τ₁), annihilator(:i, :↑, :τ₁),
        creator(:i, :↓, :τ₁), annihilator(:i, :↓, :τ₁)
    ]
    terms = wick_theorem(ops)

    @test length(terms) == 1
    @test length(terms[1].contractions) == 2
end

@testset "Wick's theorem - 8 operators (2nd order Hubbard)" begin
    # Two Hubbard vertices: 8 operators, 2 creators per spin, 2 annihilators per spin
    ops = [
        creator(:i, :↑, :τ₁), annihilator(:i, :↑, :τ₁),
        creator(:i, :↓, :τ₁), annihilator(:i, :↓, :τ₁),
        creator(:j, :↑, :τ₂), annihilator(:j, :↑, :τ₂),
        creator(:j, :↓, :τ₂), annihilator(:j, :↓, :τ₂),
    ]
    terms = wick_theorem(ops)

    # 2 creators ↑ × 2 annihilators ↑ → 2! = 2 matchings
    # 2 creators ↓ × 2 annihilators ↓ → 2! = 2 matchings
    # Total: 2 × 2 = 4 terms
    @test length(terms) == 4
    @test all(t -> length(t.contractions) == 4, terms)
end

@testset "Wick's theorem - odd number of operators" begin
    ops = [creator(:i, :↑, :τ₁)]
    terms = wick_theorem(ops)
    @test length(terms) == 0
end

@testset "Wick's theorem - spin mismatch" begin
    # One ↑ creator, one ↓ annihilator — no valid contractions
    ops = [creator(:i, :↑, :τ₁), annihilator(:j, :↓, :τ₂)]
    terms = wick_theorem(ops)
    @test length(terms) == 0
end
```

Add to `test/runtests.jl`:

```julia
@testset "FeynmanEngine" begin
    include("algebra/test_operators.jl")
    include("algebra/test_wick.jl")
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `wick_theorem` not defined.

**Step 3: Write implementation**

Write `src/algebra/wick.jl`:

```julia
using Combinatorics

export wick_theorem

"""
    wick_theorem(ops::Vector{FermionOperator}) -> Vector{WickTerm}

Apply Wick's theorem to a time-ordered product of fermion operators.
Returns all fully-contracted terms with correct fermion signs.
Enforces spin conservation (G₀ is spin-diagonal).
"""
function wick_theorem(ops::Vector{FermionOperator})::Vector{WickTerm}
    n = length(ops)
    # Odd number of operators → zero
    if n % 2 != 0
        return WickTerm[]
    end

    creators = [(i, op) for (i, op) in enumerate(ops) if op.creation]
    annihilators = [(i, op) for (i, op) in enumerate(ops) if !op.creation]

    # Must have equal numbers
    if length(creators) != length(annihilators)
        return WickTerm[]
    end

    # Group by spin for spin conservation
    spins = unique(op.spin for (_, op) in creators)
    creators_by_spin = Dict(s => [(i, op) for (i, op) in creators if op.spin == s] for s in spins)
    annihilators_by_spin = Dict(s => [(i, op) for (i, op) in annihilators if op.spin == s] for s in spins)

    # Check spin balance
    for s in spins
        nc = length(get(creators_by_spin, s, []))
        na = length(get(annihilators_by_spin, s, []))
        if nc != na
            return WickTerm[]
        end
    end

    # Check annihilator spins have no extra spins not in creators
    ann_spins = unique(op.spin for (_, op) in annihilators)
    for s in ann_spins
        if !haskey(creators_by_spin, s)
            return WickTerm[]
        end
    end

    # Generate all matchings: for each spin, permute annihilators
    # and pair with creators in order
    spin_matchings = Dict{Symbol, Vector{Vector{Tuple{Int,Int}}}}()
    for s in spins
        cs = creators_by_spin[s]
        as = annihilators_by_spin[s]
        matchings = Vector{Tuple{Int,Int}}[]
        for perm in permutations(1:length(as))
            matching = [(cs[k][1], as[perm[k]][1]) for k in 1:length(cs)]
            push!(matchings, matching)
        end
        spin_matchings[s] = matchings
    end

    # Cartesian product over spins
    sorted_spins = sort(collect(spins))
    all_matchings = _cartesian_product([spin_matchings[s] for s in sorted_spins])

    terms = WickTerm[]
    for combined in all_matchings
        # Flatten all (creator_idx, annihilator_idx) pairs
        pairs = Tuple{Int,Int}[]
        for matching in combined
            append!(pairs, matching)
        end

        # Build contractions
        contractions = [Contraction(ops[ci], ops[ai]) for (ci, ai) in pairs]

        # Compute fermion sign from the permutation
        sign = _fermion_sign(ops, pairs)

        push!(terms, WickTerm(sign, contractions))
    end

    return terms
end

"""
Compute the sign of the permutation that brings operators into contracted pairs.
"""
function _fermion_sign(ops::Vector{FermionOperator}, pairs::Vector{Tuple{Int,Int}})::Int
    n = length(ops)
    # Build the permutation: original positions → paired order
    # The target ordering is: (c1, a1, c2, a2, ...) where (ci, ai) are pairs
    target = Int[]
    for (ci, ai) in pairs
        push!(target, ci)
        push!(target, ai)
    end
    # Count transpositions to sort target into this order from 1:n
    return _permutation_sign(target)
end

"""
Compute sign of a permutation by counting inversions.
"""
function _permutation_sign(perm::Vector{Int})::Int
    n = length(perm)
    inversions = 0
    for i in 1:n
        for j in (i+1):n
            if perm[i] > perm[j]
                inversions += 1
            end
        end
    end
    return iseven(inversions) ? 1 : -1
end

"""
Cartesian product of vectors of vectors.
"""
function _cartesian_product(lists::Vector{Vector{Vector{Tuple{Int,Int}}}})
    if isempty(lists)
        return [Vector{Tuple{Int,Int}}[]]
    end
    result = [Vector{Tuple{Int,Int}}[]]
    for list in lists
        new_result = Vector{Vector{Vector{Tuple{Int,Int}}}}()
        for prev in result
            for item in list
                push!(new_result, vcat(prev, [item]))
            end
        end
        result = new_result
    end
    return result
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/algebra/wick.jl test/algebra/test_wick.jl test/runtests.jl
git commit -m "feat: implement Wick's theorem with spin conservation"
```

---

## Task 4: Lattice and Hubbard Model

**Files:**
- Create: `src/models/lattice.jl`
- Create: `src/models/hubbard.jl`
- Create: `test/models/test_lattice.jl`
- Create: `test/models/test_hubbard.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/models/test_lattice.jl`:

```julia
using Test
using FeynmanEngine

@testset "SquareLattice" begin
    lat = SquareLattice(4, 4)
    @test lat.Lx == 4
    @test lat.Ly == 4
    @test num_sites(lat) == 16

    # k-grid: should give 16 k-points
    kg = k_grid(lat)
    @test length(kg) == 16

    # Dispersion ε(k) = -2t(cos kx + cos ky)
    # At Γ point k=(0,0): ε = -4t
    @test dispersion(lat, [0.0, 0.0], 1.0) ≈ -4.0
    # At M point k=(π,π): ε = 4t
    @test dispersion(lat, [π, π], 1.0) ≈ 4.0 atol=1e-12
end
```

Write `test/models/test_hubbard.jl`:

```julia
using Test
using FeynmanEngine

@testset "HubbardModel" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    @test model.t == 1.0
    @test model.U == -2.0
    @test model.lattice.Lx == 4
end

@testset "ThermalParameters" begin
    params = ThermalParameters(β=10.0, μ=0.0, n_matsubara=128)
    @test params.β == 10.0
    @test params.μ == 0.0
    @test params.n_matsubara == 128
end
```

Add to `test/runtests.jl`:

```julia
include("models/test_lattice.jl")
include("models/test_hubbard.jl")
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `SquareLattice` not defined.

**Step 3: Write implementation**

Write `src/models/lattice.jl`:

```julia
export SquareLattice, num_sites, k_grid, dispersion

struct SquareLattice
    Lx::Int
    Ly::Int
end

num_sites(lat::SquareLattice) = lat.Lx * lat.Ly

function k_grid(lat::SquareLattice)
    kpoints = Vector{Float64}[]
    for ix in 0:(lat.Lx-1), iy in 0:(lat.Ly-1)
        kx = 2π * ix / lat.Lx
        ky = 2π * iy / lat.Ly
        push!(kpoints, [kx, ky])
    end
    return kpoints
end

function dispersion(lat::SquareLattice, k::Vector{Float64}, t::Float64)
    return -2t * (cos(k[1]) + cos(k[2]))
end
```

Write `src/models/hubbard.jl`:

```julia
export HubbardModel, ThermalParameters

struct HubbardModel
    lattice::SquareLattice
    t::Float64
    U::Float64
end

struct ThermalParameters
    β::Float64
    μ::Float64
    n_matsubara::Int
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/models/ test/models/ test/runtests.jl
git commit -m "feat: add SquareLattice, HubbardModel, ThermalParameters"
```

---

## Task 5: Free Green's Function (Numerical)

**Files:**
- Create: `src/numerical/greens_function.jl`
- Create: `test/numerical/test_greens_function.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/numerical/test_greens_function.jl`:

```julia
using Test
using FeynmanEngine

@testset "Free Green's function G₀" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)
    params = ThermalParameters(β=10.0, μ=0.0, n_matsubara=128)

    # G₀(k, iωₙ) = 1/(iωₙ - εₖ + μ)
    k = [0.0, 0.0]  # Γ point, εₖ = -4
    n = 0            # ω₀ = π/β

    g = evaluate_G0(model, params, k, n)
    ω₀ = π / params.β
    εk = dispersion(lat, k, model.t)
    expected = 1.0 / (im * ω₀ - εk + params.μ)

    @test g ≈ expected
end

@testset "G₀ sum rule" begin
    # (1/β) Σₙ G₀(k, iωₙ) e^{iωₙ 0⁺} = n_F(εₖ - μ)
    # For εₖ = 0 and μ = 0: n_F(0) = 0.5
    lat = SquareLattice(8, 8)
    model = HubbardModel(lat, t=1.0, U=0.0)
    params = ThermalParameters(β=20.0, μ=0.0, n_matsubara=2048)

    k = [π/2, π/2]  # εₖ = 0 for this k
    εk = dispersion(lat, k, model.t)
    @test abs(εk) < 1e-10  # sanity check

    # Sum over Matsubara frequencies
    nk = 0.0
    for n in -params.n_matsubara:params.n_matsubara-1
        g = evaluate_G0(model, params, k, n)
        nk += real(g) / params.β
    end
    # Should give Fermi function n_F(0) = 0.5
    @test nk ≈ 0.5 atol=0.01
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `evaluate_G0` not defined.

**Step 3: Write implementation**

Write `src/numerical/greens_function.jl`:

```julia
export evaluate_G0, matsubara_frequency

"""Fermionic Matsubara frequency ωₙ = (2n+1)π/β."""
function matsubara_frequency(n::Int, β::Float64)
    return (2n + 1) * π / β
end

"""
    evaluate_G0(model, params, k, n) -> ComplexF64

Free Matsubara Green's function G₀(k, iωₙ) = 1/(iωₙ - εₖ + μ).
"""
function evaluate_G0(model::HubbardModel, params::ThermalParameters,
                     k::Vector{Float64}, n::Int)::ComplexF64
    ωₙ = matsubara_frequency(n, params.β)
    εₖ = dispersion(model.lattice, k, model.t)
    return 1.0 / (im * ωₙ - εₖ + params.μ)
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/numerical/greens_function.jl test/numerical/test_greens_function.jl test/runtests.jl
git commit -m "feat: add free Matsubara Green's function G₀(k, iωₙ)"
```

---

## Task 6: Matsubara Summation and PP Susceptibility

**Files:**
- Create: `src/numerical/matsubara.jl`
- Create: `test/numerical/test_matsubara.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/numerical/test_matsubara.jl`:

```julia
using Test
using FeynmanEngine

@testset "PP susceptibility χ₀" begin
    lat = SquareLattice(8, 8)
    model = HubbardModel(lat, t=1.0, U=-2.0)
    params = ThermalParameters(β=5.0, μ=0.0, n_matsubara=512)

    # χ₀(q=0, iν₀=0) at half-filling should be real and negative
    # (for pp channel with our sign convention)
    q = [0.0, 0.0]
    m = 0  # bosonic frequency index
    χ₀ = compute_pp_susceptibility(model, params, q, m)

    @test imag(χ₀) ≈ 0.0 atol=1e-6  # should be real at q=0, ν=0
    @test real(χ₀) < 0.0             # negative for pp bubble
end

@testset "PP susceptibility increases with β" begin
    lat = SquareLattice(8, 8)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    # |χ₀| should increase as temperature decreases (β increases)
    q = [0.0, 0.0]
    params_high_T = ThermalParameters(β=2.0, μ=0.0, n_matsubara=256)
    params_low_T = ThermalParameters(β=10.0, μ=0.0, n_matsubara=1024)

    χ_high = abs(real(compute_pp_susceptibility(model, params_high_T, q, 0)))
    χ_low = abs(real(compute_pp_susceptibility(model, params_low_T, q, 0)))

    @test χ_low > χ_high
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `compute_pp_susceptibility` not defined.

**Step 3: Write implementation**

Write `src/numerical/matsubara.jl`:

```julia
export compute_pp_susceptibility, bosonic_matsubara_frequency

"""Bosonic Matsubara frequency νₘ = 2mπ/β."""
function bosonic_matsubara_frequency(m::Int, β::Float64)
    return 2m * π / β
end

"""
    compute_pp_susceptibility(model, params, q, m) -> ComplexF64

Particle-particle bubble:
χ₀(q, iνₘ) = -(1/Nβ) Σₖ Σₙ G₀(k, iωₙ) G₀(q-k, iνₘ - iωₙ)
"""
function compute_pp_susceptibility(model::HubbardModel, params::ThermalParameters,
                                    q::Vector{Float64}, m::Int)::ComplexF64
    kpoints = k_grid(model.lattice)
    N = length(kpoints)
    χ₀ = 0.0 + 0.0im

    for k in kpoints
        qmk = q .- k
        for n in -params.n_matsubara:params.n_matsubara-1
            G1 = evaluate_G0(model, params, k, n)
            G2 = evaluate_G0(model, params, qmk, m - 1 - n)  # iνₘ - iωₙ maps to index
            χ₀ += G1 * G2
        end
    end

    return -χ₀ / (N * params.β)
end
```

**Note:** The index for the second Green's function needs care. If iνₘ = 2mπ/β and iωₙ = (2n+1)π/β, then iνₘ - iωₙ = (2(m-n) - 1)π/β corresponds to fermionic index `m - 1 - n`. This must be verified during implementation — the test's physics checks (real, correct sign, temperature dependence) will catch errors.

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS. If the frequency index mapping is wrong, the "should be real" test will fail — adjust accordingly.

**Step 5: Commit**

```bash
git add src/numerical/matsubara.jl test/numerical/test_matsubara.jl test/runtests.jl
git commit -m "feat: add particle-particle susceptibility χ₀(q, iνₘ)"
```

---

## Task 7: Bethe-Salpeter Equation and Tc Finder

**Files:**
- Create: `src/resummation/channels.jl`
- Create: `src/resummation/bethe_salpeter.jl`
- Create: `test/resummation/test_bethe_salpeter.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/resummation/test_bethe_salpeter.jl`:

```julia
using Test
using FeynmanEngine

@testset "T-matrix" begin
    lat = SquareLattice(8, 8)
    model = HubbardModel(lat, t=1.0, U=-2.0)
    params = ThermalParameters(β=5.0, μ=0.0, n_matsubara=512)

    # T(q, iνₘ) = U / (1 - U·χ₀(q, iνₘ))
    bse = PPLadder(model)
    T = solve_tmatrix(bse, params, [0.0, 0.0], 0)

    # Manual calculation for comparison
    χ₀ = compute_pp_susceptibility(model, params, [0.0, 0.0], 0)
    T_expected = model.U / (1.0 - model.U * χ₀)

    @test T ≈ T_expected
end

@testset "Thouless criterion - find Tc" begin
    lat = SquareLattice(8, 8)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    # For attractive Hubbard (U < 0) at half-filling, Tc should exist
    Tc = find_Tc(model, T_range=(0.01, 2.0), n_matsubara=512)

    # Tc should be positive and less than bandwidth
    @test Tc > 0.0
    @test Tc < 4.0  # bandwidth = 8t, Tc should be much smaller

    # Verify: at T just above Tc, 1 - U·χ₀ > 0 (no divergence)
    # At T just below Tc, 1 - U·χ₀ < 0 (divergence crossed)
    params_above = ThermalParameters(β=1.0/(Tc * 1.05), μ=0.0, n_matsubara=512)
    params_below = ThermalParameters(β=1.0/(Tc * 0.95), μ=0.0, n_matsubara=512)

    χ_above = compute_pp_susceptibility(model, params_above, [0.0, 0.0], 0)
    χ_below = compute_pp_susceptibility(model, params_below, [0.0, 0.0], 0)

    @test real(1.0 - model.U * χ_above) > 0
    @test real(1.0 - model.U * χ_below) < 0
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `PPLadder` not defined.

**Step 3: Write implementation**

Write `src/resummation/channels.jl`:

```julia
export PPLadder

"""Particle-particle ladder resummation channel."""
struct PPLadder
    model::HubbardModel
end
```

Write `src/resummation/bethe_salpeter.jl`:

```julia
export solve_tmatrix, find_Tc

"""
    solve_tmatrix(bse::PPLadder, params, q, m) -> ComplexF64

T-matrix from ladder resummation: T(q,iνₘ) = U / (1 - U·χ₀(q,iνₘ))
"""
function solve_tmatrix(bse::PPLadder, params::ThermalParameters,
                        q::Vector{Float64}, m::Int)::ComplexF64
    U = bse.model.U
    χ₀ = compute_pp_susceptibility(bse.model, params, q, m)
    return U / (1.0 - U * χ₀)
end

"""
    find_Tc(model; T_range, q, n_matsubara, tol) -> Float64

Find superconducting Tc via Thouless criterion: 1 - U·χ₀(q=0, iν₀=0; T) = 0.
Uses bisection method.
"""
function find_Tc(model::HubbardModel;
                  T_range::Tuple{Float64,Float64} = (0.001, 2.0),
                  q::Vector{Float64} = [0.0, 0.0],
                  n_matsubara::Int = 512,
                  tol::Float64 = 1e-4)::Float64
    T_lo, T_hi = T_range

    function criterion(T)
        params = ThermalParameters(1.0 / T, model.lattice isa SquareLattice ? 0.0 : 0.0, n_matsubara)
        χ₀ = compute_pp_susceptibility(model, params, q, 0)
        return real(1.0 - model.U * χ₀)
    end

    # Verify bracket: criterion should be negative at low T, positive at high T
    f_lo = criterion(T_lo)
    f_hi = criterion(T_hi)

    if f_lo * f_hi > 0
        error("Thouless criterion not bracketed in T_range=($T_lo, $T_hi). " *
              "f(T_lo)=$f_lo, f(T_hi)=$f_hi")
    end

    # Bisection
    while (T_hi - T_lo) > tol
        T_mid = (T_lo + T_hi) / 2.0
        f_mid = criterion(T_mid)
        if f_mid * f_lo < 0
            T_hi = T_mid
            f_hi = f_mid
        else
            T_lo = T_mid
            f_lo = f_mid
        end
    end

    return (T_lo + T_hi) / 2.0
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS. The Tc test may take a few seconds due to Matsubara sums.

**Step 5: Commit**

```bash
git add src/resummation/ test/resummation/ test/runtests.jl
git commit -m "feat: add PP ladder resummation and Thouless Tc finder"
```

---

## Task 8: Diagram Graph Representation

**Files:**
- Create: `src/diagrams/graph.jl`
- Create: `test/diagrams/test_graph.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/diagrams/test_graph.jl`:

```julia
using Test
using FeynmanEngine

@testset "FeynmanDiagram construction" begin
    v1 = Vertex(1, :i, :τ₁)
    v2 = Vertex(2, :j, :τ₂)

    p1 = Propagator(v1, v2, :↑, false)  # internal
    p2 = Propagator(v2, v1, :↓, false)

    d = FeynmanDiagram(1, [v1, v2], [p1, p2], -1, 1//1)

    @test d.order == 1
    @test length(d.vertices) == 2
    @test length(d.propagators) == 2
    @test d.sign == -1
end

@testset "Loop count" begin
    # Single bubble: 2 vertices, 2 internal propagators → 1 loop
    v1 = Vertex(1, :i, :τ₁)
    v2 = Vertex(2, :j, :τ₂)
    p1 = Propagator(v1, v2, :↑, false)
    p2 = Propagator(v2, v1, :↓, false)
    d = FeynmanDiagram(1, [v1, v2], [p1, p2], 1, 1//1)

    @test count_loops(d) == 1

    # 2nd order self-energy (sunset): 3 vertices, 2 external + 3 internal → 2 loops
    # Actually let's test with known topology:
    # Tadpole: 1 vertex, 1 self-loop → 1 loop
    v = Vertex(1, :i, :τ₁)
    p_ext_in = Propagator(Vertex(0, :ext, :τ_in), v, :↑, true)
    p_ext_out = Propagator(v, Vertex(3, :ext, :τ_out), :↑, true)
    p_loop = Propagator(v, v, :↓, false)
    d_tadpole = FeynmanDiagram(1, [v], [p_ext_in, p_ext_out, p_loop], 1, 1//1)
    @test count_loops(d_tadpole) == 1
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `Vertex` not defined.

**Step 3: Write implementation**

Write `src/diagrams/graph.jl`:

```julia
export Vertex, Propagator, FeynmanDiagram, count_loops

struct Vertex
    id::Int
    site::Symbol
    time::Symbol
end

struct Propagator
    from::Vertex
    to::Vertex
    spin::Symbol
    external::Bool
end

struct FeynmanDiagram
    order::Int
    vertices::Vector{Vertex}
    propagators::Vector{Propagator}
    sign::Int
    symmetry_factor::Rational{Int}
end

"""
Count independent loops = n_internal_edges - n_vertices + n_connected_components.
For connected diagrams: loops = n_internal_edges - n_vertices + 1.
"""
function count_loops(d::FeynmanDiagram)::Int
    internal = filter(p -> !p.external, d.propagators)
    n_edges = length(internal)
    n_verts = length(d.vertices)
    # Assuming connected diagram
    return n_edges - n_verts + 1
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/diagrams/graph.jl test/diagrams/test_graph.jl test/runtests.jl
git commit -m "feat: add Vertex, Propagator, FeynmanDiagram types with loop counting"
```

---

## Task 9: Diagram Generation from Wick Contractions

**Files:**
- Create: `src/diagrams/generate.jl`
- Create: `test/diagrams/test_generate.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/diagrams/test_generate.jl`:

```julia
using Test
using FeynmanEngine

@testset "Generate 1st order self-energy diagrams" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    diagrams = generate_diagrams(model, order=1, observable=:self_energy)

    # 1st order self-energy: Hartree (tadpole) diagram
    # c†_ext c_ext × one vertex c†↑c↑c†↓c↓
    # Should produce Hartree-type diagram(s)
    @test length(diagrams) > 0
    @test all(d -> d.order == 1, diagrams)
end

@testset "Generate 2nd order self-energy diagrams" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    diagrams = generate_diagrams(model, order=2, observable=:self_energy)

    # 2nd order: should produce multiple topologically distinct diagrams
    @test length(diagrams) > 0
    @test all(d -> d.order == 2, diagrams)
end

@testset "Generate vacuum diagrams" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    # 1st order vacuum: single closed loop
    diagrams = generate_diagrams(model, order=1, observable=:vacuum)
    @test length(diagrams) > 0
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `generate_diagrams` not defined.

**Step 3: Write implementation**

Write `src/diagrams/generate.jl`:

```julia
export generate_diagrams

"""
    generate_diagrams(model, order; observable) -> Vector{FeynmanDiagram}

Generate all Feynman diagrams at given perturbation order using Wick's theorem.

`observable` can be:
- `:vacuum` — no external legs
- `:self_energy` — 2 external legs (one creator, one annihilator)
- `:susceptibility` — 4 external legs
"""
function generate_diagrams(model::HubbardModel; order::Int, observable::Symbol=:vacuum)
    # Build operator list: external operators + vertex operators
    ops = FermionOperator[]
    ext_indices = Int[]

    if observable == :self_energy
        # External legs: c†_{ext,σ}(τ_in) ... c_{ext,σ}(τ_out)
        push!(ops, creator(:ext, :↑, :τ_in))
        push!(ext_indices, length(ops))
        push!(ops, annihilator(:ext, :↑, :τ_out))
        push!(ext_indices, length(ops))
    end

    # Internal vertices: each vertex contributes c†↑ c↑ c†↓ c↓
    for v in 1:order
        site = Symbol("v$v")
        time = Symbol("τ$v")
        push!(ops, creator(site, :↑, time))
        push!(ops, annihilator(site, :↑, time))
        push!(ops, creator(site, :↓, time))
        push!(ops, annihilator(site, :↓, time))
    end

    # Apply Wick's theorem
    terms = wick_theorem(ops)

    # Convert each WickTerm to a FeynmanDiagram
    diagrams = FeynmanDiagram[]
    for term in terms
        d = _wick_term_to_diagram(term, order, ext_indices, ops)
        if !isnothing(d)
            push!(diagrams, d)
        end
    end

    return diagrams
end

function _wick_term_to_diagram(term::WickTerm, order::Int,
                                ext_indices::Vector{Int},
                                ops::Vector{FermionOperator})
    # Build vertices from unique (site, time) pairs of internal operators
    vertex_map = Dict{Tuple{Symbol,Symbol}, Vertex}()
    vid = 0
    for c in term.contractions
        for op in [c.creator, c.annihilator]
            if op.site != :ext
                key = (op.site, op.time)
                if !haskey(vertex_map, key)
                    vid += 1
                    vertex_map[key] = Vertex(vid, op.site, op.time)
                end
            end
        end
    end

    vertices = collect(values(vertex_map))
    sort!(vertices, by=v -> v.id)

    # External pseudo-vertices
    ext_in = Vertex(0, :ext, :τ_in)
    ext_out = Vertex(-1, :ext, :τ_out)

    # Build propagators from contractions
    propagators = Propagator[]
    for c in term.contractions
        from_key = (c.creator.site, c.creator.time)
        to_key = (c.annihilator.site, c.annihilator.time)

        from_v = c.creator.site == :ext ? ext_in : vertex_map[from_key]
        to_v = c.annihilator.site == :ext ? ext_out : vertex_map[to_key]

        is_ext = (c.creator.site == :ext || c.annihilator.site == :ext)
        push!(propagators, Propagator(from_v, to_v, c.creator.spin, is_ext))
    end

    return FeynmanDiagram(order, vertices, propagators, term.sign, 1 // factorial(order))
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/diagrams/generate.jl test/diagrams/test_generate.jl test/runtests.jl
git commit -m "feat: generate Feynman diagrams from Wick contractions"
```

---

## Task 10: Diagram Classification (Topology)

**Files:**
- Create: `src/diagrams/classify.jl`
- Create: `test/diagrams/test_classify.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/diagrams/test_classify.jl`:

```julia
using Test
using FeynmanEngine

@testset "Classify diagrams by topology" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    # 2nd order self-energy: multiple Wick terms → fewer unique topologies
    diagrams = generate_diagrams(model, order=2, observable=:self_energy)
    classified = classify_diagrams(diagrams)

    # Should have fewer unique topologies than total Wick terms
    @test length(classified) <= length(diagrams)
    @test length(classified) > 0

    # Each classified entry should have a weight (multiplicity × sign)
    for (diag, weight) in classified
        @test weight != 0
        @test diag isa FeynmanDiagram
    end
end

@testset "1st order self-energy topology" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    diagrams = generate_diagrams(model, order=1, observable=:self_energy)
    classified = classify_diagrams(diagrams)

    # 1st order: should be Hartree diagram only (one topology)
    @test length(classified) == 1
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `classify_diagrams` not defined.

**Step 3: Write implementation**

Write `src/diagrams/classify.jl`:

```julia
using Graphs

export classify_diagrams

"""
    classify_diagrams(diagrams) -> Vector{Tuple{FeynmanDiagram, Int}}

Group topologically equivalent diagrams. Returns unique diagrams with weights
(sum of signs × multiplicity).
"""
function classify_diagrams(diagrams::Vector{FeynmanDiagram})
    groups = Tuple{FeynmanDiagram, Int}[]

    for d in diagrams
        found = false
        for (idx, (rep, weight)) in enumerate(groups)
            if _is_topologically_equivalent(rep, d)
                groups[idx] = (rep, weight + d.sign)
                found = true
                break
            end
        end
        if !found
            push!(groups, (d, d.sign))
        end
    end

    # Remove zero-weight diagrams (cancelled)
    return filter(((_, w),) -> w != 0, groups)
end

"""
Check topological equivalence using colored graph isomorphism.
Two diagrams are equivalent if their graph structure is the same
considering spin labels and external/internal status.
"""
function _is_topologically_equivalent(d1::FeynmanDiagram, d2::FeynmanDiagram)::Bool
    if d1.order != d2.order
        return false
    end
    if length(d1.propagators) != length(d2.propagators)
        return false
    end
    if length(d1.vertices) != length(d2.vertices)
        return false
    end

    # Build adjacency structure with spin coloring
    sig1 = _diagram_signature(d1)
    sig2 = _diagram_signature(d2)

    return sig1 == sig2
end

"""
Compute a canonical signature for a diagram topology.
Uses sorted edge descriptors relative to vertex connectivity pattern.
"""
function _diagram_signature(d::FeynmanDiagram)
    # For each vertex, compute (in_degree_↑, out_degree_↑, in_degree_↓, out_degree_↓, n_self_loops)
    # Sort these vertex descriptors to get a canonical form
    verts = d.vertices
    if isempty(verts)
        return []
    end

    vert_ids = Dict(v.id => i for (i, v) in enumerate(verts))

    # Adjacency: for each pair (i,j), count propagators by spin
    n = length(verts)
    adj = zeros(Int, n, n, 2)  # (from, to, spin_idx)

    ext_connections = Dict{Int, Vector{Tuple{Symbol, Bool}}}()  # vertex_idx → [(spin, is_incoming)]

    for p in d.propagators
        if p.external
            # Track external connections per vertex
            if haskey(vert_ids, p.from.id)
                vi = vert_ids[p.from.id]
                push!(get!(ext_connections, vi, []), (p.spin, false))
            end
            if haskey(vert_ids, p.to.id)
                vi = vert_ids[p.to.id]
                push!(get!(ext_connections, vi, []), (p.spin, true))
            end
            continue
        end

        fi = get(vert_ids, p.from.id, nothing)
        ti = get(vert_ids, p.to.id, nothing)
        if isnothing(fi) || isnothing(ti)
            continue
        end
        si = p.spin == :↑ ? 1 : 2
        adj[fi, ti, si] += 1
    end

    # Vertex descriptors
    descriptors = []
    for i in 1:n
        self_loops = adj[i, i, 1] + adj[i, i, 2]
        out_up = sum(adj[i, :, 1])
        in_up = sum(adj[:, i, 1])
        out_dn = sum(adj[i, :, 2])
        in_dn = sum(adj[:, i, 2])
        ext = sort(get(ext_connections, i, []))
        push!(descriptors, (in_up, out_up, in_dn, out_dn, self_loops, ext))
    end

    return sort(descriptors)
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS. Note: the signature-based isomorphism is approximate — for higher orders, proper graph isomorphism (via Graphs.jl) may be needed. This is sufficient for low orders.

**Step 5: Commit**

```bash
git add src/diagrams/classify.jl test/diagrams/test_classify.jl test/runtests.jl
git commit -m "feat: classify diagrams by topology using graph signatures"
```

---

## Task 11: Feynman Rules (Diagram → Symbolic Expression)

**Files:**
- Create: `src/diagrams/rules.jl`
- Create: `test/diagrams/test_rules.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/diagrams/test_rules.jl`:

```julia
using Test
using FeynmanEngine
using Symbolics

@testset "Feynman rules - 1st order Hartree" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    diagrams = generate_diagrams(model, order=1, observable=:self_energy)
    classified = classify_diagrams(diagrams)

    # Apply Feynman rules to get symbolic expression
    for (diag, weight) in classified
        expr = apply_feynman_rules(diag, model)
        @test expr !== nothing
        # Should contain U as a parameter
        @test expr isa FeynmanExpression
    end
end

@testset "FeynmanExpression" begin
    # A FeynmanExpression stores the integrand, loop variables, and prefactor
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    diagrams = generate_diagrams(model, order=1, observable=:self_energy)
    classified = classify_diagrams(diagrams)

    diag, weight = classified[1]
    expr = apply_feynman_rules(diag, model)

    @test expr.n_loops == count_loops(diag)
    @test expr.prefactor != 0
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `apply_feynman_rules` not defined.

**Step 3: Write implementation**

Write `src/diagrams/rules.jl`:

```julia
export FeynmanExpression, apply_feynman_rules

"""
Symbolic representation of a Feynman diagram's mathematical expression.
"""
struct FeynmanExpression
    n_loops::Int
    prefactor::Float64           # (-U)^n / n! × sign × symmetry_factor
    propagators::Vector{Tuple{Symbol, Symbol}}  # (momentum_var, frequency_var) per internal line
    description::String          # human-readable description
end

"""
    apply_feynman_rules(d::FeynmanDiagram, model::HubbardModel) -> FeynmanExpression

Convert a Feynman diagram to its mathematical expression via Feynman rules:
1. Each internal propagator → G₀(k, iωₙ)
2. Each vertex → (-U)
3. Each loop → Σₖ Σ_ωₙ (1/Nβ)
4. Fermion sign and symmetry factor
"""
function apply_feynman_rules(d::FeynmanDiagram, model::HubbardModel)
    n_loops = count_loops(d)

    # Prefactor: (-U)^order × sign × symmetry_factor
    prefactor = (-model.U)^d.order * d.sign * Float64(d.symmetry_factor)

    # Label internal propagators with loop variables
    internal_props = filter(p -> !p.external, d.propagators)
    prop_vars = Tuple{Symbol, Symbol}[]
    for (i, p) in enumerate(internal_props)
        k_var = Symbol("k$i")
        ω_var = Symbol("ω$i")
        push!(prop_vars, (k_var, ω_var))
    end

    desc = "Order $(d.order): $(n_loops) loop(s), $(length(internal_props)) internal propagators"

    return FeynmanExpression(n_loops, prefactor, prop_vars, desc)
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/diagrams/rules.jl test/diagrams/test_rules.jl test/runtests.jl
git commit -m "feat: apply Feynman rules to convert diagrams to expressions"
```

---

## Task 12: Diagram Visualization

**Files:**
- Create: `src/visualization/diagram_plot.jl`
- Create: `test/visualization/test_diagram_plot.jl`
- Modify: `test/runtests.jl`

**Step 1: Write the failing test**

Write `test/visualization/test_diagram_plot.jl`:

```julia
using Test
using FeynmanEngine

@testset "Diagram plot - smoke test" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    diagrams = generate_diagrams(model, order=1, observable=:self_energy)
    classified = classify_diagrams(diagrams)

    diag, weight = classified[1]

    # Should return a Figure without error
    fig = plot_diagram(diag)
    @test fig !== nothing

    # Save to file
    save_path = tempname() * ".png"
    save_diagram(fig, save_path)
    @test isfile(save_path)
    rm(save_path)
end

@testset "Plot all diagrams" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    diagrams = generate_diagrams(model, order=1, observable=:self_energy)
    classified = classify_diagrams(diagrams)

    fig = plot_all_diagrams(classified)
    @test fig !== nothing
end
```

**Step 2: Run test to verify it fails**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL — `plot_diagram` not defined.

**Step 3: Write implementation**

Write `src/visualization/diagram_plot.jl`:

```julia
using CairoMakie

export plot_diagram, save_diagram, plot_all_diagrams

"""
    plot_diagram(d::FeynmanDiagram) -> Figure

Draw a Feynman diagram using Makie.
Vertices are circles, propagators are arrows, self-loops are drawn as loops.
"""
function plot_diagram(d::FeynmanDiagram)
    fig = Figure(size=(400, 300))
    ax = Axis(fig[1, 1], title="Order $(d.order), sign=$(d.sign)",
              aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)

    # Layout: place vertices on a horizontal line
    n_verts = length(d.vertices)
    positions = Dict{Int, Tuple{Float64, Float64}}()

    for (i, v) in enumerate(d.vertices)
        x = Float64(i)
        y = 0.0
        positions[v.id] = (x, y)
    end

    # External pseudo-vertices
    if n_verts > 0
        positions[0] = (0.0, 0.0)       # ext_in
        positions[-1] = (Float64(n_verts + 1), 0.0)  # ext_out
    end

    # Draw propagators
    for p in d.propagators
        from_pos = get(positions, p.from.id, (0.0, 0.0))
        to_pos = get(positions, p.to.id, (0.0, 0.0))

        color = p.spin == :↑ ? :blue : :red
        style = p.external ? :dash : :solid

        if p.from.id == p.to.id
            # Self-loop: draw as arc above vertex
            cx, cy = from_pos
            θ = range(0, 2π, length=50)
            r = 0.3
            loop_x = cx .+ r .* cos.(θ)
            loop_y = cy .+ r .+ r .* sin.(θ)
            lines!(ax, loop_x, loop_y, color=color, linewidth=2)
        else
            arrows!(ax, [from_pos[1]], [from_pos[2]],
                    [to_pos[1] - from_pos[1]], [to_pos[2] - from_pos[2]],
                    color=color, linewidth=2, arrowsize=10)
        end
    end

    # Draw vertices
    for (i, v) in enumerate(d.vertices)
        x, y = positions[v.id]
        scatter!(ax, [x], [y], color=:black, markersize=15)
        text!(ax, x, y + 0.2, text=string(v.site), align=(:center, :bottom), fontsize=12)
    end

    return fig
end

"""Save a figure to file."""
function save_diagram(fig, path::String)
    CairoMakie.save(path, fig)
end

"""
    plot_all_diagrams(classified) -> Figure

Plot all classified diagrams in a grid layout.
"""
function plot_all_diagrams(classified::Vector{Tuple{FeynmanDiagram, Int}})
    n = length(classified)
    cols = min(n, 3)
    rows = ceil(Int, n / cols)

    fig = Figure(size=(400 * cols, 300 * rows))

    for (idx, (diag, weight)) in enumerate(classified)
        row = div(idx - 1, cols) + 1
        col = mod(idx - 1, cols) + 1

        ax = Axis(fig[row, col],
                  title="Topology $idx (weight=$weight)",
                  aspect=DataAspect())
        hidedecorations!(ax)
        hidespines!(ax)

        # Simplified drawing
        n_verts = length(diag.vertices)
        positions = Dict{Int, Tuple{Float64, Float64}}()
        for (i, v) in enumerate(diag.vertices)
            positions[v.id] = (Float64(i), 0.0)
        end
        positions[0] = (0.0, 0.0)
        positions[-1] = (Float64(n_verts + 1), 0.0)

        for p in diag.propagators
            from_pos = get(positions, p.from.id, (0.0, 0.0))
            to_pos = get(positions, p.to.id, (0.0, 0.0))
            color = p.spin == :↑ ? :blue : :red

            if p.from.id == p.to.id
                cx, cy = from_pos
                θ = range(0, 2π, length=50)
                r = 0.3
                lines!(ax, cx .+ r .* cos.(θ), cy .+ r .+ r .* sin.(θ),
                       color=color, linewidth=2)
            else
                arrows!(ax, [from_pos[1]], [from_pos[2]],
                        [to_pos[1] - from_pos[1]], [to_pos[2] - from_pos[2]],
                        color=color, linewidth=2, arrowsize=10)
            end
        end

        for v in diag.vertices
            x, y = positions[v.id]
            scatter!(ax, [x], [y], color=:black, markersize=15)
        end
    end

    return fig
end
```

**Step 4: Run test to verify it passes**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

**Step 5: Commit**

```bash
git add src/visualization/diagram_plot.jl test/visualization/test_diagram_plot.jl test/runtests.jl
git commit -m "feat: add Feynman diagram visualization with Makie"
```

---

## Task 13: End-to-End Example — Ladder Resummation

**Files:**
- Create: `examples/ladder_resummation.jl`
- Create: `examples/second_order_self_energy.jl`

**Step 1: Write ladder resummation example**

Write `examples/ladder_resummation.jl`:

```julia
using FeynmanEngine
using CairoMakie

# 2D square lattice attractive Hubbard model
println("=== Attractive Hubbard Model: Ladder Resummation ===")
println()

lattice = SquareLattice(16, 16)
model = HubbardModel(lattice, t=1.0, U=-2.0)

# Find Tc via Thouless criterion
println("Searching for Tc...")
Tc = find_Tc(model, T_range=(0.01, 2.0), n_matsubara=512)
println("Tc = $Tc")
println()

# Plot 1 - U·χ₀ vs temperature
T_values = range(Tc * 0.5, 2.0, length=100)
criterion_values = Float64[]

for T in T_values
    params = ThermalParameters(1.0 / T, 0.0, 512)
    χ₀ = compute_pp_susceptibility(model, params, [0.0, 0.0], 0)
    push!(criterion_values, real(1.0 - model.U * χ₀))
end

fig = Figure(size=(600, 400))
ax = Axis(fig[1, 1], xlabel="T/t", ylabel="1 - U·χ₀(0,0)",
          title="Thouless Criterion (U=$(model.U))")
lines!(ax, collect(T_values), criterion_values)
hlines!(ax, [0.0], color=:red, linestyle=:dash)
vlines!(ax, [Tc], color=:green, linestyle=:dash, label="Tc=$(@sprintf("%.4f", Tc))")
axislegend(ax)

save("thouless_criterion.png", fig)
println("Saved: thouless_criterion.png")
```

**Step 2: Write self-energy example**

Write `examples/second_order_self_energy.jl`:

```julia
using FeynmanEngine

println("=== 2nd Order Self-Energy Diagrams ===")
println()

lattice = SquareLattice(4, 4)
model = HubbardModel(lattice, t=1.0, U=-2.0)

# Generate and classify diagrams
for order in 1:2
    diagrams = generate_diagrams(model, order=order, observable=:self_energy)
    classified = classify_diagrams(diagrams)

    println("Order $order:")
    println("  Total Wick terms: $(length(diagrams))")
    println("  Unique topologies: $(length(classified))")
    for (i, (diag, weight)) in enumerate(classified)
        expr = apply_feynman_rules(diag, model)
        println("  Topology $i: weight=$weight, loops=$(expr.n_loops)")
        println("    $(expr.description)")
    end
    println()
end

# Visualize
diagrams = generate_diagrams(model, order=2, observable=:self_energy)
classified = classify_diagrams(diagrams)
fig = plot_all_diagrams(classified)
save_diagram(fig, "self_energy_2nd_order.png")
println("Saved: self_energy_2nd_order.png")
```

**Step 3: Run examples to verify they work**

```bash
julia --project=. examples/second_order_self_energy.jl
julia --project=. examples/ladder_resummation.jl
```

Expected: Both run successfully, print results, and save PNG files.

**Step 4: Commit**

```bash
git add examples/
git commit -m "feat: add ladder resummation and self-energy example scripts"
```

---

## Task 14: Final Integration Test

**Files:**
- Create: `test/test_integration.jl`
- Modify: `test/runtests.jl`

**Step 1: Write integration test**

Write `test/test_integration.jl`:

```julia
using Test
using FeynmanEngine

@testset "End-to-end: Wick → Diagrams → Classification" begin
    lat = SquareLattice(4, 4)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    # Generate, classify, apply rules
    diagrams = generate_diagrams(model, order=2, observable=:self_energy)
    classified = classify_diagrams(diagrams)
    expressions = [apply_feynman_rules(d, model) for (d, _) in classified]

    @test length(classified) > 0
    @test all(e -> e.n_loops > 0, expressions)
end

@testset "End-to-end: Ladder resummation Tc" begin
    lat = SquareLattice(8, 8)
    model = HubbardModel(lat, t=1.0, U=-2.0)

    # Should find Tc without error
    Tc = find_Tc(model, T_range=(0.01, 2.0), n_matsubara=256)
    @test 0.0 < Tc < 4.0

    # Verify T-matrix diverges near Tc
    bse = PPLadder(model)
    params_at_Tc = ThermalParameters(1.0 / (Tc * 1.001), 0.0, 256)
    T_val = solve_tmatrix(bse, params_at_Tc, [0.0, 0.0], 0)
    @test abs(T_val) > abs(model.U) * 10  # should be large near Tc
end

@testset "End-to-end: Repulsive Hubbard has no pp instability" begin
    lat = SquareLattice(8, 8)
    model = HubbardModel(lat, t=1.0, U=2.0)  # repulsive

    # For repulsive U > 0, pp ladder should NOT diverge
    # 1 - U·χ₀ should remain > 0 for all reasonable T
    params = ThermalParameters(β=20.0, μ=0.0, n_matsubara=512)
    χ₀ = compute_pp_susceptibility(model, params, [0.0, 0.0], 0)
    @test real(1.0 - model.U * χ₀) > 0
end
```

Add to `test/runtests.jl`:

```julia
include("test_integration.jl")
```

**Step 2: Run all tests**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: ALL PASS.

**Step 3: Commit**

```bash
git add test/test_integration.jl test/runtests.jl
git commit -m "feat: add end-to-end integration tests"
```

---

## Summary of Tasks

| Task | Description | Key Deliverable |
|------|-------------|-----------------|
| 1 | Project scaffolding | Package structure, dependencies |
| 2 | FermionOperator types | Core data types |
| 3 | Wick's theorem | Automatic contraction enumeration |
| 4 | Lattice + Hubbard model | Model definitions |
| 5 | Free Green's function | G₀(k, iωₙ) evaluation |
| 6 | PP susceptibility | χ₀(q, iνₘ) computation |
| 7 | Bethe-Salpeter + Tc | Ladder resummation, Thouless criterion |
| 8 | Diagram graph types | Vertex, Propagator, FeynmanDiagram |
| 9 | Diagram generation | Wick terms → diagram objects |
| 10 | Diagram classification | Topology grouping |
| 11 | Feynman rules | Diagram → symbolic expression |
| 12 | Visualization | Diagram drawing with Makie |
| 13 | Examples | Ladder resummation + self-energy scripts |
| 14 | Integration tests | End-to-end verification |

**Dependencies:** Tasks 1→2→3, Tasks 1→4→5→6→7, Tasks 8→9→10→11→12. Tasks 3 and 7 can proceed in parallel. Task 13-14 depend on all previous tasks.
