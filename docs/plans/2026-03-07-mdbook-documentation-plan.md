# mdBook Documentation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Set up mdBook documentation with theory-first structure, English main + Japanese scaffold, MathJax support, and D3.js visualization chapter.

**Architecture:** Two independent mdBook projects (`book/` for English, `book-ja/` for Japanese) with shared structure. English content written in full; Japanese has translated intro + placeholders. MathJax enabled for LaTeX math rendering.

**Tech Stack:** mdBook, MathJax, Markdown

**Prerequisites:** Merge `feat/d3-feynman-visualization` into `master` before starting.

---

### Task 0: Merge feature branch to master

**Step 1: Merge**

```bash
git checkout master
git merge feat/d3-feynman-visualization
```

**Step 2: Verify**

```bash
cargo test
cargo build
```

Expected: All tests pass, build succeeds.

**Step 3: Commit note**

No commit needed (merge commit is automatic).

---

### Task 1: Initialize mdBook structure (English)

**Files:**
- Create: `book/book.toml`
- Create: `book/src/SUMMARY.md`

**Step 1: Install mdbook if needed**

```bash
cargo install mdbook
```

**Step 2: Create `book/book.toml`**

```toml
[book]
title = "Feynman - Perturbation Theory Engine"
authors = ["shumpei"]
language = "en"

[build]
build-dir = "../target/book"

[output.html]
mathjax-support = true
```

**Step 3: Create `book/src/SUMMARY.md`**

```markdown
# Summary

[Introduction](./introduction.md)

# Theory

- [The Hubbard Model](./theory/hubbard-model.md)
- [Perturbation Theory & Wick's Theorem](./theory/perturbation-theory.md)
- [Feynman Diagrams](./theory/feynman-diagrams.md)
- [Ladder Resummation & Tc](./theory/ladder-resummation.md)

# Getting Started

- [Installation](./getting-started/installation.md)
- [Quick Start](./getting-started/quick-start.md)

# Guide

- [Wick's Theorem](./guide/wick-theorem.md)
- [Diagram Generation](./guide/diagram-generation.md)
- [Symbolic Expressions](./guide/symbolic-expressions.md)
- [Numerical Evaluation](./guide/numerical.md)
- [Finding Tc](./guide/finding-tc.md)

# Visualization

- [Graphviz Export](./visualization/graphviz.md)
- [D3.js Interactive Diagrams](./visualization/d3-interactive.md)
```

**Step 4: Verify build**

```bash
mdbook build book/
```

Expected: Builds with warnings about missing files (expected at this stage).

**Step 5: Add `/target/book` and `/target/book-ja` to `.gitignore`**

Append:
```
/target/book
/target/book-ja
```

**Step 6: Commit**

```bash
git add book/book.toml book/src/SUMMARY.md .gitignore
git commit -m "feat: initialize mdBook structure with SUMMARY and config"
```

---

### Task 2: Write introduction.md

**Files:**
- Create: `book/src/introduction.md`

**Step 1: Write `book/src/introduction.md`**

Content: project overview adapted from README.md. Cover what feynman-engine does, two pipelines (symbolic + numerical), intended audience. Keep concise (200-300 words). Include the module map table from README.

**Step 2: Verify**

```bash
mdbook build book/
```

Expected: Builds successfully, introduction page renders.

**Step 3: Commit**

```bash
git add book/src/introduction.md
git commit -m "docs: add introduction page"
```

---

### Task 3: Write theory chapters

**Files:**
- Create: `book/src/theory/hubbard-model.md`
- Create: `book/src/theory/perturbation-theory.md`
- Create: `book/src/theory/feynman-diagrams.md`
- Create: `book/src/theory/ladder-resummation.md`

**Step 1: Write `book/src/theory/hubbard-model.md`**

Content:
- Hamiltonian: \\( H = -t \sum_{\langle i,j \rangle, \sigma} c^\dagger_{i\sigma} c_{j\sigma} + U \sum_i n_{i\uparrow} n_{i\downarrow} \\)
- Square lattice dispersion: \\( \varepsilon(\mathbf{k}) = -2t(\cos k_x + \cos k_y) \\)
- Attractive (U<0) vs repulsive (U>0) physics
- Reference to `models/` module (`SquareLattice`, `HubbardModel`)

**Step 2: Write `book/src/theory/perturbation-theory.md`**

Content:
- Matsubara formalism: fermionic frequencies \\( i\omega_n = i(2n+1)\pi/\beta \\)
- Bare Green's function: \\( G_0(k, i\omega_n) = 1/(i\omega_n - \varepsilon_k + \mu) \\)
- Wick's theorem: all possible contractions with fermion sign
- Spin conservation filtering
- Reference to `algebra/` module (`FermionOperator`, `WickTerm`, `wick_theorem()`)

**Step 3: Write `book/src/theory/feynman-diagrams.md`**

Content:
- Diagram generation from Wick contractions
- Graph representation (vertices, propagators)
- Topological classification via vertex-descriptor signatures
- Feynman rules: converting diagrams to symbolic expressions
- Reference to `diagrams/` and `symbolic/` modules

**Step 4: Write `book/src/theory/ladder-resummation.md`**

Content:
- PP bubble: \\( \chi_0(q, i\nu_m) = -\frac{1}{N\beta} \sum_k \sum_n G_0(k, i\omega_n) G_0(q-k, i\nu_m - i\omega_n) \\)
- T-matrix: \\( T(q, i\nu_m) = U / (1 - U\chi_0) \\)
- Thouless criterion: \\( 1 - U\chi_0(q=0, i\nu_0=0) = 0 \\)
- Bisection algorithm for Tc
- Reference to `resummation/` and `numerical/` modules

**Step 5: Verify**

```bash
mdbook build book/
```

Expected: All four theory pages render with MathJax equations.

**Step 6: Commit**

```bash
git add book/src/theory/
git commit -m "docs: add theory chapters (Hubbard, perturbation, diagrams, resummation)"
```

---

### Task 4: Write getting-started chapters

**Files:**
- Create: `book/src/getting-started/installation.md`
- Create: `book/src/getting-started/quick-start.md`

**Step 1: Write `book/src/getting-started/installation.md`**

Content:
- Rust toolchain requirement (edition 2021, MSRV 1.70)
- `git clone`, `cargo build`, `cargo test`
- Optional: Graphviz for DOT rendering, Chrome for PNG export

**Step 2: Write `book/src/getting-started/quick-start.md`**

Content:
- Run `cargo run --example ladder_resummation` with expected output
- Run `cargo run --example second_order_self_energy` with expected output
- Brief explanation of what each example does

**Step 3: Verify**

```bash
mdbook build book/
```

**Step 4: Commit**

```bash
git add book/src/getting-started/
git commit -m "docs: add getting-started chapters (installation, quick-start)"
```

---

### Task 5: Write guide chapters

**Files:**
- Create: `book/src/guide/wick-theorem.md`
- Create: `book/src/guide/diagram-generation.md`
- Create: `book/src/guide/symbolic-expressions.md`
- Create: `book/src/guide/numerical.md`
- Create: `book/src/guide/finding-tc.md`

**Step 1: Write `book/src/guide/wick-theorem.md`**

Content: How to use the `algebra` module.
- Create `FermionOperator` instances with `Spin::Up`/`Spin::Down`
- Call `wick_theorem()` on a slice of operators
- Interpret `WickTerm` results (contractions, sign)
- Code example from `src/lib.rs` or tests

**Step 2: Write `book/src/guide/diagram-generation.md`**

Content: How to use the `diagrams` module.
- `generate_diagrams(&model, order, Observable::SelfEnergy)`
- `classify_diagrams()` for topological grouping
- Interpreting `FeynmanDiagram` graph (vertices, propagators)
- Code example from `second_order_self_energy.rs`

**Step 3: Write `book/src/guide/symbolic-expressions.md`**

Content: How to use `apply_feynman_rules()`.
- `FeynmanExpression` fields: `expr`, `description`, `prefactor`
- `Expr` enum variants: `G0`, `U`, `Sum`, `Mul`, `Add`, `Inv`, `Scalar`, `Neg`
- Display formatting
- Code example

**Step 4: Write `book/src/guide/numerical.md`**

Content: Numerical evaluation.
- `evaluate_g0()` for bare Green's function
- `ThermalParams::new(beta, mu, n_matsubara)`
- `compute_pp_susceptibility()` for the PP bubble
- Matsubara frequency index mapping: `m-1-n` for \\( i\nu_m - i\omega_n \\)

**Step 5: Write `book/src/guide/finding-tc.md`**

Content: Full workflow to find Tc.
- `solve_tmatrix()` for T-matrix at given temperature
- `find_tc(&model, (t_low, t_high), n_matsubara, tolerance)` bisection
- Interpreting results
- Code from `ladder_resummation.rs` example

**Step 6: Verify**

```bash
mdbook build book/
```

**Step 7: Commit**

```bash
git add book/src/guide/
git commit -m "docs: add guide chapters (wick, diagrams, symbolic, numerical, tc)"
```

---

### Task 6: Write visualization chapters

**Files:**
- Create: `book/src/visualization/graphviz.md`
- Create: `book/src/visualization/d3-interactive.md`

**Step 1: Write `book/src/visualization/graphviz.md`**

Content:
- `to_dot()` / `to_dot_all()` usage
- CLI: `--dot diagrams.dot` flag
- Rendering: `dot -Tpng diagrams.dot -o diagrams.png`
- Color scheme: spin-up=blue, spin-down=red, external=dashed

**Step 2: Write `book/src/visualization/d3-interactive.md`**

Content:
- JSON export: `--json diagrams.json` flag and `to_json_all()`
- `visualize_d3.sh` script usage
- Two templates: `template_dot_vertex.html` (dot-vertex style), `template_wavy.html` (wavy interaction style)
- SVG export buttons in templates
- PNG export via headless Chrome
- Screenshots of the rendered diagrams (optional: add images to `book/src/images/`)

**Step 3: Verify**

```bash
mdbook build book/
```

**Step 4: Commit**

```bash
git add book/src/visualization/
git commit -m "docs: add visualization chapters (graphviz, d3-interactive)"
```

---

### Task 7: Initialize Japanese version

**Files:**
- Create: `book-ja/book.toml`
- Create: `book-ja/src/SUMMARY.md`
- Create: `book-ja/src/introduction.md`
- Create placeholder files for all other chapters

**Step 1: Create `book-ja/book.toml`**

```toml
[book]
title = "Feynman - 摂動論エンジン"
authors = ["shumpei"]
language = "ja"

[build]
build-dir = "../target/book-ja"

[output.html]
mathjax-support = true
```

**Step 2: Create `book-ja/src/SUMMARY.md`**

Same structure as English but with Japanese chapter titles:

```markdown
# 目次

[はじめに](./introduction.md)

# 理論

- [Hubbardモデル](./theory/hubbard-model.md)
- [摂動論とWickの定理](./theory/perturbation-theory.md)
- [ファインマンダイアグラム](./theory/feynman-diagrams.md)
- [ラダー再足し上げとTc](./theory/ladder-resummation.md)

# 導入

- [インストール](./getting-started/installation.md)
- [クイックスタート](./getting-started/quick-start.md)

# ガイド

- [Wickの定理](./guide/wick-theorem.md)
- [ダイアグラム生成](./guide/diagram-generation.md)
- [記号式](./guide/symbolic-expressions.md)
- [数値計算](./guide/numerical.md)
- [Tcの探索](./guide/finding-tc.md)

# 可視化

- [Graphviz出力](./visualization/graphviz.md)
- [D3.jsインタラクティブダイアグラム](./visualization/d3-interactive.md)
```

**Step 3: Translate `book-ja/src/introduction.md`**

Translate the English introduction into Japanese.

**Step 4: Create placeholder files for all other chapters**

Each placeholder contains:

```markdown
> このページはまだ翻訳されていません。[英語版](../../book/src/<path>.md)を参照してください。
>
> This page has not been translated yet. Please see the English version.
```

Create placeholders for:
- `book-ja/src/theory/hubbard-model.md`
- `book-ja/src/theory/perturbation-theory.md`
- `book-ja/src/theory/feynman-diagrams.md`
- `book-ja/src/theory/ladder-resummation.md`
- `book-ja/src/getting-started/installation.md`
- `book-ja/src/getting-started/quick-start.md`
- `book-ja/src/guide/wick-theorem.md`
- `book-ja/src/guide/diagram-generation.md`
- `book-ja/src/guide/symbolic-expressions.md`
- `book-ja/src/guide/numerical.md`
- `book-ja/src/guide/finding-tc.md`
- `book-ja/src/visualization/graphviz.md`
- `book-ja/src/visualization/d3-interactive.md`

**Step 5: Verify**

```bash
mdbook build book-ja/
```

Expected: Builds successfully with Japanese intro and placeholder pages.

**Step 6: Commit**

```bash
git add book-ja/
git commit -m "docs: add Japanese mdBook scaffold with translated intro"
```

---

### Task 8: Final verification and cleanup

**Step 1: Build both books**

```bash
mdbook build book/
mdbook build book-ja/
```

Expected: Both build without errors.

**Step 2: Preview and check**

```bash
mdbook serve book/
```

Open in browser, verify:
- MathJax renders equations correctly
- All chapters are accessible from sidebar
- Code examples are syntax-highlighted

**Step 3: Commit any final fixes**

```bash
git add -A
git commit -m "docs: finalize mdBook documentation"
```
