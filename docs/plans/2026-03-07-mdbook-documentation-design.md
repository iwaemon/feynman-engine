# mdBook Documentation Design

## Overview

mdBookを使ってFeynmanプロジェクトのドキュメントを生成する。対象読者はライブラリを使う開発者と物理を学ぶ学生/研究者。英語メイン、日本語版は段階的に翻訳。

## Prerequisites

- feat/d3-feynman-visualization ブランチを master にマージしてからドキュメント作業を開始

## Directory Structure

```
book/
├── book.toml
├── src/
│   ├── SUMMARY.md
│   ├── introduction.md
│   ├── theory/
│   │   ├── hubbard-model.md
│   │   ├── perturbation-theory.md
│   │   ├── feynman-diagrams.md
│   │   └── ladder-resummation.md
│   ├── getting-started/
│   │   ├── installation.md
│   │   └── quick-start.md
│   ├── guide/
│   │   ├── wick-theorem.md
│   │   ├── diagram-generation.md
│   │   ├── symbolic-expressions.md
│   │   ├── numerical.md
│   │   └── finding-tc.md
│   └── visualization/
│       ├── graphviz.md
│       └── d3-interactive.md
book-ja/
├── book.toml
├── src/
│   ├── SUMMARY.md
│   ├── introduction.md
│   └── ...  (placeholder pages)
```

## Book Configuration

### English (book/book.toml)

```toml
[book]
title = "Feynman -- Perturbation Theory Engine"
authors = ["shumpei"]
language = "en"

[build]
build-dir = "../target/book"

[output.html]
mathjax-support = true
```

### Japanese (book-ja/book.toml)

```toml
[book]
title = "Feynman -- Perturbation Theory Engine"
authors = ["shumpei"]
language = "ja"

[build]
build-dir = "../target/book-ja"

[output.html]
mathjax-support = true
```

## Content Scope

### Theory (理論編)

Structure: 物理の背景 -> コードへの対応

- **hubbard-model.md** -- Hamiltonian, square lattice dispersion e(k)=-2t(cos kx+cos ky), attractive vs repulsive physics
- **perturbation-theory.md** -- Finite-temperature Matsubara formalism, Wick's theorem, fermion sign via inversion counting, spin conservation filtering
- **feynman-diagrams.md** -- Diagram generation (Observable enum), topology classification (vertex-descriptor signatures), Feynman rules application
- **ladder-resummation.md** -- PP bubble chi_0, T-matrix U/(1-U*chi_0), Thouless criterion, bisection for Tc

Each chapter includes LaTeX math (via MathJax) and links to corresponding Rust modules.

### Getting Started (入門)

- **installation.md** -- `cargo build`, Rust toolchain requirements
- **quick-start.md** -- Run `ladder_resummation` example end-to-end to find Tc

### Guide (実践ガイド)

Module-by-module usage with code examples showing how theory maps to code:

- **wick-theorem.md** -- algebra module (FermionOperator, WickTerm, wick_theorem())
- **diagram-generation.md** -- diagrams module (FeynmanDiagram, generate_diagrams(), classify_diagrams())
- **symbolic-expressions.md** -- symbolic module (Expr tree, apply_feynman_rules())
- **numerical.md** -- numerical module (evaluate_g0, compute_pp_susceptibility, Matsubara frequencies)
- **finding-tc.md** -- resummation module (PPLadder, solve_tmatrix, find_tc with bisection)

### Visualization (可視化)

- **graphviz.md** -- to_dot()/to_dot_all() usage, DOT to PNG conversion
- **d3-interactive.md** -- D3.js templates, visualize_d3.sh usage, screenshots

## Japanese Version Strategy

- Initial release: SUMMARY.md + introduction.md translated
- Remaining pages: placeholder with "This page has not been translated yet. See the English version."
- Translation added incrementally in future tasks

## Build Commands

```bash
# English
mdbook build book/

# Japanese
mdbook build book-ja/

# Local preview
mdbook serve book/
```

## Non-Goals

- No GitHub Actions / CI deployment (manual build only)
- No changes to Cargo.toml (mdBook is independent of the Rust project)
- No `cargo doc` integration (separate concern)
