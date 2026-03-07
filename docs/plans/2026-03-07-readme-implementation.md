# README & Usage Documentation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create an English README.md and a Japanese usage guide (docs/usage-ja.md) for the feynman-engine project.

**Architecture:** Two standalone Markdown files. README.md at the project root for GitHub display; docs/usage-ja.md as a Japanese companion guide with additional physical explanations.

**Tech Stack:** Markdown only. No build step required.

---

### Task 1: Create README.md

**Files:**
- Create: `README.md`

**Step 1: Write README.md**

Create `README.md` with the following sections:

1. **Title**: `# feynman-engine` with subtitle "A perturbation theory engine for the Hubbard model in Rust"
2. **Overview**: One paragraph describing the pipeline (Wick's theorem → Feynman diagrams → symbolic expressions → ladder resummation → Tc via Thouless criterion). Mention finite-temperature Matsubara formalism.
3. **Features**: Bullet list:
   - Automatic Wick contraction with spin conservation and fermion sign
   - Feynman diagram generation and topological classification
   - Symbolic expression tree with Feynman rules
   - Particle-particle ladder resummation (Bethe-Salpeter)
   - Thouless criterion Tc finder (bisection)
   - Graphviz DOT visualization (spin-colored propagators)
4. **Installation**:
   ```
   git clone <repo-url>
   cd feynman
   cargo build
   cargo test
   ```
5. **Quick Start**: Two subsections with `cargo run --example` commands and output snippets:
   - **Ladder resummation** — show the Tc output and Thouless table (trimmed to first few lines + last line)
   - **Self-energy diagrams** — show 1st order output (trim 2nd order to summary)
6. **Architecture**: The data flow diagram from CLAUDE.md (both the algebraic and numerical pipelines) and a brief module list table.
7. **Physics Background**: 3-4 paragraphs covering:
   - Hubbard model and interaction vertex
   - Matsubara formalism (finite temperature, fermionic frequencies)
   - PP susceptibility and Thouless criterion
   - Attractive vs repulsive Hubbard
8. **License**: "TBD" or omit

**Step 2: Review README.md for accuracy**

Read through the file and verify all commands, output snippets, and physics descriptions are correct.

**Step 3: Commit**

```bash
git add README.md
git commit -m "docs: add README.md"
```

---

### Task 2: Create docs/usage-ja.md

**Files:**
- Create: `docs/usage-ja.md`

**Step 1: Write docs/usage-ja.md**

Create `docs/usage-ja.md` with the following sections (all in Japanese):

1. **タイトル**: `# feynman-engine 使い方ガイド`
2. **概要**: README.mdの Overview を日本語化
3. **機能一覧**: README.md の Features を日本語化
4. **インストール**: README.md と同じコマンド、日本語の説明
5. **使い方**:
   - **梯子型再総和 (ladder resummation)**: exampleの実行コマンド、全出力、物理的な意味の解説（Thouless判定条件が0になる温度がTcであること、引力Hubbardでpp channelにCooper不安定性が生じること）
   - **自己エネルギーダイアグラム**: exampleの実行コマンド、出力、各トポロジーの物理的意味（Hartree項、Fock項、2次の寄与）、DOT出力のGraphvizでの可視化方法
6. **アーキテクチャ**: README.mdの Architecture セクションを日本語化
7. **物理的背景**: README.md の Physics Background をより詳しく日本語で。特に：
   - Hubbard模型の定義式 H = -t Σ c†c + U Σ n↑n↓
   - 松原形式の説明（虚時間、フェルミオン松原周波数 iω_n = i(2n+1)πT）
   - PP感受率 χ₀ の定義と物理的意味
   - Thouless判定条件 1 - Uχ₀(T) = 0 の意味
8. **モジュール解説**: 各モジュール（algebra, models, diagrams, symbolic, resummation, numerical, visualization）の簡潔な説明と主要な型・関数

**Step 2: Review for accuracy**

Read through the file and verify all Japanese text, commands, and physics descriptions are correct.

**Step 3: Commit**

```bash
git add docs/usage-ja.md
git commit -m "docs: add Japanese usage guide"
```

---

### Task 3: Cross-link the two documents

**Files:**
- Modify: `README.md` (add link to Japanese guide)
- Modify: `docs/usage-ja.md` (add link back to README)

**Step 1: Add links**

- In README.md, add near the top: `[日本語ガイド](docs/usage-ja.md)`
- In docs/usage-ja.md, add near the top: `[English README](../README.md)`

**Step 2: Commit**

```bash
git add README.md docs/usage-ja.md
git commit -m "docs: cross-link README and Japanese guide"
```
