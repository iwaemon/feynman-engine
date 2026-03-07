# feynman-engine 使い方ガイド

[English README](../README.md)

## 概要

`feynman-engine` は有限温度における引力 Hubbard 模型の摂動展開を自動化する Rust ライブラリです。フェルミオン演算子から出発し、Wick の定理を適用してフェルミオン符号とスピン保存を満たすすべての縮約を生成し、Feynman 図を構成してトポロジーで分類し、Feynman 規則により記号的な式へ変換します。数値計算側では、松原形式で粒子-粒子感受率を評価し、梯子型再総和（Bethe-Salpeter 方程式）により Thouless 判定条件を通じて超伝導臨界温度 Tc を求めます。

## 機能一覧

- スピン保存とフェルミオン符号を考慮した自動 Wick 縮約
- Feynman 図の生成とトポロジー分類
- Feynman 規則に基づく記号式ツリー
- 粒子-粒子梯子型再総和（Bethe-Salpeter）
- Thouless 判定条件に基づく Tc 探索（二分法）
- Graphviz DOT による可視化（スピンで色分けした伝播線）

## インストール

### 前提条件

- [Rust](https://www.rust-lang.org/tools/install) ツールチェイン（edition 2021、MSRV 1.70）
- [Graphviz](https://graphviz.org/)（任意、DOT 図の PNG/SVG 描画に必要）
  - macOS: `brew install graphviz`
  - Ubuntu: `apt install graphviz`

```bash
git clone <repo-url>
cd feynman
cargo build
cargo test
```

リポジトリの URL を指定してクローンし、プロジェクト直下で `cargo build` と `cargo test` を実行してください。

## 使い方

### 梯子型再総和 (ladder resummation)

正方格子上の引力 Hubbard 模型で超伝導臨界温度 Tc を Thouless 判定条件により求めます。

```bash
cargo run --example ladder_resummation
```

出力例（全出力）:

```
=== Ladder Resummation: Attractive Hubbard Model ===
Lattice: 8x8 square, t = 1, U = -2

Searching for Tc via Thouless criterion (n_matsubara = 256)...
Tc = 0.223253

=== Thouless Criterion: 1 - U * chi_0(q=0, nu=0) ===
T            chi_0(q=0)           1 - U*chi_0         
----------------------------------------------------
0.178603     -0.56792543          -0.13585086         
0.200928     -0.53070584          -0.06141169         
0.212091     -0.51465041          -0.02930081         
0.223253     -0.49996228          0.00007544          
0.234416     -0.48645102          0.02709796          
0.245579     -0.47396132          0.05207736          
0.334880     -0.39902075          0.20195849          
0.446507     -0.33946158          0.32107685          

At T = Tc, the Thouless criterion 1 - U*chi_0 should vanish,
signaling the onset of superconducting instability.
```

**物理的な意味:** Thouless 判定条件は 1 − Uχ₀(q=0, iν₀=0) = 0 で与えられ、この式がゼロになる温度が超伝導臨界温度 Tc です。引力 Hubbard 模型（U < 0）では粒子-粒子（pp）チャネルに Cooper 不安定性が生じ、T を下げていくとこの条件が満たされ、対凝縮の onset が表れます。

### 自己エネルギーダイアグラム

1 次および 2 次の自己エネルギーダイアグラムを生成・分類します。

```bash
cargo run --example second_order_self_energy
```

出力例:

```
=== Attractive Hubbard Model ===
Lattice: 4x4 square, t = 1, U = -2

1st order: 2 Wick contractions
2nd order: 12 Wick contractions

1st order: 2 unique topologies
2nd order: 7 unique topologies

=== 1st Order Self-Energy Diagrams ===
Topology 1 (weight = 1):
  Order 1: 2 loop(s), 2 internal propagators
  Prefactor: 2.0000
  Expression: Σ_k1,ω1 [Σ_k2,ω2 [(2 × G₀(k1, iω1) × G₀(k2, iω2))]]

Topology 2 (weight = -1):
  Order 1: 1 loop(s), 1 internal propagators
  Prefactor: -2.0000
  Expression: Σ_k1,ω1 [(-2 × G₀(k1, iω1))]

=== 2nd Order Self-Energy Diagrams ===
... 7 topologies (4 internal propagators × 3, or 3 internal × 4) ...
```

続いて Graphviz 用の DOT 形式で全トポロジーの図が標準出力に出力されます。

**各トポロジーの物理的意味:**

- **1 次 Topology 1（weight=1）**: 2 本のループと 2 本の内線。Hartree 型の寄与（2 つの伝播線のループ）。
- **1 次 Topology 2（weight=-1）**: 1 本のループと 1 本の内線。Fock 型の寄与（交換で符号が反転）。
- **2 次の寄与**: 4 頂点・複数ループのダイアグラム群。繰り込みや散乱の 2 次効果を表します。

**図のエクスポート:** `--dot` と `--json` フラグでファイルに直接出力できます。

```bash
# DOT を出力して PNG に変換
cargo run --example second_order_self_energy -- --dot diagrams.dot
dot -Tpng diagrams.dot -o diagrams.png

# D3.js 用 JSON を出力
cargo run --example second_order_self_energy -- --json diagrams.json

# 両方同時に出力
cargo run --example second_order_self_energy -- --dot diagrams.dot --json diagrams.json

# DOT を標準出力に書き出す（テキスト要約は抑制）
cargo run --example second_order_self_energy -- --dot -
```

スピン上向き伝播線は青、下向きは赤、外線は破線で描画されます。

JSON 出力には `topology_id`、`weight`、`vertices`、`propagators`（`from`、`to`、`spin`、`external`）、`sign`、`symmetry_factor` フィールドが含まれます。

## アーキテクチャ

データは 2 本のパイプラインで下から上に流れます。

**記号パイプライン**（図の生成）:

```
FermionOperator → [Wick's theorem] → WickTerm (縮約 + 符号)
    → [generate_diagrams] → FeynmanDiagram (グラフ表現)
    → [classify_diagrams] → 重み付きユニークトポロジー
    → [apply_feynman_rules] → FeynmanExpression (記号 Expr ツリー)
```

**数値パイプライン**（再総和）:

```
HubbardModel + ThermalParams → evaluate_g0(k, iωₙ)
    → compute_pp_susceptibility(q, iνₘ)  [χ₀ = pp バブル]
    → solve_tmatrix / find_tc            [T = U/(1-Uχ₀), Thouless 判定]
```

### モジュール一覧


| モジュール            | 説明                                                   |
| ---------------- | ---------------------------------------------------- |
| `algebra/`       | フェルミオン演算子、スピン、縮約、Wick の定理                            |
| `models/`        | 正方格子（k グリッド、分散）、Hubbard 模型パラメータ                      |
| `diagrams/`      | Feynman 図グラフ、生成、トポロジー分類、Feynman 規則                   |
| `symbolic/`      | 式ツリー（`G0`, `U`, `Sum`, `Mul`, `Add`, ...）と `Display` |
| `resummation/`   | PP 梯子チャネル、T 行列、Thouless 判定 Tc 探索                     |
| `numerical/`     | グリーン関数 G0(k, iωₙ)、松原和、PP 感受率                         |
| `visualization/` | Graphviz DOT 出力（スピン↑=青、スピン↓=赤）                       |


## 物理的背景

### Hubbard 模型

Hubbard 模型は格子上の電子のホッピングとオンサイト相互作用を記述します。ハミルトニアンは

**H = −t Σ_{****σ} c†*iσ c_jσ + U Σ_i n*{i↑} n_{i↓}**

で定義されます。2 次元正方格子では分散は ε(k) = −2t(cos kx + cos ky) です。相互作用頂点は同一サイトの 4 つのフェルミオン演算子 c†↑ c↑ c†↓ c↓ を結合し、結合定数が U です。U < 0（引力）のとき s 波超伝導が可能で、U > 0（斥力）のときは粒子-粒子チャネルは発散しません。

### 松原形式

有限温度 T = 1/β では、時間順序の摂動論は虚時間で定式化されます。フェルミオンの松原周波数は **iωₙ = i(2n+1)πT** です。裸のグリーン関数は G₀(k, iωₙ) = 1/(iωₙ − εₖ + μ) で、すべての内線の運動量・周波数和はブリルアンゾーンと離散松原周波数で実行されます。

### PP 感受率と梯子型再総和

粒子-粒子バブル（裸の感受率）は

**χ₀(q, iνₘ) = −(1/Nβ) Σₖ Σₙ G₀(k, iωₙ) G₀(q−k, iνₘ−iωₙ)**

で定義されます。物理的には、2 本の伝播線でつながった粒子-粒子の「泡」で、Cooper 対の揺らぎの裸の応答です。梯子型再総和は、粒子-粒子散乱の繰り返しを T 行列 T(q, iνₘ) = U / (1 − U χ₀(q, iνₘ)) によりすべての次数で足し上げます。これは粒子-粒子チャネルでの Bethe-Salpeter 方程式を解くことに相当します。

### Thouless 判定条件

超伝導不安定性は T 行列が発散するとき、すなわち **1 − U χ₀(q=0, iν₀=0) = 0** のときに起こります。これが Thouless 判定条件です。臨界温度 Tc は、この条件が満たされるまで温度を走査する二分法で求めます。引力 Hubbard（U < 0）では χ₀ < 0 なので Uχ₀ > 0 となり、1 − |U||χ₀| = 0 で発散し、Cooper 不安定性と対凝縮の onset を表します。

## モジュール解説

- **algebra**: `FermionOperator`, `Spin`, `Contraction`, `WickTerm`, `wick_theorem()`。スピン保存でマッチングを制限（↑同士・↓同士のみ）。フェルミオン符号は反転数で計算。
- **models**: `SquareLattice`（k グリッド、分散 ε(k)=-2t(cos kx+cos ky)）、`HubbardModel`（t, U）、`ThermalParams`（β, μ, n_matsubara）。
- **diagrams**: 頂点・伝播線からなる `FeynmanDiagram` グラフ。`generate_diagrams()` は `Observable`（Vacuum, SelfEnergy）を引数に取る。`classify_diagrams()` は頂点記述子のソートされた署名でトポロジーを判定。`apply_feynman_rules()` で記号 `Expr` ツリーの `FeynmanExpression` を生成。
- **symbolic**: `Expr` 列挙型（G0, U, Sum, Mul, Add, Inv, Scalar, Neg）と数式表示用 `Display`。
- **resummation**: PP 梯子チャネル `PPLadder`。`solve_tmatrix()` で T = U/(1−Uχ₀) を計算。`find_tc()` で Thouless 条件 1−Uχ₀(T)=0 を二分法で解く。
- **numerical**: `evaluate_g0()` で G₀(k,iωₙ)=1/(iωₙ−εₖ+μ)。`compute_pp_susceptibility()` で k と松原周波数の和。ボゾン周波数 iνₘ−iωₙ のインデックスはフェルミオン側で m−1−n に対応。
- **visualization**: `to_dot()` / `to_dot_all()` で Graphviz DOT を出力。スピン↑=青、スピン↓=赤、外線=破線。

