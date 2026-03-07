# feynman-engine 使い方ガイド

## 概要

`feynman-engine` は、引力Hubbard模型の有限温度における摂動展開を自動化するRustライブラリである。フェルミオン演算子からWick定理を適用して全ての縮約をスピン保存則とフェルミオン符号付きで生成し、Feynmanダイアグラムを構築してトポロジーごとに分類した上で、Feynman rulesにより記号式（式木）へ変換する。数値計算側では、松原形式で粒子-粒子（PP）感受率を評価し、梯子型再総和（Bethe-Salpeter方程式）を通じてT行列を求め、Thouless判定条件 1 - U*chi_0 = 0 を二分法で解くことにより超伝導転移温度 Tc を決定する。

## 機能一覧

- Wick縮約の自動生成（スピン保存則によるフィルタリング、置換の反転数によるフェルミオン符号計算）
- Feynmanダイアグラムの生成とトポロジー分類（頂点記述子の正規化による一意な同定）
- 記号式ツリー（Expr）とFeynman rulesによる自動変換
- PP梯子型再総和（Bethe-Salpeter方程式によるT行列計算）
- Thouless判定条件による Tc ファインダー（二分法）
- Graphviz DOT形式での可視化（スピン上=青、スピン下=赤、外線=破線）

## インストール

Rust 2021 editionが必要である。依存クレート: `num-complex`, `rayon`, `itertools`。

```bash
# リポジトリをクローン
git clone <repo-url>
cd feynman

# ライブラリとバイナリをビルド
cargo build

# 全テスト（ユニットテスト＋結合テスト）を実行（約0.5秒）
cargo test
```

特定モジュールのテストのみ実行する場合:

```bash
# algebra::wick モジュールのテストを実行
cargo test algebra::wick

# 結合テストのみ実行
cargo test --test test_integration
```

## 使い方

### 梯子型再総和 (ladder resummation)

2D正方格子上の引力Hubbard模型に対し、PP梯子型再総和により超伝導転移温度 Tc を求める。

```bash
cargo run --example ladder_resummation
```

出力:

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

**物理的意味**: Thouless判定条件 1 - U*chi_0(q=0, i*nu_0=0) = 0 が満たされる温度が超伝導転移温度 Tc である。引力Hubbard模型（U < 0）では、PP感受率 chi_0 が負の値を取るため、U*chi_0 = |U|*|chi_0| > 0 となる。温度を下げると |chi_0| が増大し、|U|*|chi_0| = 1 に達した時点でT行列が発散する。これはCooper不安定性の発現、すなわちs波超伝導秩序の形成を示す。上の出力では、T = 0.223253 付近で 1 - U*chi_0 がほぼ0に近づいていることが確認できる。

### 自己エネルギーダイアグラム

Hubbard模型の1次および2次自己エネルギーダイアグラムを生成し、トポロジーごとに分類して記号式を出力する。

```bash
cargo run --example second_order_self_energy
```

出力（DOT出力部分は省略）:

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
Topology 1 (weight = 1):
  Order 2: 3 loop(s), 4 internal propagators
  ...
(7つのトポロジーが出力される)
```

**各トポロジーの物理的意味**:

- **1次 Topology 1**: Hartree項に対応する。2つの内部伝播関数からなるループを含み、自己エネルギーへの平均場的な寄与を表す。
- **1次 Topology 2**: Fock項（交換項）に対応する。1つの内部伝播関数によるタドポール型のダイアグラムであり、交換相互作用による自己エネルギー補正を表す。
- **2次トポロジー（7種類）**: 2次の摂動寄与であり、2つの相互作用頂点を含む。粒子-粒子チャンネル、粒子-ホールチャンネルなど異なる散乱過程に対応するダイアグラムが分類されている。

**Graphvizによる可視化**: DOT出力をファイルに保存し、Graphvizで画像化できる。

```bash
# 例: DOT出力をファイルにリダイレクト（プログラム出力のDOT部分を保存）
cargo run --example second_order_self_energy > output.txt

# DOTファイルをPNG画像に変換
dot -Tpng output.dot -o diagram.png

# PDF形式で出力する場合
dot -Tpdf output.dot -o diagram.pdf
```

可視化されたダイアグラムでは、スピン上の伝播関数が青、スピン下が赤で表示され、外部脚は破線で描画される。

## アーキテクチャ

データは以下の2つのパイプラインを通じてボトムアップに処理される。

**記号式パイプライン**（ダイアグラム生成）:

```
フェルミオン演算子 → [Wick定理] → WickTerm (縮約 + 符号)
    → [generate_diagrams] → FeynmanDiagram (グラフ表現)
    → [classify_diagrams] → 一意なトポロジー + 重み
    → [apply_feynman_rules] → FeynmanExpression (記号式ツリー)
```

**数値パイプライン**（再総和）:

```
HubbardModel + ThermalParams → evaluate_g0(k, iωₙ)
    → compute_pp_susceptibility(q, iνₘ)  [χ₀ = pp bubble]
    → solve_tmatrix / find_tc            [T = U/(1-Uχ₀), Thouless判定条件]
```

### モジュール一覧

| モジュール | 説明 |
|---|---|
| `algebra/` | フェルミオン演算子、スピン、縮約、Wick定理の実装 |
| `models/` | 正方格子（k点グリッド、分散関係）、Hubbard模型のパラメータ |
| `diagrams/` | Feynmanダイアグラムのグラフ表現、生成、トポロジー分類、Feynman rules |
| `symbolic/` | 式木（`G0`, `U`, `Sum`, `Mul`, `Add` 等）と `Display` による可読な数式出力 |
| `resummation/` | PP梯子型チャンネル、T行列、Thouless判定条件による Tc ファインダー |
| `numerical/` | Green関数 G0(k, iwn) の評価、松原周波数和、PP感受率の計算 |
| `visualization/` | Graphviz DOT形式でのエクスポート（スピン上=青、スピン下=赤） |

## 物理的背景

### Hubbard模型

Hubbard模型は、格子上の電子のホッピングとオンサイト相互作用を記述する最も基本的な強相関電子系の模型である。ハミルトニアンは以下の通り:

```
H = -t Σ_{<i,j>,σ} c†_{i,σ} c_{j,σ} + U Σ_i n_{i,↑} n_{i,↓}
```

ここで t はホッピング積分、U はオンサイト相互作用、c†/c は生成・消滅演算子、n = c†c は数演算子である。

2D正方格子上では、分散関係は以下のようになる:

```
ε(k) = -2t(cos kx + cos ky)
```

相互作用頂点は同一サイト・同一時間の4つのフェルミオン演算子 c†_↑ c_↑ c†_↓ c_↓ を結合度 U で結びつける。

- **U < 0（引力）**: s波超伝導が発現する。PP梯子型チャンネルにCooper不安定性が生じる。
- **U > 0（斥力）**: PPチャンネルには不安定性は生じない（反強磁性などの秩序が別チャンネルに現れる）。

### 松原形式

有限温度 T = 1/beta における摂動論は、虚時間形式（松原形式）で定式化される。

フェルミオンの松原周波数は離散的であり:

```
iωₙ = i(2n+1)π/β   (n = 0, ±1, ±2, ...)
```

裸のGreen関数は:

```
G₀(k, iωₙ) = 1 / (iωₙ - εₖ + μ)
```

ここで μ は化学ポテンシャルである。全ての内部運動量和はブリルアンゾーン上で、周波数和は離散松原周波数の集合上で実行される。

### PP感受率（粒子-粒子バブル）

粒子-粒子感受率（裸の感受率、PPバブル）は以下で定義される:

```
χ₀(q, iνₘ) = -1/(Nβ) Σₖ Σₙ G₀(k, iωₙ) G₀(q-k, iνₘ - iωₙ)
```

ここで N は格子点数、iνₘ はボソン型松原周波数である。この量は、運動量 q、周波数 iνₘ を持つCooper対が形成される確率振幅に対応する。

### 梯子型再総和とT行列

PPチャンネルにおける繰り返し散乱を全次数にわたって足し上げると、T行列が得られる:

```
T(q, iνₘ) = U / (1 - U × χ₀(q, iνₘ))
```

これは粒子-粒子チャンネルにおけるBethe-Salpeter方程式の解に等しい。

### Thouless判定条件

超伝導不安定性は、T行列が発散する条件:

```
1 - U × χ₀(q=0, iν₀=0) = 0
```

によって決定される。これがThouless判定条件である。

引力Hubbard模型（U < 0）では χ₀ < 0 であるため:

```
U × χ₀ = |U| × |χ₀| > 0
```

となり、温度を下げて |χ₀| が増大し |U| × |χ₀| = 1 に達すると T行列が発散する。この発散はCooper不安定性を示し、対凝縮（超伝導秩序）の発現に対応する。

転移温度 Tc は、1 - U*χ₀(T) = 0 を温度 T について二分法で解くことにより数値的に決定される。

## モジュール解説

### algebra

Wick定理に基づく代数的操作を実装する。

- **`FermionOperator`**: フェルミオンの生成・消滅演算子を表す構造体。サイトインデックス、スピン、生成/消滅の種別を保持する。
- **`Spin`**: スピンの上（Up）・下（Down）を表す列挙型。
- **`Contraction`**: 2つのフェルミオン演算子間の縮約（伝播関数に対応）を表す。
- **`WickTerm`**: 縮約の集合とフェルミオン符号（+1 または -1）の組。Wick定理の各項に対応する。
- **`wick_theorem()`**: フェルミオン演算子の列に対してWick定理を適用し、全ての可能な縮約パターンを生成する。スピン保存則により、同じスピンの演算子同士のみが縮約される。フェルミオン符号は置換の反転数から計算される。

### models

物理模型のパラメータを定義する。

- **`SquareLattice`**: 2D正方格子。k点グリッドと分散関係 ε(k) = -2t(cos kx + cos ky) を提供する。
- **`HubbardModel`**: Hubbard模型のパラメータ（ホッピング t、相互作用 U）と格子を保持する。
- **`ThermalParams`**: 有限温度計算のパラメータ。逆温度 β、化学ポテンシャル μ、松原周波数の打ち切り数 n_matsubara を指定する。

### diagrams

Feynmanダイアグラムの生成・分類・Feynman rules適用を行う。

- **`FeynmanDiagram`**: 頂点（`Vertex`）と伝播関数（`Propagator`）からなるグラフ構造でダイアグラムを表現する。
- **`generate_diagrams()`**: `Observable` 列挙型（`Vacuum` または `SelfEnergy`）に応じて、指定次数のFeynmanダイアグラムを生成する。
- **`classify_diagrams()`**: 頂点記述子の正規化された署名を用いて、トポロジーが等価なダイアグラムをグループ化し、各トポロジーの重みを計算する。
- **`apply_feynman_rules()`**: ダイアグラムにFeynman rulesを適用し、`FeynmanExpression`（記号式ツリーと前因子）を生成する。

### symbolic

記号式を表現するための式木。

- **`Expr`**: 記号式を表す列挙型。以下のバリアントを持つ:
  - `G0` -- 裸のGreen関数
  - `U` -- 相互作用頂点
  - `Sum` -- 運動量・周波数に関する和
  - `Mul` -- 積
  - `Add` -- 和
  - `Inv` -- 逆数
  - `Scalar` -- スカラー定数
  - `Neg` -- 符号反転

  `Display` トレイトの実装により、人間が読みやすい数式形式で出力される。

### resummation

梯子型再総和を実装する。

- **`PPLadder`**: 粒子-粒子梯子型チャンネルを表す。
- **`solve_tmatrix()`**: T行列 T = U/(1 - U*chi_0) を計算する。
- **`find_tc()`**: Thouless判定条件 1 - U*chi_0(T) = 0 を二分法で解き、転移温度 Tc を返す。温度範囲と松原周波数の打ち切り数、収束判定の閾値を引数に取る。

### numerical

数値計算の基本関数を提供する。

- **`evaluate_g0()`**: 裸のGreen関数 G₀(k, iωₙ) = 1/(iωₙ - εₖ + μ) を評価する。
- **`compute_pp_susceptibility()`**: PP感受率 χ₀(q, iνₘ) = -1/(Nβ) Σₖ Σₙ G₀(k, iωₙ) G₀(q-k, iνₘ-iωₙ) を計算する。ボソン松原周波数 iνₘ に対するフェルミオン周波数のインデックス変換 m-1-n を内部で処理する。

### visualization

ダイアグラムの可視化機能を提供する。

- **`to_dot()`**: 単一のFeynmanダイアグラムをGraphviz DOT形式の文字列に変換する。
- **`to_dot_all()`**: 分類済みの全ダイアグラムをまとめてDOT形式に出力する。スピン上の伝播関数は青、スピン下は赤で色分けされ、外部脚は破線で表示される。
