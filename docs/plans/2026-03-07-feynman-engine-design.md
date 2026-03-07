# FeynmanEngine.jl 設計ドキュメント

## 概要

Julia製の量子多体系シンボリック摂動計算エンジン。Wickの定理に基づくファインマンダイアグラムの自動生成・分類・数式導出に加え、梯子型ダイアグラム等の無限部分和（リサメーション）を等比級数として取り込む機能を持つ。

## 動機と目標

- **対象模型:** 単バンドハバード模型（まず引力ハバード）
- **摂動展開:** 弱結合展開（自由電子 G₀ を基準、U で展開）
- **温度:** 有限温度（松原形式）を基本とし、ゼロ温度はその極限
- **物理的目標:** 粒子-粒子チャンネルの梯子ダイアグラムを足し上げ、クーパー対不安定性・超伝導転移温度 Tc を求める

## アプローチ

ボトムアップ方式：Wickの定理をシンボリックに実装し、任意次数のダイアグラムを自動生成。繰り返し構造を検出してBethe-Salpeter方程式として閉じた形に変換する。

## 出力

1. 解析的なシンボリック式（ギャップ方程式等）
2. 数値評価可能な関数（格子上での Tc 計算等）
3. ダイアグラムの可視化

## モジュール構成

```
feynman/
├── src/
│   ├── FeynmanEngine.jl          # パッケージエントリポイント
│   │
│   ├── algebra/                   # コア層：演算子代数
│   │   ├── operators.jl           # 生成・消滅演算子の型定義
│   │   ├── normal_ordering.jl     # 正規順序化
│   │   └── wick.jl                # Wickの定理（縮約の自動列挙）
│   │
│   ├── diagrams/                  # ダイアグラム層
│   │   ├── graph.jl               # ファインマンダイアグラムのグラフ表現
│   │   ├── generate.jl            # Wick縮約 → ダイアグラム変換
│   │   ├── classify.jl            # トポロジー分類（同型判定）
│   │   └── rules.jl               # ファインマンルール（ダイアグラム → 数式）
│   │
│   ├── resummation/               # リサメーション層
│   │   ├── patterns.jl            # 繰り返し構造のパターンマッチ
│   │   ├── bethe_salpeter.jl      # Bethe-Salpeter方程式の構築
│   │   └── channels.jl            # pp/phチャンネルの定義
│   │
│   ├── models/                    # 模型定義
│   │   ├── lattice.jl             # 格子構造（正方格子など）
│   │   └── hubbard.jl             # ハバード模型のハミルトニアン・G₀
│   │
│   ├── numerical/                 # 数値評価層
│   │   ├── greens_function.jl     # G₀(k, iωₙ) の数値評価
│   │   ├── matsubara.jl           # 松原周波数の和
│   │   └── evaluate.jl            # シンボリック式 → 数値関数変換
│   │
│   └── visualization/             # 可視化層
│       ├── diagram_plot.jl        # ダイアグラムの描画
│       └── spectral.jl            # スペクトル関数などの物理量プロット
│
├── test/
├── docs/
├── examples/
│   ├── second_order_self_energy.jl
│   └── ladder_resummation.jl
└── Project.toml
```

## セクション1：コア層 — 演算子代数とWickの定理

### 型定義

```julia
struct FermionOperator
    creation::Bool       # true=生成(c†), false=消滅(c)
    site::Symbol         # サイトラベル :i, :j, ...
    spin::Symbol         # :↑ or :↓
    time::Symbol         # 虚時間ラベル :τ₁, :τ₂, ...
end

struct Contraction
    creator::FermionOperator
    annihilator::FermionOperator
end

struct WickTerm
    sign::Int
    contractions::Vector{Contraction}
end
```

### Wickの定理アルゴリズム

1. 演算子列から生成・消滅演算子を分離
2. 全完全マッチングを列挙（スピン保存則でフィルタ）
3. 各マッチングのフェルミオン符号を計算
4. WickTermのリストを返す

ハバード模型ではスピン保存により組み合わせ数が大幅削減。

### 摂動項の生成

n次摂動 → (-U)^n/n! × (4n個の演算子のWick縮約)。各頂点は同一サイト上の c†↑ c↑ c†↓ c↓。

## セクション2：ダイアグラム層

### グラフ表現

```julia
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
    symmetry_factor::Rational
end
```

### トポロジー分類

Graphs.jl の彩色グラフ同型判定により、同じトポロジーのダイアグラムをグループ化。各グループの重み（出現回数 × 符号）を計算。

### ファインマンルール

ダイアグラムから数式を生成:
1. 各内部伝播線 → G₀(k, iωₙ)
2. 各頂点 → U
3. 各内部ループに Σₖ Σ_ωₙ
4. フェルミオン符号 × (-1)^(ループ数)
5. 対称性因子

独立な内部変数の数 = ループ数 = 辺数 - 頂点数 + 1（連結の場合）。

## セクション3：リサメーション層

### 梯子ダイアグラムの無限和

n次の梯子は [U·χ₀]^(n-1)·U の構造を持つ。等比級数として:

```
T = U + U·χ₀·U + U·χ₀·U·χ₀·U + ...
  = U / (1 - U·χ₀)
```

### パターンマッチ

```julia
struct RepeatUnit
    kernel::FeynmanDiagram
    channel::Symbol          # :pp, :ph, :ph_crossed
    external_legs::Tuple{Int,Int}
end
```

n次とn+1次のダイアグラムを比較し、繰り返し挿入される部分グラフ（repeat unit）を抽出。

### チャンネル

- **ppチャンネル（粒子-粒子）:** χ₀(q,iνₘ) = -(1/Nβ) Σₖ Σ_ωₙ G₀(k,iωₙ) G₀(q-k, iνₘ-iωₙ)
- **phチャンネル（粒子-正孔）:** 将来の拡張用

### Bethe-Salpeter方程式

```julia
struct BetheSalpeterEquation
    channel::Symbol
    irreducible_vertex::Any   # 最低次は U
    bare_susceptibility::Any  # χ₀
end
```

T = Γ + Γ·χ₀·T を解く。クーパー対不安定性は 1 - U·χ₀(T) = 0（Thouless criterion）で判定。

### 実装フェーズ

- **フェーズ1:** チャンネルをユーザーが指定し、その梯子を自動リサム
- **フェーズ2:** Wick定理で生成したダイアグラム群から繰り返し構造を自動検出

## セクション4：数値評価層

### 松原グリーン関数

G₀(k, iωₙ) = 1/(iωₙ - εₖ + μ)、ωₙ = (2n+1)π/β

```julia
struct ThermalParameters
    β::Float64
    μ::Float64
    n_matsubara::Int
end
```

### 高速化

- **高周波テイル補正:** 大きい |ωₙ| での解析的処理で収束加速
- **FFT畳み込み:** χ₀ を虚時間・実空間で計算し FFT 変換。O(N²) → O(N log N)

### シンボリック → 数値変換

Symbolics.build_function() でシンボリック式をJulia関数にコンパイル。

## セクション5：可視化層

- **ダイアグラム描画:** Makie.jl で頂点（●）、伝播線（→）、相互作用線（波線）を描画
- **物理量プロット:** スペクトル関数 A(k,ω)、感受率の温度依存性、Thouless criterionの図示

## 依存パッケージ

- Symbolics.jl — シンボリック計算
- Graphs.jl — グラフ表現・同型判定
- Makie.jl — 可視化
- FFTW.jl — 波数空間のFFT計算

## 使用例イメージ

```julia
using FeynmanEngine

lattice = SquareLattice(16, 16)
model = HubbardModel(lattice, t=1.0, U=-2.0)
params = ThermalParameters(β=10.0, μ=0.0, n_matsubara=512)

# 2次までのダイアグラムを自動生成・表示
diagrams = generate_diagrams(model, order=2, observable=:self_energy)
plot_all_diagrams(diagrams, 2)

# 梯子ダイアグラムのリサメーション → Tc 探索
bse = build_pp_ladder(model)
Tc = find_Tc(model, lattice, T_range=(0.01, 1.0))
println("Tc = $Tc")

# T行列と感受率の温度依存性を可視化
plot_susceptibility_vs_temperature(model, 0.01:0.01:1.0)
```
