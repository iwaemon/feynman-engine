# FeynmanEngine 設計ドキュメント (Rust版)

## 概要

Rust製の量子多体系摂動計算エンジン。Wickの定理に基づくファインマンダイアグラムの自動生成・分類・数式導出に加え、梯子型ダイアグラム等の無限部分和（リサメーション）を等比級数として取り込む機能を持つ。

## 動機と目標

- **対象模型:** 単バンドハバード模型（まず引力ハバード）
- **摂動展開:** 弱結合展開（自由電子 G₀ を基準、U で展開）
- **温度:** 有限温度（松原形式）を基本とし、ゼロ温度はその極限
- **物理的目標:** 粒子-粒子チャンネルの梯子ダイアグラムを足し上げ、クーパー対不安定性・超伝導転移温度 Tc を求める

## アプローチ

ボトムアップ方式：Wickの定理をシンボリックに実装し、任意次数のダイアグラムを自動生成。繰り返し構造を検出してBethe-Salpeter方程式として閉じた形に変換する。

## 出力

1. 解析的なシンボリック式（式木AST → 人間可読な数式表示）
2. 数値評価可能な関数（格子上での Tc 計算等、rayonで並列化）
3. ダイアグラムの可視化（Graphviz DOT出力 + plottersによるプロット）

## モジュール構成

```
feynman/
├── src/
│   ├── main.rs                    # CLI エントリポイント
│   ├── lib.rs                     # ライブラリエントリポイント
│   │
│   ├── algebra/                   # コア層：演算子代数
│   │   ├── mod.rs
│   │   ├── operators.rs           # 生成・消滅演算子の型定義
│   │   └── wick.rs                # Wickの定理（縮約の自動列挙）
│   │
│   ├── diagrams/                  # ダイアグラム層
│   │   ├── mod.rs
│   │   ├── graph.rs               # ファインマンダイアグラムのグラフ表現
│   │   ├── generate.rs            # Wick縮約 → ダイアグラム変換
│   │   ├── classify.rs            # トポロジー分類（同型判定）
│   │   └── rules.rs               # ファインマンルール（ダイアグラム → 式木）
│   │
│   ├── symbolic/                  # シンボリック式木
│   │   ├── mod.rs
│   │   ├── expr.rs                # 式木（AST）の定義
│   │   ├── simplify.rs            # 式の簡約化
│   │   └── display.rs             # 人間可読な数式表示
│   │
│   ├── resummation/               # リサメーション層
│   │   ├── mod.rs
│   │   ├── channels.rs            # pp/phチャンネルの定義
│   │   └── bethe_salpeter.rs      # Bethe-Salpeter方程式の構築
│   │
│   ├── models/                    # 模型定義
│   │   ├── mod.rs
│   │   ├── lattice.rs             # 格子構造（正方格子など）
│   │   └── hubbard.rs             # ハバード模型のハミルトニアン・G₀
│   │
│   ├── numerical/                 # 数値評価層
│   │   ├── mod.rs
│   │   ├── greens_function.rs     # G₀(k, iωₙ) の数値評価
│   │   └── matsubara.rs           # 松原周波数の和、χ₀計算
│   │
│   └── visualization/             # 可視化層
│       ├── mod.rs
│       ├── dot.rs                 # Graphviz DOT出力
│       └── plot.rs                # plottersによる物理量プロット
│
├── tests/                         # 統合テスト
│   ├── test_wick.rs
│   ├── test_diagrams.rs
│   ├── test_resummation.rs
│   └── test_integration.rs
│
├── examples/                      # 使用例
│   ├── second_order_self_energy.rs
│   └── ladder_resummation.rs
│
├── Cargo.toml
└── README.md
```

## セクション1：コア層 — 演算子代数とWickの定理

### 型定義

```rust
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FermionOperator {
    pub creation: bool,     // true=生成(c†), false=消滅(c)
    pub site: String,       // サイトラベル "i", "j", ...
    pub spin: Spin,         // Up or Down
    pub time: String,       // 虚時間ラベル "τ₁", "τ₂", ...
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Spin {
    Up,
    Down,
}

#[derive(Debug, Clone)]
pub struct Contraction {
    pub creator: FermionOperator,
    pub annihilator: FermionOperator,
}

#[derive(Debug, Clone)]
pub struct WickTerm {
    pub sign: i32,          // +1 or -1
    pub contractions: Vec<Contraction>,
}
```

### Wickの定理アルゴリズム

1. 演算子列から生成・消滅演算子を分離
2. スピン保存則でフィルタし全完全マッチングを列挙
3. 各マッチングのフェルミオン符号を計算（置換の符号）
4. WickTermのVecを返す

## セクション2：シンボリック式木（AST）

Julia版の Symbolics.jl の代わりに、限定的だが十分な式木を自前実装。

```rust
#[derive(Debug, Clone)]
pub enum Expr {
    /// 自由グリーン関数 G₀(k, iωₙ)
    G0 { k: String, omega: String },
    /// 相互作用頂点 U
    U,
    /// 松原周波数の和 (1/Nβ) Σ_{var}
    Sum { var: String, body: Box<Expr> },
    /// 積
    Mul(Vec<Expr>),
    /// 和
    Add(Vec<Expr>),
    /// 逆数 1/x
    Inv(Box<Expr>),
    /// スカラー値
    Scalar(f64),
    /// 負符号
    Neg(Box<Expr>),
}

impl std::fmt::Display for Expr {
    // 人間可読な数式表示
    // 例: "(-U)² × (1/Nβ) Σₖ Σ_ωₙ G₀(k,iωₙ) G₀(k+q,iωₙ+iνₘ)"
}
```

## セクション3：ダイアグラム層

```rust
#[derive(Debug, Clone)]
pub struct Vertex {
    pub id: usize,
    pub site: String,
    pub time: String,
}

#[derive(Debug, Clone)]
pub struct Propagator {
    pub from: usize,        // vertex id
    pub to: usize,          // vertex id
    pub spin: Spin,
    pub external: bool,
}

#[derive(Debug, Clone)]
pub struct FeynmanDiagram {
    pub order: usize,
    pub vertices: Vec<Vertex>,
    pub propagators: Vec<Propagator>,
    pub sign: i32,
    pub symmetry_factor: f64,
}
```

トポロジー分類はソート済み頂点記述子による正規形で判定。

## セクション4：リサメーション層

梯子ダイアグラムの無限和:
```
T = U + U·χ₀·U + U·χ₀·U·χ₀·U + ... = U / (1 - U·χ₀)
```

```rust
pub struct PPLadder {
    pub model: HubbardModel,
}

impl PPLadder {
    pub fn solve_tmatrix(&self, params: &ThermalParams, q: &[f64; 2], m: i32) -> Complex64;
}

pub fn find_tc(model: &HubbardModel, t_range: (f64, f64), n_matsubara: usize) -> f64;
```

フェーズ1: チャンネルをユーザーが指定し梯子を自動リサム
フェーズ2: ダイアグラム群から繰り返し構造を自動検出

## セクション5：数値評価層

```rust
pub fn evaluate_g0(model: &HubbardModel, params: &ThermalParams,
                   k: &[f64; 2], n: i32) -> Complex64 {
    let omega_n = matsubara_freq(n, params.beta);
    let eps_k = model.dispersion(k);
    Complex64::new(0.0, omega_n) - eps_k + params.mu
    // → 1.0 / above
}
```

**並列化:** rayonでk点ループを自動並列化。FFTはrustfftクレート。

## セクション6：可視化層

- ダイアグラム → Graphviz DOT形式で出力、`dot -Tsvg` で描画
- 物理量プロット → plottersクレートでPNG出力、またはCSV出力

## 依存クレート

- `num-complex` — 複素数
- `rayon` — 並列計算
- `rustfft` — FFT
- `plotters` — プロット（オプション）
- `itertools` — 組み合わせ列挙

## 使用例イメージ

```rust
use feynman_engine::prelude::*;

fn main() {
    let lattice = SquareLattice::new(16, 16);
    let model = HubbardModel::new(lattice, 1.0, -2.0); // t=1, U=-2

    // 2次までのダイアグラムを自動生成
    let diagrams = generate_diagrams(&model, 2, Observable::SelfEnergy);
    let classified = classify_diagrams(&diagrams);
    for (diag, weight) in &classified {
        let expr = apply_feynman_rules(diag, &model);
        println!("weight={}, expr={}", weight, expr);
    }

    // ダイアグラム可視化
    let dot = to_dot(&classified[0].0);
    std::fs::write("diagram.dot", dot).unwrap();

    // 梯子リサメーション → Tc 探索
    let tc = find_tc(&model, (0.01, 2.0), 512);
    println!("Tc = {:.4}", tc);
}
```
