# はじめに

## feynman-engine の概要

`feynman-engine` は、有限温度における引力Hubbardモデルの摂動的多体計算を自動化するRustライブラリです。裸のフェルミオン演算子から出発し、Wickの定理（Wick's theorem）を適用してフェルミオン符号とスピン保存を正しく考慮した全ての縮約を生成し、ファインマンダイアグラムを構成し、トポロジーによって分類し、標準的なファインマンルールを用いて記号的な数式へと変換します。

数値計算の面では、松原形式（Matsubara formalism）における粒子-粒子感受率を評価し、ラダー再足し上げ（Bethe-Salpeter方程式の解法）を行うことで、Thouless判定条件（Thouless criterion）を通じて超伝導転移温度 Tc を求めます。

## 2つのパイプライン

本エンジンは、2つの相補的なパイプラインで構成されています。

**記号パイプライン**（ダイアグラム生成）:

```text
FermionOperator -> [Wick's theorem] -> WickTerm (contractions + sign)
    -> [generate_diagrams] -> FeynmanDiagram (graph representation)
    -> [classify_diagrams] -> unique topologies with weights
    -> [apply_feynman_rules] -> FeynmanExpression (symbolic Expr tree)
```

**数値パイプライン**（再足し上げ）:

```text
HubbardModel + ThermalParams -> evaluate_g0(k, iwn)
    -> compute_pp_susceptibility(q, inu_m)  [chi_0 = pp bubble]
    -> solve_tmatrix / find_tc              [T = U/(1-U*chi_0), Thouless criterion]
```

## モジュール一覧

| モジュール | 説明 |
|---|---|
| `algebra/` | フェルミオン演算子、スピン、縮約、Wickの定理 |
| `models/` | 正方格子（k点グリッド、分散関係）、Hubbardモデルのパラメータ |
| `diagrams/` | ファインマンダイアグラムのグラフ表現、生成、トポロジー分類、ファインマンルール |
| `symbolic/` | 式木（`G0`, `U`, `Sum`, `Mul`, `Add`, ...）と `Display` による表示 |
| `resummation/` | 粒子-粒子ラダーチャンネル、T行列、Thouless判定条件による Tc の探索 |
| `numerical/` | グリーン関数 G0(k, iwn)、松原周波数の和、粒子-粒子感受率 |
| `visualization/` | Graphviz DOT および D3.js JSON エクスポート |

## 本書の対象読者

本書は2つの読者層を対象としています。**開発者** — `feynman-engine` をライブラリとして使い、ダイアグラムの自動生成や数値的な多体計算を行いたい方。そして **物理学の学生・研究者** — 摂動論、Wickの定理、ラダー再足し上げが実際にどのように実装されるかを理解したい方です。物理の章ではRustの事前知識を前提とせず、コード解説の章では多体理論の事前知識を前提としません。
