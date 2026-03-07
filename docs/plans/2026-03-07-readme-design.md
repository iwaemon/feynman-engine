# README / 使い方ドキュメント設計

## 目的

物理学研究者向けに、feynman-engineの概要・使い方・背景を提供する。

## ファイル構成

- `README.md`（英語）— GitHubトップに表示
- `docs/usage-ja.md`（日本語）— 日本語の使い方ガイド

## README.md セクション構成

1. **タイトル + サブタイトル** — `feynman-engine`: A perturbation theory engine for the Hubbard model
2. **Overview** — 1段落。Wick定理→ダイアグラム生成→記号式→梯子型再総和→Tc探索のパイプライン
3. **Features** — 箇条書き（Wick定理、ダイアグラム生成・分類、記号式、PP梯子型再総和、Graphviz DOT出力）
4. **Installation** — `cargo build` / `cargo test`
5. **Quick Start** — 2つのexample（`ladder_resummation`, `second_order_self_energy`）の実行コマンドと出力例
6. **Architecture** — データフロー図（CLAUDE.mdベース）とモジュール一覧
7. **Physics Background** — Hubbard模型、Matsubara形式、Thouless判定条件の簡潔な説明
8. **License** — 省略またはTBD

## docs/usage-ja.md セクション構成

README.mdの日本語版に加え:
- exampleの出力の物理的意味の解説
- 各モジュールの使い方の簡単な説明

## 設計方針

- 想定読者: 摂動論やHubbard模型に馴染みのある物理学研究者・学生
- exampleの実際の出力を含める
- YAGNI: API リファレンス、Contributing ガイド、ベンチマークは現時点では含めない
