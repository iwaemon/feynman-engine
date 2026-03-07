# ファインマンダイアグラム可視化の改善 設計書

## 概要

既存のGraphviz DOT出力はそのまま残し、D3.jsを用いた物理的に正しいファインマンダイアグラム可視化を追加する。

## バージョン構成

| バージョン | 相互作用の描き方 | 出力形式 | ツール |
|-----------|----------------|---------|-------|
| original (既存) | Graphviz DOTスタイル | PNG | `visualize.sh` (変更なし) |
| dot-vertex | 点のみ（接触相互作用） | HTML + SVG | D3.js (CDN) |
| wavy | 短い波線 | HTML + SVG | D3.js (CDN) |

## 共通仕様

- レイアウト: 左→右の横配置。外線in(左端)→out(右端)が水平に流れる
- スピン色: 青=↑、赤=↓
- 伝搬線: ベジェ曲線、矢印付き
- 自己ループ（タドポール）: 涙型ループ（ベジェ曲線）
- ラベル: 各伝搬線にk, iωの運動量・周波数ラベル
- 頂点: 黒い点（●）

## アーキテクチャ

```
Rustバイナリ (既存example)
  → DOT出力 (既存、変更なし)
  → JSON出力 (新規: ダイアグラムのグラフ構造をJSON化)

visualize.sh        → diagram.png (既存、変更なし)
visualize_d3.sh     → diagram_dot_vertex.html, diagram_wavy.html (新規)
                    → diagram_dot_vertex.svg, diagram_wavy.svg (静的エクスポート)
```

### データフロー

1. Rustのexampleを拡張してJSON形式でもダイアグラムデータを出力
2. D3.jsのHTMLテンプレートがJSONを読み込んで描画
3. 頂点の座標は左→右レイアウトで固定配置

## D3.js HTML仕様

- スタンドアロンHTML（CDNから`d3@7`を読み込み）
- SVGでベジェ曲線の伝搬線を描画
- `d3-shape`のcurve関数を活用
- マーカー要素で矢印を定義
- 波線は`d3.line`でsin波パスを生成（wavyのみ）
- 各トポロジーをサブグラフとして横に並べる

## ファイル構成（新規）

```
src/visualization/json.rs   — ダイアグラムのJSON出力
vis/template_dot_vertex.html — D3.js テンプレート（点頂点版）
vis/template_wavy.html       — D3.js テンプレート（波線版）
visualize_d3.sh              — JSON生成→HTMLに埋め込み→SVGエクスポート
```

## 決定事項

- 既存のvisualizer.sh/diagram.pngは変更しない
- 相互作用: 点版と波線版の2バージョン
- スピン: 色で区別（青=↑、赤=↓）
- 出力: インタラクティブHTML + 静的SVG/PNG
- レイアウト: 左→右の横配置
- ラベル: 運動量・周波数あり
- 曲線: ベジェ曲線（滑らか）、タドポールは涙型
