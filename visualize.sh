#!/usr/bin/env bash
# 例の出力から DOT 部分だけを取り出して diagram.dot に保存
# DOT は「=== DOT Visualization ===」から末尾まで（空行で止めると不完全になる）
cargo run --example second_order_self_energy 2>/dev/null | sed -n '/=== DOT Visualization ===/,$p' | tail -n +2 > diagram.dot

# Graphviz で PNG に変換（dot は brew install graphviz でインストール）
dot -Tpng -o diagram.png diagram.dot

echo "Generated diagram.png"
