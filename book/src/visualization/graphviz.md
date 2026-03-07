# Graphviz Export

The `visualization::dot` module converts Feynman diagrams into
[Graphviz](https://graphviz.org/) DOT format for static rendering to PNG, SVG,
or PDF.

## API

Two public functions are available:

| Function | Input | Output |
|----------|-------|--------|
| `to_dot(diagram)` | a single `&FeynmanDiagram` | DOT string for one diagram |
| `to_dot_all(classified)` | `&[(FeynmanDiagram, i32)]` from `classify_diagrams` | DOT string with one `subgraph cluster` per topology |

### `to_dot` -- single diagram

Produces a `digraph` with:

- **Internal vertices** rendered as filled black circles, labelled by site name.
- **External pseudo-vertices** (`ext_in`, `ext_out`) rendered as plain text
  labels `"in"` and `"out"`.
- **Propagators** as directed edges with spin-dependent styling.

### `to_dot_all` -- classified diagrams

Wraps each topology in a numbered `subgraph cluster_N` with a label showing the
topology index and its combinatorial weight.  Vertex names are prefixed with
`dN_` to avoid collisions across subgraphs.

## Color scheme

| Element | Color | Style |
|---------|-------|-------|
| Spin-up (\\(\uparrow\\)) propagator | blue | solid |
| Spin-down (\\(\downarrow\\)) propagator | red | solid |
| External leg | (spin color) | dashed |

Propagator edge labels use the notation \\(G_0\uparrow\\) or
\\(G_0\downarrow\\).

## CLI usage

The `second_order_self_energy` example supports `--dot` for DOT export:

```bash
# Write DOT to a file
cargo run --example second_order_self_energy -- --dot diagrams.dot

# Write DOT to stdout
cargo run --example second_order_self_energy -- --dot -
```

## Rendering with Graphviz

Once you have a `.dot` file, render it with the `dot` command
(`brew install graphviz` on macOS, `apt install graphviz` on Debian/Ubuntu):

```bash
# PNG output
dot -Tpng diagrams.dot -o diagrams.png

# SVG output (scalable, good for papers)
dot -Tsvg diagrams.dot -o diagrams.svg

# PDF output
dot -Tpdf diagrams.dot -o diagrams.pdf
```

## Quick rendering with `visualize.sh`

The repository includes a convenience script that builds the example, extracts
the DOT output, and renders it to PNG in one step:

```bash
./visualize.sh
# => Generates diagram.dot and diagram.png
```

The script runs:

```bash
cargo run --example second_order_self_energy 2>/dev/null \
  | sed -n '/=== DOT Visualization ===/,$p' \
  | tail -n +2 > diagram.dot

dot -Tpng -o diagram.png diagram.dot
```

## Using the API directly

You can call the DOT export functions from your own Rust code:

```rust
use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::HubbardModel;
use feynman_engine::diagrams::generate::{generate_diagrams, Observable};
use feynman_engine::diagrams::classify::classify_diagrams;
use feynman_engine::visualization::dot::{to_dot, to_dot_all};

fn main() {
    let lattice = SquareLattice::new(4, 4);
    let model = HubbardModel::new(lattice, 1.0, -2.0);

    let diagrams = generate_diagrams(&model, 2, Observable::SelfEnergy);
    let classified = classify_diagrams(&diagrams);

    // Export all classified topologies as a single DOT file
    let dot_all = to_dot_all(&classified);
    std::fs::write("all_topologies.dot", &dot_all).unwrap();

    // Export a single diagram
    if let Some((diag, _weight)) = classified.first() {
        let dot_single = to_dot(diag);
        std::fs::write("single.dot", &dot_single).unwrap();
    }
}
```

## DOT output structure

A typical `to_dot_all` output looks like:

```dot
digraph AllDiagrams {
  compound=true;

  subgraph cluster_0 {
    label="Topology 1 (weight=2)"
    style=rounded;
    node [shape=circle, style=filled, fillcolor=black,
          fontcolor=white, width=0.3];

    d0_v1 [label="v1"];
    d0_v2 [label="v2"];
    d0_ext_in  [shape=none, label="in",  ...];
    d0_ext_out [shape=none, label="out", ...];

    d0_ext_in -> d0_v1 [color=blue, style=dashed];
    d0_v1 -> d0_v2     [color=blue, style=solid];
    d0_v2 -> d0_ext_out [color=blue, style=dashed];
    d0_v1 -> d0_v1     [color=red,  style=solid];
    d0_v2 -> d0_v2     [color=red,  style=solid];
  }

  // ... more subgraphs ...
}
```
