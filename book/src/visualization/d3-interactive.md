# D3.js Interactive Diagrams

For interactive, browser-based visualization, the engine can export diagram data
as JSON and render it with D3.js-powered HTML templates.  Two rendering styles
are provided: a simple dot-vertex view and a textbook-style wavy-line view.

## JSON export

### API

The `visualization::json` module provides:

```rust
pub fn to_json_all(classified: &[(FeynmanDiagram, i32)]) -> String
```

This serializes all classified diagrams into a pretty-printed JSON string via
`serde_json`.

### CLI

The `second_order_self_energy` example supports `--json`:

```bash
# Write JSON to a file
cargo run --example second_order_self_energy -- --json diagrams.json

# Write JSON to stdout
cargo run --example second_order_self_energy -- --json -
```

### JSON schema

The output is an object with a single `diagrams` array.  Each element has:

```json
{
  "diagrams": [
    {
      "topology_id": 1,
      "weight": -1,
      "diagram": {
        "order": 2,
        "vertices": [
          { "id": 1, "site": "v1", "time": "τ1" },
          { "id": 2, "site": "v2", "time": "τ2" }
        ],
        "propagators": [
          { "from": 0, "to": 1, "spin": "Up",   "external": true  },
          { "from": 1, "to": 2, "spin": "Up",   "external": false },
          { "from": 2, "to": 3, "spin": "Up",   "external": true  },
          { "from": 1, "to": 1, "spin": "Down", "external": false },
          { "from": 2, "to": 2, "spin": "Down", "external": false }
        ],
        "sign": -1,
        "symmetry_factor": 1.0
      }
    }
  ]
}
```

| Field | Type | Description |
|-------|------|-------------|
| `topology_id` | `usize` | 1-indexed topology identifier |
| `weight` | `i32` | Combinatorial weight (number of equivalent Wick contractions) |
| `diagram.order` | `usize` | Perturbation order |
| `diagram.vertices` | array | Internal vertices with `id`, `site`, `time` |
| `diagram.propagators` | array | Edges with `from`, `to` (vertex ids), `spin` (`"Up"` / `"Down"`), `external` flag |
| `diagram.sign` | `i32` | Fermion sign from Wick contraction |
| `diagram.symmetry_factor` | `f64` | Symmetry factor |

For external propagators, `from = 0` represents the incoming external point and
`to = max_vertex_id + 1` represents the outgoing external point.

## The `visualize_d3.sh` script

The main entry point for D3 visualization is the `visualize_d3.sh` script.  It
performs the full pipeline in one command:

```bash
./visualize_d3.sh
```

The script:

1. **Builds and runs** the `second_order_self_energy` example.
2. **Extracts JSON** from the example's stdout (after the `=== JSON Data ===`
   marker).
3. **Injects the JSON** into two HTML templates by replacing the
   `__DIAGRAM_JSON__` placeholder, producing:
   - `diagram_dot_vertex.html` -- dot-vertex style
   - `diagram_wavy.html` -- wavy interaction style
4. **(Optional) Exports PNG** via headless Chrome, with ImageMagick whitespace
   trimming.

Open the generated HTML files in a browser:

```bash
open diagram_dot_vertex.html
open diagram_wavy.html
```

## HTML templates

Both templates live in `vis/` and are self-contained single-file HTML
applications using inline D3.js.

### `vis/template_dot_vertex.html` -- Dot-vertex style

Renders diagrams with:

- Black filled circles for interaction vertices
- Straight directed lines for propagators
- Blue lines for spin-up, red lines for spin-down
- Dashed lines for external legs
- Each topology displayed in a card with its weight label

This style is closest to the Graphviz DOT output and is useful for quickly
inspecting diagram topology.

### `vis/template_wavy.html` -- Wavy interaction style

Renders textbook-style Feynman diagrams with:

- Wavy lines for interaction vertices (the Hubbard \\(U\\))
- Solid directed lines for fermion propagators
- Momentum-space (\\(k\\)-space) labelling
- Blue for spin-up, red for spin-down
- External legs shown as dashed lines entering/leaving the diagram

This style produces publication-quality diagrams similar to those found in
condensed matter physics textbooks.

### SVG export buttons

Both templates include an export bar at the top with two buttons:

- **Export All SVG** -- downloads each diagram as a separate `.svg` file.
- **Export Combined SVG** -- downloads a single `.svg` containing all diagrams
  stacked vertically.

These work entirely in the browser with no server-side dependencies.

## PNG export via headless Chrome

The `visualize_d3.sh` script automatically attempts PNG export if Chrome (or
Chromium) is detected:

```bash
# The script searches for Chrome in these locations:
# - /Applications/Google Chrome.app/Contents/MacOS/Google Chrome  (macOS)
# - google-chrome  (Linux, in PATH)
# - chromium       (Linux, in PATH)
```

The headless Chrome rendering uses `--screenshot` with a fixed window size.  If
[ImageMagick](https://imagemagick.org/) is available (`magick` or `convert`),
the script automatically trims whitespace from the resulting PNG:

```bash
magick diagram_wavy.png -trim +repage diagram_wavy.png
```

If Chrome is not installed, the script prints a message suggesting you open the
HTML files in a browser and use the built-in SVG export buttons instead.

## Using the JSON API directly

```rust
use feynman_engine::models::lattice::SquareLattice;
use feynman_engine::models::hubbard::HubbardModel;
use feynman_engine::diagrams::generate::{generate_diagrams, Observable};
use feynman_engine::diagrams::classify::classify_diagrams;
use feynman_engine::visualization::json::to_json_all;

fn main() {
    let lattice = SquareLattice::new(4, 4);
    let model = HubbardModel::new(lattice, 1.0, -2.0);

    let diagrams = generate_diagrams(&model, 2, Observable::SelfEnergy);
    let classified = classify_diagrams(&diagrams);

    let json = to_json_all(&classified);
    std::fs::write("diagrams.json", &json).unwrap();

    println!("Wrote {} topologies to diagrams.json", classified.len());
}
```

## Workflow summary

| Step | Command | Output |
|------|---------|--------|
| Generate JSON | `cargo run --example second_order_self_energy -- --json diagrams.json` | `diagrams.json` |
| One-step HTML + PNG | `./visualize_d3.sh` | `diagram_dot_vertex.html`, `diagram_wavy.html`, `.png` files |
| Open in browser | `open diagram_wavy.html` | Interactive D3.js visualization |
| Export SVG | Click "Export All SVG" button in browser | `.svg` file downloads |
