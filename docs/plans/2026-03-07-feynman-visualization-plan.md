# Feynman Diagram D3.js Visualization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add physics-style Feynman diagram visualization using D3.js (dot-vertex and wavy versions) alongside the existing Graphviz output.

**Architecture:** Rust example outputs JSON (diagram graph structure) in addition to existing DOT. D3.js standalone HTML reads the JSON and renders SVG with Bezier curves, arrows, and proper Feynman diagram layout. Shell script orchestrates the pipeline.

**Tech Stack:** Rust + serde_json (JSON output), D3.js v7 (CDN, SVG rendering), Bash (orchestration)

---

### Task 1: Add serde dependency and derive Serialize

**Files:**
- Modify: `Cargo.toml`
- Modify: `src/algebra/operators.rs:1-5` (add Serialize derive to Spin)
- Modify: `src/diagrams/graph.rs:1-25` (add Serialize derive to Vertex, Propagator, FeynmanDiagram)

**Step 1: Add serde and serde_json to Cargo.toml**

```toml
[dependencies]
num-complex = "0.4"
rayon = "1.10"
itertools = "0.13"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
```

**Step 2: Add Serialize derive to Spin**

In `src/algebra/operators.rs`, change:
```rust
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Spin {
```
to:
```rust
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize)]
pub enum Spin {
```

**Step 3: Add Serialize derive to graph structs**

In `src/diagrams/graph.rs`, add `serde::Serialize` derive to `Vertex`, `Propagator`, `FeynmanDiagram`:
```rust
#[derive(Debug, Clone, serde::Serialize)]
pub struct Vertex {
    pub id: usize,
    pub site: String,
    pub time: String,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct Propagator {
    pub from: usize,
    pub to: usize,
    pub spin: Spin,
    pub external: bool,
}

#[derive(Debug, Clone, serde::Serialize)]
pub struct FeynmanDiagram {
    pub order: usize,
    pub vertices: Vec<Vertex>,
    pub propagators: Vec<Propagator>,
    pub sign: i32,
    pub symmetry_factor: f64,
}
```

**Step 4: Run tests to verify nothing breaks**

Run: `cargo test`
Expected: All existing tests PASS

**Step 5: Commit**

```bash
git add Cargo.toml src/algebra/operators.rs src/diagrams/graph.rs
git commit -m "feat: add serde Serialize derive to diagram structs"
```

---

### Task 2: Create JSON serialization module

**Files:**
- Create: `src/visualization/json.rs`
- Modify: `src/visualization/mod.rs`

**Step 1: Write the test for to_json_all**

In `src/visualization/json.rs`:
```rust
use crate::diagrams::graph::{FeynmanDiagram, Vertex, Propagator};
use crate::algebra::operators::Spin;

/// Serializable wrapper for classified diagrams
#[derive(serde::Serialize)]
pub struct ClassifiedDiagramSet {
    pub diagrams: Vec<ClassifiedDiagram>,
}

#[derive(serde::Serialize)]
pub struct ClassifiedDiagram {
    pub topology_id: usize,
    pub weight: i32,
    pub diagram: FeynmanDiagram,
}

/// Convert classified diagrams to JSON string
pub fn to_json_all(classified: &[(FeynmanDiagram, i32)]) -> String {
    let set = ClassifiedDiagramSet {
        diagrams: classified
            .iter()
            .enumerate()
            .map(|(i, (diag, weight))| ClassifiedDiagram {
                topology_id: i + 1,
                weight: *weight,
                diagram: diag.clone(),
            })
            .collect(),
    };
    serde_json::to_string_pretty(&set).expect("Failed to serialize diagrams to JSON")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_to_json_all_basic() {
        let d = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex {
                id: 1,
                site: "v1".into(),
                time: "τ1".into(),
            }],
            propagators: vec![
                Propagator { from: 0, to: 1, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 2, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 1, spin: Spin::Down, external: false },
            ],
            sign: -1,
            symmetry_factor: 1.0,
        };

        let classified = vec![(d, -1)];
        let json = to_json_all(&classified);

        assert!(json.contains("\"topology_id\": 1"));
        assert!(json.contains("\"weight\": -1"));
        assert!(json.contains("\"spin\": \"Up\""));
        assert!(json.contains("\"spin\": \"Down\""));
        assert!(json.contains("\"external\": true"));
        assert!(json.contains("\"external\": false"));
    }

    #[test]
    fn test_to_json_all_multiple() {
        let d1 = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex { id: 1, site: "v1".into(), time: "τ1".into() }],
            propagators: vec![
                Propagator { from: 1, to: 1, spin: Spin::Up, external: false },
            ],
            sign: 1,
            symmetry_factor: 1.0,
        };
        let d2 = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex { id: 1, site: "v1".into(), time: "τ1".into() }],
            propagators: vec![
                Propagator { from: 0, to: 1, spin: Spin::Up, external: true },
                Propagator { from: 1, to: 2, spin: Spin::Up, external: true },
            ],
            sign: -1,
            symmetry_factor: 1.0,
        };

        let classified = vec![(d1, 1), (d2, -1)];
        let json = to_json_all(&classified);
        let parsed: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed["diagrams"].as_array().unwrap().len(), 2);
        assert_eq!(parsed["diagrams"][0]["topology_id"], 1);
        assert_eq!(parsed["diagrams"][1]["topology_id"], 2);
    }
}
```

**Step 2: Register the module**

In `src/visualization/mod.rs`:
```rust
pub mod dot;
pub mod json;
```

**Step 3: Run tests**

Run: `cargo test visualization::json`
Expected: Both tests PASS

**Step 4: Commit**

```bash
git add src/visualization/json.rs src/visualization/mod.rs
git commit -m "feat: add JSON serialization for classified diagrams"
```

---

### Task 3: Extend example to output JSON

**Files:**
- Modify: `examples/second_order_self_energy.rs`

**Step 1: Add JSON output after DOT output**

Add `use feynman_engine::visualization::json::to_json_all;` at imports, then after the DOT section add:

```rust
    // 6. Output JSON for D3.js visualization
    let json = to_json_all(&all_classified);
    println!("=== JSON Data ===");
    println!("{}", json);
```

**Step 2: Run the example to verify JSON output**

Run: `cargo run --example second_order_self_energy 2>/dev/null | tail -20`
Expected: JSON output with diagrams array containing 9 topologies

**Step 3: Commit**

```bash
git add examples/second_order_self_energy.rs
git commit -m "feat: add JSON output to self-energy example"
```

---

### Task 4: Create D3.js dot-vertex template

**Files:**
- Create: `vis/template_dot_vertex.html`

This is the main D3.js visualization. Key rendering logic:
- Left-to-right layout with external legs as horizontal lines
- Interaction vertices as black dots (●)
- Propagators as Bezier curves with arrows
- Self-loops (tadpoles) as teardrop shapes
- Spin colors: blue=↑, red=↓
- Momentum/frequency labels on each propagator line

**Step 1: Create the HTML template**

Create `vis/template_dot_vertex.html`:

```html
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Feynman Diagrams — Dot Vertex Style</title>
<style>
  body { font-family: 'Helvetica Neue', Arial, sans-serif; background: #fff; margin: 20px; }
  h1 { font-size: 18px; color: #333; }
  .diagram-container { display: flex; flex-wrap: wrap; gap: 20px; }
  .diagram-card { border: 1px solid #ccc; border-radius: 8px; padding: 10px; background: #fafafa; }
  .diagram-card h3 { margin: 0 0 8px 0; font-size: 13px; color: #555; }
  .propagator-up { stroke: #2266cc; fill: none; }
  .propagator-down { stroke: #cc2222; fill: none; }
  .propagator-up-dashed { stroke: #2266cc; fill: none; stroke-dasharray: 6,3; }
  .propagator-down-dashed { stroke: #cc2222; fill: none; stroke-dasharray: 6,3; }
  .vertex { fill: #222; }
  .label { font-size: 11px; fill: #333; }
  .ext-label { font-size: 12px; fill: #666; font-style: italic; }
</style>
</head>
<body>
<h1>Self-Energy Feynman Diagrams (Dot Vertex)</h1>
<div class="diagram-container" id="container"></div>

<script type="module">
import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

// DATA_PLACEHOLDER will be replaced by the shell script
const data = __DIAGRAM_JSON__;

const CARD_W = 260;
const CARD_H = 180;
const MARGIN = { top: 30, right: 30, bottom: 20, left: 30 };
const W = CARD_W - MARGIN.left - MARGIN.right;
const H = CARD_H - MARGIN.top - MARGIN.bottom;
const VERTEX_R = 5;

function layoutDiagram(diag) {
  const verts = diag.diagram.vertices;
  const props = diag.diagram.propagators;
  const n = verts.length;

  // Position vertices evenly along horizontal axis
  const positions = {};
  verts.forEach((v, i) => {
    positions[v.id] = {
      x: n === 1 ? W / 2 : (W * i) / (n - 1),
      y: H / 2
    };
  });

  // External leg positions
  const extOutId = Math.max(...verts.map(v => v.id)) + 1;
  positions[0] = { x: -MARGIN.left + 10, y: H / 2 };         // ext_in (id=0)
  positions[extOutId] = { x: W + MARGIN.right - 10, y: H / 2 }; // ext_out

  return { positions, extOutId };
}

function drawArrowMarker(svg, id, color) {
  svg.append("defs").append("marker")
    .attr("id", id)
    .attr("viewBox", "0 -4 8 8")
    .attr("refX", 7)
    .attr("refY", 0)
    .attr("markerWidth", 6)
    .attr("markerHeight", 6)
    .attr("orient", "auto")
    .append("path")
    .attr("d", "M0,-3L7,0L0,3Z")
    .attr("fill", color);
}

function drawPropagator(g, p, positions, extOutId, propIndex) {
  const color = p.spin === "Up" ? "#2266cc" : "#cc2222";
  const spinLabel = p.spin === "Up" ? "↑" : "↓";
  const isDashed = p.external;

  const fromId = (p.external && p.from === 0) ? 0 : p.from;
  const toId = (p.external && p.to === extOutId) ? extOutId : p.to;
  const from = positions[fromId];
  const to = positions[toId];

  if (!from || !to) return;

  const markerId = `arrow-${color.replace('#','')}-${propIndex}`;
  drawArrowMarker(g, markerId, color);

  if (fromId === toId) {
    // Self-loop (tadpole): teardrop shape
    const cx = from.x;
    const cy = from.y;
    const loopR = 22;
    // Count how many self-loops we've already drawn on this vertex for offset
    const path = d3.path();
    path.moveTo(cx, cy - VERTEX_R);
    path.bezierCurveTo(cx - loopR, cy - loopR * 1.8, cx + loopR, cy - loopR * 1.8, cx, cy - VERTEX_R);

    g.append("path")
      .attr("d", path.toString())
      .attr("stroke", color)
      .attr("stroke-width", 1.8)
      .attr("fill", "none")
      .attr("marker-mid", `url(#${markerId})`);

    // Label
    g.append("text")
      .attr("class", "label")
      .attr("x", cx)
      .attr("y", cy - loopR * 1.6)
      .attr("text-anchor", "middle")
      .attr("fill", color)
      .text(`k${propIndex+1},iω${propIndex+1} ${spinLabel}`);
  } else {
    // Curved propagator between different vertices
    const dx = to.x - from.x;
    const dy = to.y - from.y;
    const dist = Math.sqrt(dx * dx + dy * dy);

    // Offset for multiple edges between same pair
    const curveOffset = 20 * (propIndex % 2 === 0 ? 1 : -1);
    const mx = (from.x + to.x) / 2;
    const my = (from.y + to.y) / 2 + curveOffset;

    // Bezier control points for smooth curve
    const path = d3.path();
    path.moveTo(from.x, from.y);
    path.quadraticCurveTo(mx, my, to.x, to.y);

    g.append("path")
      .attr("d", path.toString())
      .attr("stroke", color)
      .attr("stroke-width", 1.8)
      .attr("fill", "none")
      .attr("stroke-dasharray", isDashed ? "6,3" : null)
      .attr("marker-end", `url(#${markerId})`);

    // Label at midpoint
    const labelOffset = curveOffset > 0 ? -8 : 14;
    if (!p.external) {
      g.append("text")
        .attr("class", "label")
        .attr("x", mx)
        .attr("y", my + labelOffset)
        .attr("text-anchor", "middle")
        .attr("fill", color)
        .text(`k${propIndex+1},iω${propIndex+1} ${spinLabel}`);
    }
  }
}

// Render each topology
const container = d3.select("#container");

data.diagrams.forEach((diag, idx) => {
  const card = container.append("div").attr("class", "diagram-card");
  card.append("h3").text(`Topology ${diag.topology_id} (weight = ${diag.weight})`);

  const svg = card.append("svg")
    .attr("width", CARD_W)
    .attr("height", CARD_H);

  const g = svg.append("g")
    .attr("transform", `translate(${MARGIN.left},${MARGIN.top})`);

  const { positions, extOutId } = layoutDiagram(diag);

  // Draw propagators
  diag.diagram.propagators.forEach((p, pi) => {
    drawPropagator(g, p, positions, extOutId, pi);
  });

  // Draw vertices
  diag.diagram.vertices.forEach(v => {
    const pos = positions[v.id];
    g.append("circle")
      .attr("class", "vertex")
      .attr("cx", pos.x)
      .attr("cy", pos.y)
      .attr("r", VERTEX_R);
  });

  // External labels
  g.append("text").attr("class", "ext-label")
    .attr("x", positions[0].x - 5).attr("y", positions[0].y + 4)
    .attr("text-anchor", "end").text("in");
  g.append("text").attr("class", "ext-label")
    .attr("x", positions[extOutId].x + 5).attr("y", positions[extOutId].y + 4)
    .attr("text-anchor", "start").text("out");
});
</script>
</body>
</html>
```

**Step 2: Verify the HTML is valid**

Run: `head -5 vis/template_dot_vertex.html`
Expected: Shows DOCTYPE and head

**Step 3: Commit**

```bash
git add vis/template_dot_vertex.html
git commit -m "feat: add D3.js dot-vertex Feynman diagram template"
```

---

### Task 5: Create D3.js wavy template

**Files:**
- Create: `vis/template_wavy.html`

This is a copy of the dot-vertex template with one key change: interaction vertices are connected by a short wavy line instead of just a dot. The wavy line is rendered as a sine-wave SVG path between connecting propagators at each vertex.

**Step 1: Create the wavy template**

Copy `vis/template_dot_vertex.html` to `vis/template_wavy.html` and modify:

1. Change title to "Feynman Diagrams — Wavy Interaction Style"
2. Change h1 text
3. Add a `drawWavyLine` function that draws a sine-wave path vertically from each vertex:

```javascript
function drawWavyLine(g, cx, cy, length, amplitude, periods) {
  const points = [];
  const steps = periods * 20;
  for (let i = 0; i <= steps; i++) {
    const t = i / steps;
    const x = cx + amplitude * Math.sin(2 * Math.PI * periods * t);
    const y = cy - VERTEX_R - t * length;
    points.push([x, y]);
  }
  const line = d3.line().curve(d3.curveBasis);
  g.append("path")
    .attr("d", line(points))
    .attr("stroke", "#666")
    .attr("stroke-width", 1.5)
    .attr("fill", "none");
}
```

4. After drawing each vertex, call `drawWavyLine(g, pos.x, pos.y, 25, 4, 3)` to draw a wavy interaction stub.

**Step 2: Commit**

```bash
git add vis/template_wavy.html
git commit -m "feat: add D3.js wavy interaction Feynman diagram template"
```

---

### Task 6: Create visualize_d3.sh script

**Files:**
- Create: `visualize_d3.sh`

**Step 1: Write the shell script**

```bash
#!/usr/bin/env bash
set -euo pipefail

# 1. Run the Rust example and extract JSON
echo "Generating diagram data..."
JSON=$(cargo run --example second_order_self_energy 2>/dev/null | sed -n '/=== JSON Data ===/,/^$/{ /=== JSON Data ===/d; p; }')

# If sed didn't capture properly, try alternate extraction (everything after JSON marker to end)
if [ -z "$JSON" ]; then
  JSON=$(cargo run --example second_order_self_energy 2>/dev/null | awk '/=== JSON Data ===/{found=1; next} found{print}')
fi

# 2. Generate dot-vertex HTML
echo "Generating dot-vertex HTML..."
sed "s|__DIAGRAM_JSON__|${JSON//$'\n'/\\n}|" vis/template_dot_vertex.html > diagram_dot_vertex.html 2>/dev/null || \
  python3 -c "
import sys
template = open('vis/template_dot_vertex.html').read()
json_data = '''$JSON'''
print(template.replace('__DIAGRAM_JSON__', json_data))
" > diagram_dot_vertex.html

# 3. Generate wavy HTML
echo "Generating wavy HTML..."
python3 -c "
import sys
template = open('vis/template_wavy.html').read()
json_data = '''$JSON'''
print(template.replace('__DIAGRAM_JSON__', json_data))
" > diagram_wavy.html

echo "Generated:"
echo "  diagram_dot_vertex.html  (open in browser)"
echo "  diagram_wavy.html        (open in browser)"

# 4. Optional: generate static SVG if Chrome/Chromium is available
if command -v "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" &>/dev/null; then
  echo "Exporting SVG via headless Chrome..."
  # Note: SVG export requires additional JS in template; skip for now
  echo "  (SVG export: not yet implemented)"
fi
```

**Step 2: Make executable**

Run: `chmod +x visualize_d3.sh`

**Step 3: Run end-to-end test**

Run: `./visualize_d3.sh`
Expected: Creates `diagram_dot_vertex.html` and `diagram_wavy.html`

**Step 4: Verify HTML files open in browser**

Run: `open diagram_dot_vertex.html` (macOS)
Expected: Browser shows Feynman diagrams with Bezier curves, arrows, labels

**Step 5: Commit**

```bash
git add visualize_d3.sh
git commit -m "feat: add D3.js visualization pipeline script"
```

---

### Task 7: Polish and adjust layout

After seeing the initial render, iterate on:

**Step 1: Adjust self-loop rendering**

- Multiple self-loops on the same vertex need vertical offset (one above, one below)
- Each self-loop at a vertex should get a different angle/direction

**Step 2: Adjust propagator curve offsets**

- When multiple propagators connect the same vertex pair, they need distinct curvature
- Use alternating positive/negative offsets with increasing magnitude

**Step 3: Adjust label placement**

- Labels should not overlap with curves
- External leg labels ("in", "out") should be clear

**Step 4: Run and verify**

Run: `./visualize_d3.sh && open diagram_dot_vertex.html`
Expected: All 9 topologies render clearly without overlapping elements

**Step 5: Commit**

```bash
git add vis/template_dot_vertex.html vis/template_wavy.html
git commit -m "fix: polish diagram layout and label placement"
```

---

### Task 8: Add SVG export support

**Files:**
- Modify: `vis/template_dot_vertex.html`
- Modify: `vis/template_wavy.html`

**Step 1: Add "Export SVG" button to each template**

Add a button and JS function that serializes the SVG elements to standalone SVG files:

```javascript
function exportSVG() {
  const svgs = document.querySelectorAll('svg');
  svgs.forEach((svg, i) => {
    const serializer = new XMLSerializer();
    const svgStr = serializer.serializeToString(svg);
    const blob = new Blob([svgStr], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `diagram_topology_${i+1}.svg`;
    a.click();
  });
}
```

**Step 2: Commit**

```bash
git add vis/template_dot_vertex.html vis/template_wavy.html
git commit -m "feat: add SVG export button to D3.js templates"
```
