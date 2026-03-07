use crate::algebra::operators::Spin;
use crate::diagrams::graph::FeynmanDiagram;

/// Convert a Feynman diagram to Graphviz DOT format
pub fn to_dot(diagram: &FeynmanDiagram) -> String {
    let mut dot = String::new();
    dot.push_str("digraph FeynmanDiagram {\n");
    dot.push_str("  rankdir=LR;\n");
    dot.push_str("  node [shape=circle, style=filled, fillcolor=black, fontcolor=white, width=0.3];\n");
    dot.push_str("  edge [arrowsize=0.8];\n\n");

    // External pseudo-vertices
    let has_external = diagram.propagators.iter().any(|p| p.external);
    if has_external {
        dot.push_str("  ext_in [shape=none, label=\"in\", fillcolor=white, fontcolor=black];\n");
        dot.push_str("  ext_out [shape=none, label=\"out\", fillcolor=white, fontcolor=black];\n");
    }

    // Internal vertices
    for v in &diagram.vertices {
        dot.push_str(&format!("  v{} [label=\"{}\"];\n", v.id, v.site));
    }
    dot.push('\n');

    // Propagators
    for p in &diagram.propagators {
        let color = match p.spin {
            Spin::Up => "blue",
            Spin::Down => "red",
        };
        let style = if p.external { "dashed" } else { "solid" };
        let spin_label = match p.spin {
            Spin::Up => "\u{2191}",
            Spin::Down => "\u{2193}",
        };

        let from_name = if p.external && p.from == 0 {
            "ext_in".to_string()
        } else {
            format!("v{}", p.from)
        };

        // Determine ext_out id: it's the max vertex id + 1
        let ext_out_id = diagram.vertices.iter().map(|v| v.id).max().unwrap_or(0) + 1;
        let to_name = if p.external && p.to == ext_out_id {
            "ext_out".to_string()
        } else {
            format!("v{}", p.to)
        };

        dot.push_str(&format!(
            "  {} -> {} [color={}, style={}, label=\"G\u{2080}{}\"];\n",
            from_name, to_name, color, style, spin_label
        ));
    }

    dot.push_str("}\n");
    dot
}

/// Convert multiple classified diagrams to DOT with subgraphs
pub fn to_dot_all(classified: &[(FeynmanDiagram, i32)]) -> String {
    let mut dot = String::new();
    dot.push_str("digraph AllDiagrams {\n");
    dot.push_str("  compound=true;\n\n");

    for (idx, (diag, weight)) in classified.iter().enumerate() {
        dot.push_str(&format!("  subgraph cluster_{} {{\n", idx));
        dot.push_str(&format!(
            "    label=\"Topology {} (weight={})\"\n",
            idx + 1,
            weight
        ));
        dot.push_str("    style=rounded;\n");
        dot.push_str("    node [shape=circle, style=filled, fillcolor=black, fontcolor=white, width=0.3];\n\n");

        let prefix = format!("d{}_", idx);

        for v in &diag.vertices {
            dot.push_str(&format!("    {prefix}v{} [label=\"{}\"];\n", v.id, v.site));
        }

        let has_ext = diag.propagators.iter().any(|p| p.external);
        if has_ext {
            dot.push_str(&format!(
                "    {prefix}ext_in [shape=none, label=\"in\", fillcolor=white, fontcolor=black];\n"
            ));
            dot.push_str(&format!(
                "    {prefix}ext_out [shape=none, label=\"out\", fillcolor=white, fontcolor=black];\n"
            ));
        }

        let ext_out_id = diag.vertices.iter().map(|v| v.id).max().unwrap_or(0) + 1;

        for p in &diag.propagators {
            let color = match p.spin {
                Spin::Up => "blue",
                Spin::Down => "red",
            };
            let style = if p.external { "dashed" } else { "solid" };

            let from_name = if p.external && p.from == 0 {
                format!("{prefix}ext_in")
            } else {
                format!("{prefix}v{}", p.from)
            };

            let to_name = if p.external && p.to == ext_out_id {
                format!("{prefix}ext_out")
            } else {
                format!("{prefix}v{}", p.to)
            };

            dot.push_str(&format!(
                "    {} -> {} [color={}, style={}];\n",
                from_name, to_name, color, style
            ));
        }

        dot.push_str("  }\n\n");
    }

    dot.push_str("}\n");
    dot
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagrams::graph::*;

    #[test]
    fn test_to_dot_basic() {
        let d = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex {
                id: 1,
                site: "v1".into(),
                time: "\u{03c4}1".into(),
            }],
            propagators: vec![
                Propagator {
                    from: 0,
                    to: 1,
                    spin: Spin::Up,
                    external: true,
                },
                Propagator {
                    from: 1,
                    to: 2,
                    spin: Spin::Up,
                    external: true,
                },
                Propagator {
                    from: 1,
                    to: 1,
                    spin: Spin::Down,
                    external: false,
                },
            ],
            sign: 1,
            symmetry_factor: 1.0,
        };

        let dot = to_dot(&d);
        assert!(dot.contains("digraph"));
        assert!(dot.contains("v1"));
        assert!(dot.contains("ext_in"));
        assert!(dot.contains("ext_out"));
        assert!(dot.contains("blue"));
        assert!(dot.contains("red"));
    }

    #[test]
    fn test_to_dot_all() {
        let d1 = FeynmanDiagram {
            order: 1,
            vertices: vec![Vertex {
                id: 1,
                site: "v1".into(),
                time: "\u{03c4}1".into(),
            }],
            propagators: vec![Propagator {
                from: 1,
                to: 1,
                spin: Spin::Up,
                external: false,
            }],
            sign: 1,
            symmetry_factor: 1.0,
        };

        let classified = vec![(d1, 1)];
        let dot = to_dot_all(&classified);
        assert!(dot.contains("subgraph cluster_0"));
        assert!(dot.contains("weight=1"));
    }
}
