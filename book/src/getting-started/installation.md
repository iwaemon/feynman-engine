# Installation

## Rust Toolchain

feynman-engine requires Rust edition 2021 with a minimum supported Rust version (MSRV) of **1.70**.

If you do not have Rust installed, follow the official instructions at
[rustup.rs](https://rustup.rs/):

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Verify your toolchain version:

```bash
rustc --version   # must be >= 1.70
```

## Clone and Build

```bash
git clone https://github.com/shumpei/feynman-engine.git
cd feynman-engine
cargo build
cargo test
```

All tests should pass in under a second. A successful `cargo test` confirms that
the library, examples, and integration tests are working correctly.

## Dependencies

The following crates are pulled in automatically by Cargo:

| Crate | Purpose |
|---|---|
| `num-complex` | Complex number arithmetic for Matsubara Green's functions |
| `rayon` | Data parallelism for momentum and frequency sums |
| `itertools` | Combinatorial utilities used in Wick contractions |
| `serde` / `serde_json` | JSON serialization for diagram export |

Dev-only (used in tests):

| Crate | Purpose |
|---|---|
| `approx` | Floating-point approximate equality assertions |

## Optional External Tools

These are **not** required to build or run the library, but enable additional
visualization workflows:

### Graphviz

Renders `.dot` diagram files to PNG or SVG.

```bash
# macOS
brew install graphviz

# Ubuntu / Debian
apt install graphviz
```

Usage:

```bash
cargo run --example second_order_self_energy -- --dot diagrams.dot
dot -Tpng diagrams.dot -o diagrams.png
```

### Google Chrome (headless)

Used for headless PNG export of the D3.js HTML visualization. Any
Chromium-based browser works. This is only needed if you want to generate
static PNG images from the interactive HTML diagrams.

### ImageMagick

Used to auto-trim whitespace from exported PNG images:

```bash
# macOS
brew install imagemagick

# Ubuntu / Debian
apt install imagemagick
```
