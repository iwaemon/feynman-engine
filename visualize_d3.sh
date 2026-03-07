#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "Building and running example..."
OUTPUT=$(cargo run --example second_order_self_energy 2>/dev/null)

# Extract JSON: everything after "=== JSON Data ===" line
JSON=$(echo "$OUTPUT" | awk '/^=== JSON Data ===$/{found=1; next} found{print}')

if [ -z "$JSON" ]; then
  echo "ERROR: No JSON data found in example output" >&2
  exit 1
fi

# Use Python to safely inject JSON into templates (handles special characters)
echo "Generating dot-vertex HTML..."
python3 -c "
import sys
template = open('vis/template_dot_vertex.html').read()
json_data = sys.stdin.read()
print(template.replace('__DIAGRAM_JSON__', json_data))
" <<< "$JSON" > diagram_dot_vertex.html

echo "Generating wavy HTML..."
python3 -c "
import sys
template = open('vis/template_wavy.html').read()
json_data = sys.stdin.read()
print(template.replace('__DIAGRAM_JSON__', json_data))
" <<< "$JSON" > diagram_wavy.html

echo ""
echo "Generated:"
echo "  diagram_dot_vertex.html  (dot-vertex style)"
echo "  diagram_wavy.html        (wavy interaction style)"
echo ""
echo "Open in browser:"
echo "  open diagram_dot_vertex.html"
echo "  open diagram_wavy.html"
