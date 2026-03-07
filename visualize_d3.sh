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
echo "Generated HTML:"
echo "  diagram_dot_vertex.html  (dot-vertex style)"
echo "  diagram_wavy.html        (wavy interaction style)"

# PNG export via headless Chrome (optional)
CHROME=""
if [ -x "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" ]; then
  CHROME="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
elif command -v google-chrome &>/dev/null; then
  CHROME="google-chrome"
elif command -v chromium &>/dev/null; then
  CHROME="chromium"
fi

if [ -n "$CHROME" ]; then
  echo ""
  echo "Exporting PNG via headless Chrome..."
  # Convert to absolute file:// URLs for headless Chrome
  ABS_DOT="file://$(pwd)/diagram_dot_vertex.html"
  ABS_WAVY="file://$(pwd)/diagram_wavy.html"

  "$CHROME" --headless --disable-gpu --no-sandbox \
    --screenshot="diagram_dot_vertex.png" --window-size=760,4200 \
    "$ABS_DOT" 2>/dev/null || true
  "$CHROME" --headless --disable-gpu --no-sandbox \
    --screenshot="diagram_wavy.png" --window-size=760,4200 \
    "$ABS_WAVY" 2>/dev/null || true

  if [ -f diagram_dot_vertex.png ] && [ -f diagram_wavy.png ]; then
    echo "  diagram_dot_vertex.png"
    echo "  diagram_wavy.png"
  else
    echo "  (PNG export failed — open HTML in browser and use Export buttons)"
  fi
else
  echo ""
  echo "Chrome not found. Skipping PNG export."
  echo "Open HTML in browser and use Export SVG buttons, or install Chrome."
fi

echo ""
echo "Open in browser:"
echo "  open diagram_dot_vertex.html"
echo "  open diagram_wavy.html"
