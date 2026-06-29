#!/bin/bash
# Run the full scrape → featurize → replace_bfactors pipeline for a given config.
# Usage: ./run_pipeline.sh [config.toml]

set -euo pipefail

CONFIG=${1:-config.toml}

echo "=== scrape.py ==="
pixi run python scripts/scrape.py "$CONFIG"

echo ""
echo "=== featurize.py ==="
pixi run python scripts/featurize.py "$CONFIG"

echo ""
echo "=== replace_bfactors.py ==="
pixi run python scripts/replace_bfactors.py "$CONFIG"

echo ""
echo "Pipeline complete."
