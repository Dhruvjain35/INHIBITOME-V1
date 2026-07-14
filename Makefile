.PHONY: setup token pilot test lint clean gate

setup:            ## Install deps into a uv-managed venv
	uv sync --extra dev

token:            ## One-time CAVE auth setup + verify
	uv run python scripts/00_setup_cave_token.py

pilot:            ## Run the full 10-day pilot in order (stops on a failed gate)
	uv run python scripts/01_freeze_and_join.py
	uv run python scripts/02_sample_accounting.py
	uv run python scripts/03_functional_targets.py
	uv run python scripts/04_fingerprints.py
	uv run python scripts/05_blinded_comparison.py

gate:             ## Re-run just the Day-3 data gate
	uv run python scripts/02_sample_accounting.py

test:             ## Run unit tests (pure logic; no network)
	uv run pytest -q

lint:
	uv run ruff check src scripts tests

clean:            ## Drop the CAVE query cache (forces re-download)
	rm -rf data/cache

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
	  awk 'BEGIN {FS = ":.*?## "}; {printf "  %-10s %s\n", $$1, $$2}'
