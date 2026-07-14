"""Load the frozen pilot config (config/pilot.yaml).

Single source of truth for every reproducibility-critical string: materialization version, table
names, gates, seeds. Import `CFG` everywhere; never hard-code a table name in a module.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

# Repo root = two levels up from this file (src/inhibitome/config.py).
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "config" / "pilot.yaml"


@dataclass(frozen=True)
class Config:
    raw: dict[str, Any] = field(repr=False)

    # --- convenience accessors (fail loudly if a key is missing) ---
    @property
    def datastack(self) -> str:
        return self.raw["data"]["datastack"]

    @property
    def materialization_version(self) -> int:
        return int(self.raw["data"]["materialization_version"])

    @property
    def tables(self) -> dict[str, str]:
        return self.raw["data"]["tables"]

    def table(self, key: str) -> str:
        try:
            return self.raw["data"]["tables"][key]
        except KeyError as e:  # noqa: PERF203
            raise KeyError(
                f"No table '{key}' in config/pilot.yaml. Known: {sorted(self.tables)}"
            ) from e

    @property
    def gates(self) -> dict[str, Any]:
        return self.raw["gates"]

    @property
    def validation(self) -> dict[str, Any]:
        return self.raw["validation"]

    @property
    def seed(self) -> int:
        return int(self.raw["validation"]["random_seed"])

    @property
    def compartments(self) -> list[str]:
        return list(self.raw["compartments"])

    def path(self, key: str) -> Path:
        """Resolve a configured path relative to the repo root and ensure it exists."""
        p = REPO_ROOT / self.raw["paths"][key]
        p.mkdir(parents=True, exist_ok=True)
        return p


@lru_cache(maxsize=1)
def load_config(path: str | Path = DEFAULT_CONFIG) -> Config:
    with open(path) as f:
        raw = yaml.safe_load(f)
    return Config(raw=raw)


CFG = load_config()
