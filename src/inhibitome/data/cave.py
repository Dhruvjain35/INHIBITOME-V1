"""CAVEclient wrapper: pinned materialization + on-disk query cache.

Grounded in the verified MICrONS API (docs/01_DATA_ACCESS.md):
  - datastack 'minnie65_public', a token IS required (even for public data);
  - materialization pinned via `client.version = <int>`;
  - synapses via `client.materialize.synapse_query(post_ids=...)`, 500k-row cap per query;
  - annotation tables via `client.materialize.query_table(name)`.

Every query is cached to data/cache/ keyed by (table, version, filter-hash) so re-runs are cheap and
CAVE isn't hammered. Delete data/cache/ to force a refresh.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from inhibitome.config import CFG, REPO_ROOT

CACHE_DIR = REPO_ROOT / "data" / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

_SYNAPSE_QUERY_CAP = 500_000  # CAVE hard cap per query


def _cache_key(kind: str, **parts: Any) -> Path:
    blob = json.dumps({"kind": kind, "v": CFG.materialization_version, **parts}, sort_keys=True,
                      default=str)
    h = hashlib.sha1(blob.encode()).hexdigest()[:16]
    safe = kind.replace("/", "_")
    return CACHE_DIR / f"{safe}__{h}.parquet"


class Cave:
    """Thin, cached CAVEclient facade pinned to the config materialization."""

    def __init__(self, datastack: str | None = None, version: int | None = None):
        self.datastack = datastack or CFG.datastack
        self.version = version or CFG.materialization_version
        self._client = None  # lazy: importing caveclient shouldn't be required to import this module

    @property
    def client(self):
        if self._client is None:
            from caveclient import CAVEclient  # local import keeps the package importable offline

            c = CAVEclient(self.datastack)
            c.version = self.version  # pin materialization for reproducibility
            self._client = c
        return self._client

    def available_versions(self) -> list[int]:
        return sorted(self.client.materialize.get_versions())

    def query_table(self, table_key_or_name: str, *, filter_in: dict | None = None,
                    use_cache: bool = True, **kwargs) -> pd.DataFrame:
        """Query an annotation table by config key (preferred) or raw name.

        `filter_in={'pt_root_id': [...]}` maps to CAVE's `filter_in_dict`.
        """
        name = CFG.tables.get(table_key_or_name, table_key_or_name)
        cache = _cache_key("table", name=name, filter_in=filter_in, kwargs=kwargs)
        if use_cache and cache.exists():
            return pd.read_parquet(cache)

        df = self.client.materialize.query_table(
            name,
            filter_in_dict=filter_in,
            materialization_version=self.version,
            **kwargs,
        )
        df.to_parquet(cache)
        return df

    def synapses_onto(self, root_ids: Iterable[int], *, use_cache: bool = True) -> pd.DataFrame:
        """Incoming synapses for post-synaptic `root_ids` (chunked under the 500k cap).

        Returns synapses_pni_2 rows: pre_pt_root_id, post_pt_root_id, positions, size.
        """
        root_ids = list(dict.fromkeys(int(r) for r in root_ids))  # de-dup, keep order
        cache = _cache_key("syn_onto", n=len(root_ids),
                          head=root_ids[:5], tail=root_ids[-5:])
        if use_cache and cache.exists():
            return pd.read_parquet(cache)

        frames: list[pd.DataFrame] = []
        # Query per batch of post ids; if any single neuron approaches the cap, split it out.
        for batch in _batched(root_ids, 200):
            df = self.client.materialize.synapse_query(
                post_ids=batch, materialization_version=self.version
            )
            if len(df) >= _SYNAPSE_QUERY_CAP:
                # Fell into the cap: re-query these post ids one at a time to be safe.
                df = pd.concat(
                    [self.client.materialize.synapse_query(
                        post_ids=[r], materialization_version=self.version) for r in batch],
                    ignore_index=True,
                )
            frames.append(df)
        out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        out.to_parquet(cache)
        return out


def _batched(seq: list, n: int):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]
