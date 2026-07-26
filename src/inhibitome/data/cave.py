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
# Cohort neurons average ~3,000 incoming synapses (measured 2026-07-26 over 10 cohort cells:
# mean 3,037, median 2,294, max 5,971), so the cap binds at ~164 post_ids. 100 leaves headroom
# for hub neurons without making the pull needlessly chatty.
_SYNAPSE_BATCH = 100


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

    def synapses_onto(self, root_ids: Iterable[int], *, use_cache: bool = True,
                      batch_size: int = _SYNAPSE_BATCH, progress: bool = True) -> pd.DataFrame:
        """Incoming synapses for post-synaptic `root_ids`, chunked under the 500k cap.

        Returns synapses_pni_2 rows: pre_pt_root_id, post_pt_root_id, positions, size.

        Cohort neurons carry ~3,000 incoming synapses each (measured: median 2,294, max ~6,000), so
        the cap binds at ~160 post_ids. We batch well under that and, on the rare batch that still
        caps out, bisect rather than dropping to one-id-at-a-time — a 15k-neuron cohort would
        otherwise become 15k round trips.

        Each batch is cached separately, so a multi-hour pull resumes where it stopped instead of
        restarting from zero.
        """
        root_ids = list(dict.fromkeys(int(r) for r in root_ids))  # de-dup, keep order
        frames: list[pd.DataFrame] = []
        batches = list(_batched(root_ids, batch_size))
        for i, batch in enumerate(batches):
            cache = _cache_key("syn_onto_batch", ids=batch)
            if use_cache and cache.exists():
                frames.append(pd.read_parquet(cache))
                continue
            df = self._synapse_query_bisect(batch)
            df.to_parquet(cache)
            frames.append(df)
            if progress:
                done = sum(len(f) for f in frames)
                print(f"  synapses: batch {i + 1}/{len(batches)}  (+{len(df):,} -> {done:,} total)",
                      flush=True)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def _synapse_query_bisect(self, ids: list[int]) -> pd.DataFrame:
        """Query these post_ids; a result at the row cap is truncated, so split and retry."""
        df = self.client.materialize.synapse_query(
            post_ids=ids, materialization_version=self.version
        )
        if len(df) < _SYNAPSE_QUERY_CAP or len(ids) == 1:
            return df
        mid = len(ids) // 2
        return pd.concat(
            [self._synapse_query_bisect(ids[:mid]), self._synapse_query_bisect(ids[mid:])],
            ignore_index=True,
        )

    def compartment_for(self, synapse_ids: Iterable[int], *,
                        use_cache: bool = True, chunk: int = 200_000) -> pd.DataFrame:
        """Compartment predictions for specific synapse ids.

        `synapse_target_predictions_ssa_v2` has 208,644,969 rows. Querying it unfiltered does not
        complete — it must always be restricted to the synapse ids we actually pulled.
        """
        ids = [int(i) for i in pd.unique(pd.Series(list(synapse_ids), dtype="int64"))]
        name = CFG.table("synapse_compartment")
        frames: list[pd.DataFrame] = []
        chunks = list(_batched(ids, chunk))
        for i, part in enumerate(chunks):
            cache = _cache_key("syn_comp", name=name, n=len(part), head=part[:3], tail=part[-3:])
            if use_cache and cache.exists():
                frames.append(pd.read_parquet(cache))
                continue
            df = self.client.materialize.query_table(
                name, filter_in_dict={"target_id": part},
                materialization_version=self.version,
            )
            df.to_parquet(cache)
            frames.append(df)
            print(f"  compartments: chunk {i + 1}/{len(chunks)} (+{len(df):,})", flush=True)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _batched(seq: list, n: int):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]
