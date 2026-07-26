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
import time
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
from requests.adapters import HTTPAdapter

from inhibitome.config import CFG, REPO_ROOT

CACHE_DIR = REPO_ROOT / "data" / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

_SYNAPSE_QUERY_CAP = 500_000  # CAVE hard cap per query
# Cohort neurons average ~3,000 incoming synapses (measured 2026-07-26 over 10 cohort cells:
# mean 3,037, median 2,294, max 5,971), so the cap binds at ~164 post_ids. Measured timings on
# cohort ids: 5 -> 5.9s, 20 -> 16.1s, 50 -> 29.4s (~0.6s/neuron). 50 keeps each request short
# enough to survive a gateway hiccup; a 100-id request stalled indefinitely in a real run.
_SYNAPSE_BATCH = 50
# No timeout means one hung connection blocks a six-hour job forever — which is exactly what
# happened: batch 2 sat for 59 minutes on 2.6s of CPU while the same ids answered fine in
# smaller chunks moments later.
_REQUEST_TIMEOUT = (30, 300)  # (connect, read) seconds
_MAX_ATTEMPTS = 3


class _TimeoutAdapter(HTTPAdapter):
    """requests has no default-timeout setting; inject one on every call through this adapter."""

    def send(self, request, **kwargs):
        if kwargs.get("timeout") is None:
            kwargs["timeout"] = _REQUEST_TIMEOUT
        return super().send(request, **kwargs)


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
            _mount_timeouts(c)
            self._client = c
        return self._client

    def available_versions(self) -> list[int]:
        return sorted(self.client.materialize.get_versions())

    def query_table(self, table_key_or_name: str, *, filter_in: dict | None = None,
                    use_cache: bool = True, **kwargs) -> pd.DataFrame:
        """Query an annotation table by config key (preferred) or raw name.

        `filter_in={'pt_root_id': [...]}` is applied CLIENT-SIDE by default. The cell-type and
        coregistration tables are *reference* tables: `pt_root_id` belongs to the referenced
        nucleus table, not the annotation model, so a server-side filter on it 500s with
        `KeyError: 'pt_root_id'`. Every table we read this way is under 150k rows, so fetching the
        whole thing once and filtering in pandas is both correct and better cached — one copy per
        table, reused by every caller, instead of one cache entry per id list.

        Columns that DO live on the annotation model (e.g. `target_id`) can still be pushed to the
        server with `server_filter_in=`.
        """
        name = CFG.tables.get(table_key_or_name, table_key_or_name)
        server_filter = kwargs.pop("server_filter_in", None)
        cache = _cache_key("table", name=name, server_filter=server_filter, kwargs=kwargs)

        if use_cache and cache.exists():
            df = pd.read_parquet(cache)
        else:
            df = self.client.materialize.query_table(
                name,
                filter_in_dict=server_filter,
                materialization_version=self.version,
                **kwargs,
            )
            df.to_parquet(cache)

        if filter_in:
            for col, values in filter_in.items():
                if col not in df.columns:
                    raise KeyError(f"'{col}' not in {name}; have {list(df.columns)}")
                df = df[df[col].isin(set(values))]
            df = df.reset_index(drop=True)
        return df

    def synapses_onto(self, root_ids: Iterable[int], *, keep_pre_ids: Iterable[int] | None = None,
                      use_cache: bool = True, batch_size: int = _SYNAPSE_BATCH,
                      progress: bool = True) -> tuple[pd.DataFrame, pd.Series]:
        """Incoming synapses for `root_ids`, plus each neuron's TOTAL in-degree.

        Returns `(synapses, degree)`:
          synapses — synapses_pni_2 rows, filtered to `keep_pre_ids` if given;
          degree   — total incoming count per root id, counting *every* partner.

        Why one pass and not two: 91.4% of incoming synapses come from orphan axon fragments with
        no soma, which no table can type and which therefore carry no source or compartment
        identity. We do not want them in the stored table (~47M rows vs ~4M), but we do need them
        in the *denominator* — otherwise `inh_fraction` silently encodes local reconstruction
        density instead of inhibition. So each batch is fetched once, counted in full, and only
        the typed rows are kept. The transfer is unavoidable; the storage and every downstream
        groupby shrink ~10x.

        There is no cheaper exact route: `synapse_query` has no `count` parameter, and the
        `synapses_pni_2_in_out_degree` view rejects `pt_root_id` as a filter (500) while an
        unfiltered read scans 337M rows.

        Cohort neurons carry ~3,000 incoming synapses each (measured: median 2,294, max ~6,000), so
        the 500k cap binds at ~160 post_ids. We batch under that and bisect a batch that still caps
        out, rather than dropping to one-id-at-a-time — a 15k cohort would become 15k round trips.

        Each batch is cached separately, so a long pull resumes where it stopped.
        """
        root_ids = list(dict.fromkeys(int(r) for r in root_ids))
        keep = set(int(p) for p in keep_pre_ids) if keep_pre_ids is not None else None

        frames: list[pd.DataFrame] = []
        degrees: list[pd.Series] = []
        batches = list(_batched(root_ids, batch_size))
        for i, batch in enumerate(batches):
            syn_cache = _cache_key("syn_typed_batch", ids=batch, filtered=keep is not None)
            deg_cache = _cache_key("syn_degree_batch", ids=batch)
            if use_cache and syn_cache.exists() and deg_cache.exists():
                frames.append(pd.read_parquet(syn_cache))
                d = pd.read_parquet(deg_cache)
                degrees.append(d.set_index(d.columns[0])[d.columns[1]])
                continue

            raw = self._synapse_query_bisect(batch)
            deg = raw.groupby("post_pt_root_id").size()
            kept = raw[raw["pre_pt_root_id"].isin(keep)] if keep is not None else raw

            kept.to_parquet(syn_cache)
            deg.rename("n_total_input").reset_index().to_parquet(deg_cache)
            frames.append(kept)
            degrees.append(deg)
            if progress:
                n_all, n_keep = len(raw), len(kept)
                print(f"  synapses: batch {i + 1}/{len(batches)}  "
                      f"{n_all:,} seen -> {n_keep:,} typed ({n_keep / max(n_all, 1):.1%})",
                      flush=True)

        syn = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        degree = (pd.concat(degrees).groupby(level=0).sum() if degrees
                  else pd.Series(dtype="int64"))
        return syn, degree

    def _synapse_query_bisect(self, ids: list[int]) -> pd.DataFrame:
        """Query these post_ids, splitting on either the row cap or a failure.

        Two reasons to bisect. A result AT the cap is silently truncated, so it must be split to be
        correct. And a request that times out or errors is retried a few times, then split — a
        smaller request usually succeeds where a larger one hung, and halving beats collapsing to
        one-id-at-a-time (a 15k cohort would become 15k round trips).
        """
        last: Exception | None = None
        for attempt in range(_MAX_ATTEMPTS):
            try:
                df = self.client.materialize.synapse_query(
                    post_ids=ids, materialization_version=self.version
                )
                if len(df) < _SYNAPSE_QUERY_CAP or len(ids) == 1:
                    return df
                break  # hit the cap: fall through and split
            except Exception as e:  # noqa: BLE001 — transport/server errors are all retryable here
                last = e
                if attempt < _MAX_ATTEMPTS - 1:
                    wait = 5 * 2 ** attempt
                    print(f"    retry {attempt + 1}/{_MAX_ATTEMPTS - 1} on {len(ids)} ids after "
                          f"{type(e).__name__}; waiting {wait}s", flush=True)
                    time.sleep(wait)

        if len(ids) == 1:
            raise RuntimeError(f"synapse_query failed for post_id {ids[0]}") from last
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
            # target_id IS on the annotation model, so this one filters server-side.
            df = self.client.materialize.query_table(
                name, filter_in_dict={"target_id": part},
                materialization_version=self.version,
            )
            df.to_parquet(cache)
            frames.append(df)
            print(f"  compartments: chunk {i + 1}/{len(chunks)} (+{len(df):,})", flush=True)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _mount_timeouts(client) -> None:
    """Give every CAVE sub-client's session a request timeout.

    CAVEclient builds several sub-clients, each with its own requests.Session and no default
    timeout, so any of them can hang forever on a dropped connection.
    """
    seen = set()
    for attr in ("materialize", "chunkedgraph", "annotation", "auth", "info", "l2cache", "state"):
        sub = getattr(client, attr, None)
        session = getattr(sub, "session", None)
        if session is None or id(session) in seen:
            continue
        seen.add(id(session))
        adapter = _TimeoutAdapter(max_retries=0)  # retries handled by the caller, with bisection
        session.mount("https://", adapter)
        session.mount("http://", adapter)


def _batched(seq: list, n: int):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]
