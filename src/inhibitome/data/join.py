"""Days 1-2 — build the master table: one row per unique EM neuron.

Joins the manual coregistration cohort to cell identity, area/depth, proofreading, and incoming
synapses (with pre-synaptic E/I labels + compartment predictions). Restricts to EXCITATORY
post-synaptic neurons (the study cohort). See docs/03_TEN_DAY_PILOT.md Days 1-2.

Synapse-level data is kept in a separate long table (data/processed/incoming_synapses.parquet) and
only aggregated into fingerprints in Aim 2 — we don't want 0.5B-row semantics leaking into the neuron
table.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from inhibitome.config import CFG
from inhibitome.data.cave import Cave


def build_master(cave: Cave | None = None, *, cohort: str = "coreg_manual") -> dict[str, Path]:
    """Assemble and persist the master neuron table + incoming-synapse long table.

    Returns the paths written. Idempotent: relies on the CAVE query cache.
    """
    cave = cave or Cave()
    proc = CFG.path("processed")

    # 1) Cohort: human-verified EM<->functional matches (the bridge to DANDI function).
    coreg = cave.query_table(cohort)
    # Expected columns include: pt_root_id, session, scan_idx, unit_id, residual, score.
    root_ids = coreg["pt_root_id"].dropna().astype("int64").unique().tolist()

    # 2) Cell identity — keep excitatory post-synaptic neurons only.
    ei = cave.query_table("ei_coarse", filter_in={"pt_root_id": root_ids})
    mtypes = cave.query_table("mtypes", filter_in={"pt_root_id": root_ids})
    exc_ids = _excitatory_root_ids(ei)

    # 3) Area assignment + soma position (depth via standard-transform downstream).
    area = cave.query_table("functional_area", filter_in={"pt_root_id": exc_ids})
    proof = cave.query_table("proofreading", filter_in={"pt_root_id": exc_ids})

    # 4) Incoming synapses onto the excitatory cohort, with pre-synaptic E/I + compartment.
    syn = cave.synapses_onto(exc_ids)
    syn = _annotate_synapses(cave, syn)

    # 5) One row per EM neuron (functional ROIs stay in coreg; a neuron may map to several).
    master = (
        coreg[coreg["pt_root_id"].isin(exc_ids)]
        .merge(_reduce(mtypes, "pt_root_id"), on="pt_root_id", how="left")
        .merge(_reduce(area, "pt_root_id"), on="pt_root_id", how="left")
        .merge(_reduce(proof, "pt_root_id"), on="pt_root_id", how="left")
    )

    master_path = proc / "master_neurons.parquet"
    syn_path = proc / "incoming_synapses.parquet"
    master.to_parquet(master_path)
    syn.to_parquet(syn_path)
    return {"master": master_path, "synapses": syn_path}


def _excitatory_root_ids(ei: pd.DataFrame) -> list[int]:
    """Root ids classified excitatory. Column literal ('classification_system'/'cell_type') varies
    by table version — resolve defensively and record what we used."""
    col = _first_present(ei, ["cell_type", "classification_system", "pred_cell_type", "class"])
    if col is None:
        raise ValueError(f"Could not find an E/I class column in {list(ei.columns)}")
    exc_mask = ei[col].astype(str).str.lower().str.startswith(("exc", "e", "pyr", "23p", "4p",
                                                              "5p", "6p"))
    return ei.loc[exc_mask, "pt_root_id"].astype("int64").unique().tolist()


def _annotate_synapses(cave: Cave, syn: pd.DataFrame) -> pd.DataFrame:
    """Attach pre-synaptic E/I label and post-synaptic compartment prediction to each synapse."""
    if syn.empty:
        return syn
    pre_ids = syn["pre_pt_root_id"].dropna().astype("int64").unique().tolist()
    pre_ei = cave.query_table("ei_coarse", filter_in={"pt_root_id": pre_ids})
    col = _first_present(pre_ei, ["cell_type", "classification_system", "pred_cell_type", "class"])
    pre_ei = pre_ei.rename(columns={"pt_root_id": "pre_pt_root_id", col: "pre_ei"})[
        ["pre_pt_root_id", "pre_ei"]
    ]
    syn = syn.merge(pre_ei, on="pre_pt_root_id", how="left")

    # Compartment prediction table is keyed by synapse id; merge on the shared id column.
    comp = cave.query_table("synapse_compartment")
    id_col = _first_present(comp, ["target_id", "id_ref", "synapse_id", "id"])
    syn_id = _first_present(syn, ["id", "synapse_id"])
    if id_col and syn_id:
        keep = [id_col] + [c for c in comp.columns if "compartment" in c.lower()
                          or "target" in c.lower()]
        syn = syn.merge(comp[keep].rename(columns={id_col: syn_id}), on=syn_id, how="left")
    return syn


def _reduce(df: pd.DataFrame, key: str) -> pd.DataFrame:
    """Collapse to one row per key (first non-null), so merges stay one-row-per-neuron."""
    if df.empty:
        return df
    return df.sort_values(key).groupby(key, as_index=False).first()


def _first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None
