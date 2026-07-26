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

import numpy as np
import pandas as pd

from inhibitome.config import CFG
from inhibitome.data.cave import Cave


def build_master(cave: Cave | None = None, *, cohort: str = "coreg_manual") -> dict[str, Path]:
    """Assemble and persist the master neuron table + incoming-synapse long table.

    Returns the paths written. Idempotent: relies on the CAVE query cache.
    """
    cave = cave or Cave()
    proc = CFG.path("processed")
    typing = CFG.typing

    # 1) Cohort: human-verified EM<->functional matches (the bridge to DANDI function).
    coreg = cave.query_table(cohort)
    # Verified columns: pt_root_id, session, scan_idx, unit_id, field, residual, score.
    root_ids = coreg["pt_root_id"].dropna().astype("int64").unique().tolist()
    print(f"cohort: {len(coreg):,} ROIs / {len(root_ids):,} roots")

    # 2) Cell identity — keep excitatory post-synaptic neurons only.
    ei = cave.query_table(typing["primary"], filter_in={"pt_root_id": root_ids})
    mtypes = cave.query_table("mtypes", filter_in={"pt_root_id": root_ids})
    exc_ids = _excitatory_root_ids(ei)
    print(f"excitatory cohort neurons: {len(exc_ids):,}")

    # 3) Area assignment + soma position (depth via standard-transform downstream).
    area = cave.query_table("functional_area", filter_in={"pt_root_id": exc_ids})
    proof = cave.query_table("proofreading", filter_in={"pt_root_id": exc_ids})

    # 4) Incoming synapses. 91.4% of them come from orphan axon fragments with no soma, which no
    #    table can ever type — pulling those costs ~10x the rows and yields no source or
    #    compartment identity. We restrict the pull to typed presynaptic partners and recover the
    #    true total degree separately (step 5), so the untyped mass is still *counted*, just not
    #    materialized row by row.
    typed_pre = _typed_presynaptic_ids(cave) if typing["restrict_pull_to_typed"] else None
    if typed_pre is not None:
        print(f"typed presynaptic universe: {len(typed_pre):,} cells")
    syn, degree = cave.synapses_onto(exc_ids, keep_pre_ids=typed_pre)
    syn = _annotate_synapses(cave, syn)

    # 5) Reconstruction completeness per neuron. `frac_typed` goes into the M0 technical block so
    #    every fingerprint increment is measured *beyond* how well the neighbourhood reconstructed.
    typed_counts = syn.groupby("post_pt_root_id").size().rename("n_typed_input")
    completeness = (
        pd.concat([degree.rename("n_total_input"), typed_counts], axis=1)
        .fillna(0).astype("int64")
        .rename_axis("pt_root_id").reset_index()
    )
    completeness["frac_typed"] = (
        completeness["n_typed_input"] / completeness["n_total_input"].replace(0, pd.NA)
    ).astype(float)
    print(f"synapses kept: {len(syn):,} typed of {int(degree.sum()):,} total "
          f"(median frac_typed {completeness['frac_typed'].median():.1%})")

    # 6) One row per functional ROI (a neuron may map to several; coreg keeps them distinct).
    #    Every annotation table carries the same generic columns (id, created, valid, pt_position,
    #    pt_supervoxel_id), so merging them whole collides on the second merge:
    #    MergeError: 'suffixes' which cause duplicate columns {'id_x', ...}. Select what each table
    #    is actually for, and give it a name that says so.
    master = (
        coreg.loc[coreg["pt_root_id"].isin(exc_ids),
                  ["pt_root_id", "session", "scan_idx", "unit_id", "field",
                   "residual", "score", "pt_position"]]
        .rename(columns={"residual": "coreg_residual", "score": "coreg_score"})
        .merge(_pick(mtypes, {"cell_type": "mtype"}), on="pt_root_id", how="left")
        .merge(_pick(area, {"tag": "area"}), on="pt_root_id", how="left")
        .merge(_pick(proof, {"status_dendrite": "proof_dendrite",
                             "status_axon": "proof_axon",
                             "strategy_dendrite": "proof_strategy_dendrite"}),
               on="pt_root_id", how="left")
        .merge(completeness, on="pt_root_id", how="left")
    )
    master = _add_depth(master)
    print(f"master: {len(master):,} ROI rows / {master['pt_root_id'].nunique():,} neurons; "
          f"areas {master['area'].value_counts().to_dict()}")

    master_path = proc / "master_neurons.parquet"
    syn_path = proc / "incoming_synapses.parquet"
    master.to_parquet(master_path)
    syn.to_parquet(syn_path)
    return {"master": master_path, "synapses": syn_path}


def _typed_presynaptic_ids(cave: Cave) -> list[int]:
    """Every root id the primary cell-type table can label — the pull's presynaptic universe."""
    t = cave.query_table(CFG.typing["primary"])
    return t["pt_root_id"].dropna().astype("int64").unique().tolist()


def _excitatory_root_ids(ei: pd.DataFrame) -> list[int]:
    """Root ids classified excitatory by the primary cell-type table.

    Verified schema of aibs_metamodel_celltypes_v661: `classification_system` holds
    excitatory_neuron / inhibitory_neuron / nonneuron, and `cell_type` holds the fine label
    (23P, 4P, 5P-IT, BC, MC, NGC, BPC, astrocyte, oligo, microglia). Read the coarse column
    directly — the old prefix heuristic over `cell_type` would have swept in 'astrocyte' via its
    leading 'a'... and, worse, matched almost anything against the bare 'e' prefix.
    """
    col = CFG.typing["ei_column"]
    if col not in ei.columns:
        raise ValueError(f"Expected '{col}' in the cell-type table; got {list(ei.columns)}")
    exc = ei[col].astype(str).str.lower().str.startswith("excitatory")
    return ei.loc[exc, "pt_root_id"].astype("int64").unique().tolist()


def _annotate_synapses(cave: Cave, syn: pd.DataFrame) -> pd.DataFrame:
    """Attach pre-synaptic class labels and the post-synaptic compartment prediction per synapse.

    THREE presynaptic labels, kept separate on purpose (docs/00 Aim 2: "broad *and* fine labels
    kept separate"):
      pre_ei    — excitatory_neuron / inhibitory_neuron / nonneuron; decides what counts as
                  inhibitory, and lets non-neuronal partners (astrocyte, oligo, microglia) be
                  excluded rather than silently counted as input;
      pre_type  — fine cell type from the same table (BC, MC, NGC, BPC, 23P, 4P, ...);
      pre_mtype — m-type from the census table, whose inhibitory labels are *targeting* classes
                  (PTC perisomatic-targeting, DTC dendrite-targeting) — the vocabulary docs/00
                  Aim 2 is actually written in.
    The fine labels are not optional: source entropy over the coarse E/I label, within the
    inhibitory synapses, is identically zero for every neuron — M5 would be untestable.
    """
    if syn.empty:
        return syn
    pre_ids = syn["pre_pt_root_id"].dropna().astype("int64").unique().tolist()
    typing = CFG.typing

    pre_ei = cave.query_table(typing["primary"], filter_in={"pt_root_id": pre_ids})
    pre_ei = _reduce(pre_ei, "pt_root_id")[
        ["pt_root_id", typing["ei_column"], typing["type_column"]]
    ].rename(columns={"pt_root_id": "pre_pt_root_id",
                      typing["ei_column"]: "pre_ei",
                      typing["type_column"]: "pre_type"})
    syn = syn.merge(pre_ei, on="pre_pt_root_id", how="left")

    # Fine presynaptic m-type (the targeting classes M5 source diversity is computed over).
    pre_mt = cave.query_table(typing["source_label"], filter_in={"pt_root_id": pre_ids})
    mt_col = _first_present(pre_mt, ["cell_type", "pred_cell_type", "mtype", "class"])
    if mt_col is not None:
        pre_mt = pre_mt.rename(columns={"pt_root_id": "pre_pt_root_id", mt_col: "pre_mtype"})[
            ["pre_pt_root_id", "pre_mtype"]
        ]
        syn = syn.merge(_reduce(pre_mt, "pre_pt_root_id"), on="pre_pt_root_id", how="left")

    # Compartment predictions are keyed by synapse id. The table has 208.6M rows, so it is fetched
    # ONLY for the synapse ids we actually pulled — an unfiltered query does not complete.
    syn_id = _first_present(syn, ["id", "synapse_id"])
    if syn_id:
        comp = cave.compartment_for(syn[syn_id].dropna().astype("int64"))
        id_col = _first_present(comp, ["target_id", "id_ref", "synapse_id", "id"])
        if id_col:
            keep = [id_col] + [c for c in comp.columns
                               if "compartment" in c.lower() or "tag" in c.lower()]
            comp = comp[keep].rename(columns={id_col: syn_id})
            syn = syn.merge(_reduce(comp, syn_id), on=syn_id, how="left")
    return syn


def _pick(df: pd.DataFrame, cols: dict[str, str]) -> pd.DataFrame:
    """One row per pt_root_id, keeping only `cols` (mapped old->new). Avoids merge collisions."""
    present = {old: new for old, new in cols.items() if old in df.columns}
    missing = set(cols) - set(present)
    if missing:
        raise KeyError(f"expected {sorted(missing)} in table; have {list(df.columns)}")
    return (_reduce(df, "pt_root_id")[["pt_root_id", *present]]
            .rename(columns=present))


def _add_depth(master: pd.DataFrame) -> pd.DataFrame:
    """Cortical depth of each soma, in microns from the pia.

    `depth` is named in the M0 technical block (docs/02 §2) but nothing produced it, and
    `_columns_for` drops absent columns silently — so M0 was being fit almost empty, which inflates
    every increment measured against it. Depth also carries the imaging-depth artifact that null N6
    exists to rule out, so an absent depth column quietly removes that control too.

    minnie65 is a slanted volume; `standard-transform` supplies the published nm -> depth transform.
    """
    if "pt_position" not in master.columns:
        return master
    pos = np.stack(master["pt_position"].apply(
        lambda p: np.asarray(p, dtype=float) if p is not None else np.array([np.nan] * 3)
    ).to_numpy())
    try:
        from standard_transform import minnie_transform

        # pt_position is in voxels at the segmentation's 4x4x40 nm resolution.
        tform = minnie_transform(resolution=[4, 4, 40])
        master["depth"] = tform.apply(pos)[:, 1]   # column 1 is depth below pia, in microns
    except Exception as e:  # noqa: BLE001
        # Never silently continue without depth — M0 depends on it.
        raise RuntimeError(
            f"Could not compute cortical depth via standard-transform ({type(e).__name__}: {e}). "
            "M0 needs it; fix the transform rather than dropping the column."
        ) from e

    d = master["depth"]
    if not (d.between(0, 1200).mean() > 0.95):
        raise RuntimeError(
            f"Depths look wrong (median {d.median():.0f}um, range {d.min():.0f}-{d.max():.0f}). "
            "Mouse V1 somas sit roughly 50-700um below pia; check the input resolution."
        )
    return master


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
