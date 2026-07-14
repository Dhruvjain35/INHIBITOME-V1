"""Build inhibitory fingerprints Z_i from the incoming-synapse long table.

Sections mirror docs/00 Aim 2:
  amount     — how much inhibition (counts, fractions, per-dendrite, multisynaptic)
  location   — where it lands (soma/proximal/distal_basal/apical compartment fractions)
  source     — who provides it (presynaptic inhibitory class composition)
  diversity  — how concentrated/varied (entropy, dominant fraction, effective #classes)

RULE: raw counts are never used un-normalized (bigger dendrites collect more synapses). The PILOT
fingerprint (docs/03 Days 6-7) uses only amount + a coarse perisomatic-vs-distal split + source-neuron
count + broad source-class composition; motifs and fine classes are full-development.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inhibitome.config import CFG


# --- column resolution (schema literals vary by table version) ---------------
_COMPARTMENT_COL = ["compartment", "pred_compartment", "target_compartment", "label"]
_PRECLASS_COL = ["pre_ei", "pre_cell_type", "pre_class"]


def _col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _map_compartment(raw: pd.Series) -> pd.Series:
    """Map raw spine/shaft/soma-style labels to our 5 compartments (config.compartments).

    synapse_target_predictions_ssa gives spine/shaft/soma; combined with distance-to-soma we refine
    into proximal/distal/apical. Here we do the coarse, robust mapping; the proximal/distal/apical
    split is validated against allen_v1_column_types_slanted_ref before being trusted (docs/03).
    """
    s = raw.astype(str).str.lower()
    out = pd.Series("unknown", index=raw.index)
    out[s.str.contains("soma")] = "soma"
    out[s.str.contains("shaft") | s.str.contains("proximal")] = "proximal"
    out[s.str.contains("distal") | s.str.contains("basal")] = "distal_basal"
    out[s.str.contains("apical") | s.str.contains("tuft")] = "apical"
    out[s.str.contains("spine")] = "distal_basal"  # inhibitory-on-spine is rare; treat as distal
    return out


def _is_inhibitory(pre_ei: pd.Series) -> pd.Series:
    s = pre_ei.astype(str).str.lower()
    return s.str.startswith(("inh", "i")) | s.isin({"inhibitory", "gaba"})


def build_fingerprints(
    synapses: pd.DataFrame,
    *,
    dendrite_length: pd.Series | None = None,
    pilot: bool = True,
) -> pd.DataFrame:
    """One row of inhibitory-fingerprint features per post-synaptic root id.

    `synapses` is data/processed/incoming_synapses.parquet (pre_ei + compartment already attached in
    data/join.py). `dendrite_length` (optional, indexed by root id) enables per-unit-length
    normalization; if absent, we fall back to total-input normalization only.
    """
    df = synapses.copy()
    pre_col = _col(df, _PRECLASS_COL)
    comp_col = _col(df, _COMPARTMENT_COL)
    if pre_col is None:
        raise ValueError(f"No presynaptic class column found in {list(df.columns)}")

    df["is_inh"] = _is_inhibitory(df[pre_col])
    df["compartment"] = (
        _map_compartment(df[comp_col]) if comp_col else "unknown"
    )

    rows = []
    for root_id, g in df.groupby("post_pt_root_id"):
        inh = g[g["is_inh"]]
        n_total = len(g)
        n_inh = len(inh)
        feat: dict[str, float] = {
            "pt_root_id": int(root_id),
            # --- amount ---
            "inh_synapse_count": n_inh,
            "inh_fraction": n_inh / n_total if n_total else np.nan,
            "n_inh_source_neurons": inh[
                inh.columns[inh.columns.str.contains("pre_pt_root_id")][0]
            ].nunique() if n_inh else 0,
        }
        if dendrite_length is not None and root_id in dendrite_length.index:
            L = float(dendrite_length.loc[root_id])
            feat["inh_per_um"] = n_inh / L if L > 0 else np.nan

        # multisynaptic connections (same presyn neuron -> >1 synapse)
        if n_inh:
            per_source = inh["pre_pt_root_id"].value_counts()
            feat["inh_multisynaptic_frac"] = float((per_source > 1).mean())
        else:
            feat["inh_multisynaptic_frac"] = np.nan

        # --- location (compartment fractions) ---
        comp_frac = _fractions(inh["compartment"], CFG.compartments)
        for c in CFG.compartments:
            feat[f"inh_frac_{c}"] = comp_frac.get(c, 0.0)
        feat["inh_frac_perisomatic"] = comp_frac.get("soma", 0.0) + comp_frac.get("proximal", 0.0)
        feat["inh_frac_dendritic"] = comp_frac.get("distal_basal", 0.0) + comp_frac.get("apical", 0.0)

        # --- source composition + diversity ---
        if n_inh and not pilot:
            src = _fractions(inh[pre_col], sorted(inh[pre_col].dropna().unique()))
            feat["inh_source_entropy"] = _entropy(list(src.values()))
            feat["inh_dominant_source_frac"] = max(src.values()) if src else np.nan
            feat["inh_effective_n_classes"] = _effective_n(list(src.values()))
        comp_vals = list(comp_frac.values())
        feat["inh_compartment_entropy"] = _entropy(comp_vals)
        rows.append(feat)

    return pd.DataFrame(rows)


def _fractions(labels: pd.Series, classes: list) -> dict:
    if len(labels) == 0:
        return {}
    counts = labels.value_counts(normalize=True)
    return {c: float(counts.get(c, 0.0)) for c in classes}


def _entropy(p: list[float]) -> float:
    p = np.array([x for x in p if x > 0])
    return float(-(p * np.log(p)).sum()) if len(p) else 0.0


def _effective_n(p: list[float]) -> float:
    """exp(entropy) = effective number of classes."""
    return float(np.exp(_entropy(p)))
