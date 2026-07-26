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
# `tag` FIRST: that is the real column on synapse_target_predictions_ssa_v2 (verified 2026-07-26),
# carrying 'soma' / 'shaft' / 'spine'. It was absent from this list, so _col() returned None and
# every synapse fell through to "unknown" — silently zeroing every M4 compartment fraction.
_COMPARTMENT_COL = ["compartment", "tag", "pred_compartment", "target_compartment", "label"]
_PRECLASS_COL = ["pre_ei", "pre_cell_type", "pre_class"]
# Fine presynaptic identity, for SOURCE composition/diversity (M5). Distinct from the coarse E/I
# label above: entropy over pre_ei within inhibitory synapses is 0 for everyone (see data/join.py).
_PRESOURCE_COL = ["pre_mtype", "pre_cell_type", "pre_class"]


def _col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _map_compartment(raw: pd.Series) -> pd.Series:
    """Map raw spine/shaft/soma-style labels to our 5 compartments (config.compartments).

    ⚠️ HONEST LIMITATION — read before interpreting any M4 result.
    `synapse_target_predictions_ssa_v2` supplies exactly three labels: 'soma', 'shaft', 'spine'
    (verified 2026-07-26). That is a *structural-target* axis, NOT the radial soma→apical axis that
    docs/02 M4 and the docs/00 §5 mechanistic hypothesis are written in terms of. Recovering
    perisomatic / proximal / distal-basal / apical requires **distance-to-soma along the skeleton**
    plus a basal/apical branch call — i.e. pcg-skel / meshparty skeletonization, which is not yet in
    this pipeline.

    Until that exists, the mapping below is a stand-in: shaft→proximal and spine→distal_basal are
    naming conventions, not measured radial positions, and 'apical' will never be populated. An M4
    built on it tests "soma vs shaft vs spine", which is a real and defensible question — but it is
    not the pre-registered placement hypothesis, and must not be reported as though it were.
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
    """True for inhibitory presynaptic partners.

    aibs_metamodel_celltypes_v661.classification_system emits 'inhibitory_neuron' /
    'excitatory_neuron' / 'nonneuron', so a plain 'inhibitory' prefix suffices. The old bare 'i'
    prefix was a hazard waiting to fire: any future label starting with i — 'ITC', the
    interneuron-targeting class — would have been swept in as inhibitory.
    """
    s = pre_ei.astype(str).str.lower()
    return s.str.startswith("inhibitory") | s.isin({"inh", "gaba"})


def build_fingerprints(
    synapses: pd.DataFrame,
    *,
    dendrite_length: pd.Series | None = None,
    pilot: bool = True,
) -> pd.DataFrame:
    """One row of inhibitory-fingerprint features per post-synaptic root id.

    `synapses` is data/processed/incoming_synapses.parquet (pre_ei, pre_mtype and compartment
    already attached in data/join.py). `dendrite_length` (optional, indexed by root id) enables
    per-unit-length normalization; if absent, we fall back to total-input normalization only.

    `pilot=True` reserves the *motif* block (shared/convergent/disinhibitory structure) for full
    development (docs/03 Days 6-7). It does NOT gate the amount/location/source/diversity blocks —
    all of those are needed for the M0..M5 ladder to mean what docs/02 says it means.
    """
    df = synapses.copy()
    pre_col = _col(df, _PRECLASS_COL)
    src_col = _col(df, _PRESOURCE_COL)
    comp_col = _col(df, _COMPARTMENT_COL)
    if pre_col is None:
        raise ValueError(f"No presynaptic class column found in {list(df.columns)}")
    if src_col == pre_col:
        # Only the coarse E/I label is available: source diversity would be a constant 0. Emit the
        # columns as NaN rather than as a fake zero, so the ladder can't read "no source signal".
        src_col = None

    df["is_inh"] = _is_inhibitory(df[pre_col])
    # "Typed neuronal input". Note astype("string") not astype(str): the latter renders NaN as the
    # literal "none", which passes a not-nonneuron test and would quietly count every untypable
    # orphan fragment as a neuron — restoring the all-incoming denominator we just removed.
    _lab = df[pre_col].astype("string").str.lower()
    df["is_neuron"] = _lab.notna() & ~_lab.str.startswith("nonneuron").fillna(False)
    df["compartment"] = (
        _map_compartment(df[comp_col]) if comp_col else "unknown"
    )

    rows = []
    for root_id, g in df.groupby("post_pt_root_id"):
        # Denominator = TYPED NEURONAL inputs. Not all incoming: 91.4% of a neuron's synapses come
        # from untypable orphan fragments, so n_inh/all_incoming (median 0.049) largely measures
        # local reconstruction density. Over typed inputs the median is 0.588 — the biologically
        # meaningful quantity. Reconstruction completeness is carried separately as `frac_typed`
        # in the M0 technical block, so the ladder controls for it explicitly.
        # Non-neuronal partners (astrocyte/oligo/microglia) are excluded from the denominator
        # rather than counted as input.
        typed = g[g["is_neuron"]]
        inh = typed[typed["is_inh"]]
        n_typed = len(typed)
        n_inh = len(inh)
        feat: dict[str, float] = {
            "pt_root_id": int(root_id),
            # --- amount ---
            "inh_synapse_count": n_inh,
            "n_typed_input": n_typed,
            "inh_fraction": n_inh / n_typed if n_typed else np.nan,
            "n_inh_source_neurons": inh["pre_pt_root_id"].nunique() if n_inh else 0,
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

        # --- source composition + diversity (M5) ---
        # Always emitted, including in the pilot: M5 is DEFINED by these columns, and _columns_for
        # silently skips absent ones, so omitting them would collapse M5 onto M4 and make the
        # pre-registered "does source identity matter?" test read as a null result (docs/02 §2).
        src = (_fractions(inh[src_col], sorted(inh[src_col].dropna().unique()))
               if n_inh and src_col else {})
        feat["inh_source_entropy"] = _entropy(list(src.values())) if src else np.nan
        feat["inh_dominant_source_frac"] = max(src.values()) if src else np.nan
        feat["inh_effective_n_classes"] = _effective_n(list(src.values())) if src else np.nan
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
