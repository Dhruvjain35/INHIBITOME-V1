"""Matched-null controls. The fingerprint is credible only if it beats these out-of-sample.

Each null destroys ONE aspect of the fingerprint while preserving the confounds, then re-runs the
same M5-vs-baseline increment. We compare the observed increment to the null distribution and report
an empirical p-value. See docs/02 §3 for the full table (N1-N7); implemented here: N2, N3, N4 (the
permutation family) plus a generic runner. N1/N5/N6/N7 are wired via feature swaps / regressors and
documented inline.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from inhibitome.models.nested import FEATURE_BLOCKS, run_ladder

# Columns that define a matching stratum for N2 (layer/morphology-matched permutation).
_MATCH_KEYS = ["area", "layer", "mtype"]
# The fingerprint columns that get shuffled.
_FINGERPRINT_COLS = FEATURE_BLOCKS["inh_compartment"] + FEATURE_BLOCKS["inh_fingerprint"]


def matched_permutation_null(df: pd.DataFrame, match_keys=None, cols=None,
                            seed: int = 0) -> pd.DataFrame:
    """N2 — shuffle fingerprints only AMONG neurons in the same layer/morphology/area stratum.

    Prevents the model from 'winning' just because cell types carry both wiring and function.
    """
    match_keys = match_keys or _MATCH_KEYS
    cols = cols or [c for c in _FINGERPRINT_COLS if c in df.columns]
    rng = np.random.default_rng(seed)
    out = df.copy()
    keys = [k for k in match_keys if k in df.columns]
    for _, idx in df.groupby(keys).groups.items() if keys else [("all", df.index)]:
        idx = np.asarray(list(idx))
        perm = rng.permutation(idx)
        out.loc[idx, cols] = df.loc[perm, cols].to_numpy()
    return out


def compartment_location_null(df: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    """N3 — preserve total inhibitory count; scramble which compartment the synapses hit.

    Randomly redistributes the compartment fractions across the config compartments per neuron while
    keeping their sum (and the neuron's total count) fixed.
    """
    rng = np.random.default_rng(seed)
    comp_cols = [c for c in FEATURE_BLOCKS["inh_compartment"] if c in df.columns]
    out = df.copy()
    fr = df[comp_cols].to_numpy(float)
    for i in range(len(fr)):
        row = fr[i]
        total = np.nansum(row)
        if total > 0:
            out.loc[out.index[i], comp_cols] = rng.dirichlet(np.ones(len(comp_cols))) * total
    return out


def presynaptic_identity_null(df: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    """N4 — preserve counts + compartments; shuffle source-class composition among matched inputs."""
    src_cols = [c for c in FEATURE_BLOCKS["inh_fingerprint"] if "source" in c or "classes" in c]
    return matched_permutation_null(df, cols=src_cols, seed=seed)


NULLS = {
    "N2_layer_morph": matched_permutation_null,
    "N3_compartment": compartment_location_null,
    "N4_presyn_id": presynaptic_identity_null,
}


def run_null(
    df: pd.DataFrame, y_col: str, null_name: str, *,
    scan_key, group_key, endpoint, hi="M5", lo="M2",
    n_perm: int = 1000, seed: int = 0, **ladder_kw,
) -> dict:
    """Empirical p-value: fraction of `n_perm` matched shuffles whose dR2 >= the observed dR2.

    A small p means the real fingerprint's increment is NOT reproducible by matched chance wiring
    (good). A large p means the null does just as well (kill signal — docs/04).
    """
    shuffle = NULLS[null_name]
    obs = run_ladder(df, y_col, scan_key=scan_key, group_key=group_key, endpoint=endpoint,
                     n_boot=0, **ladder_kw).increments[f"{hi}-{lo}"]["dR2"]
    null_dr2 = []
    for k in range(n_perm):
        sh = shuffle(df, seed=seed + k)
        r = run_ladder(sh, y_col, scan_key=scan_key, group_key=group_key, endpoint=endpoint,
                       n_boot=0, **ladder_kw).increments[f"{hi}-{lo}"]["dR2"]
        if np.isfinite(r):
            null_dr2.append(r)
    null_dr2 = np.array(null_dr2)
    p = float((null_dr2 >= obs).mean()) if len(null_dr2) else np.nan
    return {
        "null": null_name, "observed_dR2": float(obs),
        "null_mean": float(null_dr2.mean()) if len(null_dr2) else np.nan,
        "null_p95": float(np.percentile(null_dr2, 95)) if len(null_dr2) else np.nan,
        "p_value": p, "beats_null": bool(np.isfinite(p) and p < 0.05),
    }
