"""The nested model hierarchy M0..M5 and the incremental dR2 test (docs/02 §2).

Each model = a set of feature-column groups. We fit a regularized linear model (ridge) by default —
interpretable, and a flexible learner is only a *control* (see models/gbt via `estimator=`). The
headline is the held-out increment between nested models, with scan-clustered bootstrap CIs.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from inhibitome.models.validation import Fold, leave_one_scan_out, r2_oos

# ---- Feature blocks (column names produced by the join + fingerprint steps). ----
# Adjust the concrete column lists to what data/join.py + fingerprints emit; the STRUCTURE is fixed.
FEATURE_BLOCKS: dict[str, list[str]] = {
    "technical": ["depth", "imaging_quality", "em_boundary_dist"],           # M0
    "cellular": ["area", "layer", "mtype", "dendrite_length"],               # M1 adds
    "functional": ["baseline_activity", "tuning", "response_amplitude", "reliability"],  # M2 adds
    "total_input": ["total_exc_input", "total_inh_input"],                   # M3 adds
    "inh_compartment": ["inh_frac_soma", "inh_frac_proximal", "inh_frac_distal_basal",
                        "inh_frac_apical"],                                  # M4 adds
    "inh_fingerprint": ["inh_source_entropy", "inh_dominant_source_frac",
                        "inh_effective_n_classes", "inh_multisynaptic_frac",
                        "n_inh_source_neurons"],                             # M5 adds
}

# Ordered cumulative membership: model -> blocks included.
LADDER: dict[str, list[str]] = {
    "M0": ["technical"],
    "M1": ["technical", "cellular"],
    "M2": ["technical", "cellular", "functional"],
    "M3": ["technical", "cellular", "functional", "total_input"],
    "M4": ["technical", "cellular", "functional", "total_input", "inh_compartment"],
    "M5": ["technical", "cellular", "functional", "total_input", "inh_compartment",
           "inh_fingerprint"],
}


@dataclass
class LadderResult:
    endpoint: str
    r2: dict[str, float] = field(default_factory=dict)          # model -> R2_oos
    per_fold: dict[str, dict] = field(default_factory=dict)
    increments: dict[str, dict] = field(default_factory=dict)  # "M5-M2" -> {dR2, ci_lo, ci_hi, p}


def _columns_for(model: str, df: pd.DataFrame, drop: frozenset[str] = frozenset()) -> list[str]:
    cols: list[str] = []
    for block in LADDER[model]:
        cols += [c for c in FEATURE_BLOCKS[block] if c in df.columns and c not in drop]
    return cols


def _make_fit_predict(cols: list[str], estimator=None, alpha: float = 1.0):
    def fit_predict(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
        Xtr = pd.get_dummies(train[cols], drop_first=True)
        Xte = pd.get_dummies(test[cols], drop_first=True).reindex(columns=Xtr.columns, fill_value=0)
        est = estimator or make_pipeline(StandardScaler(with_mean=False), Ridge(alpha=alpha))
        est.fit(Xtr.fillna(0), train["__y__"])
        return est.predict(Xte.fillna(0))
    return fit_predict


def run_ladder(
    df: pd.DataFrame,
    y_col: str,
    *,
    scan_key: list[str],
    group_key: str,
    endpoint: str,
    exclude_reliability: bool = False,
    estimator=None,
    n_boot: int = 2000,
    seed: int = 0,
) -> LadderResult:
    """Fit M0..M5 with leave-one-scan-out and compute key increments with scan-clustered bootstrap.

    `exclude_reliability=True` for the reliability endpoint (drop 'reliability' from the functional
    block, since it is the target).
    """
    df = df.copy()
    df["__y__"] = df[y_col]
    # For the reliability endpoint the target must not also be a predictor (docs/02 §2, "M2'").
    # Scoped per call — never mutate FEATURE_BLOCKS, or the next endpoint inherits the exclusion.
    drop = frozenset({"reliability"}) if exclude_reliability else frozenset()

    res = LadderResult(endpoint=endpoint)
    for model in LADDER:
        cols = _columns_for(model, df, drop)
        folds = list(leave_one_scan_out(df, scan_key, group_key))
        out = r2_oos(df, "__y__", _make_fit_predict(cols, estimator), iter(folds))
        res.r2[model] = out["r2_oos"]
        res.per_fold[model] = out["per_fold"]

    # Pre-registered increments (docs/02 §2).
    for hi, lo in [("M5", "M2"), ("M5", "M3"), ("M4", "M3"), ("M5", "M4")]:
        res.increments[f"{hi}-{lo}"] = _increment(
            df, scan_key, group_key, hi, lo, estimator, n_boot, seed, drop
        )
    return res


def increment_point(df, scan_key, group_key, hi, lo, estimator=None,
                    drop: frozenset[str] = frozenset()) -> float:
    """Point estimate of dR2 = R2_oos(hi) - R2_oos(lo) on the data as given (no resampling).

    This is the *observed statistic* — what nulls and bootstrap CIs are compared against. It must
    be available independently of `n_boot`: permutation tests run with the bootstrap switched off.
    """
    def _r2_for(model: str) -> float:
        return r2_oos(df, "__y__",
                      _make_fit_predict(_columns_for(model, df, drop), estimator),
                      leave_one_scan_out(df, scan_key, group_key))["r2_oos"]

    return float(_r2_for(hi) - _r2_for(lo))


def _increment(df, scan_key, group_key, hi, lo, estimator, n_boot, seed,
               drop: frozenset[str] = frozenset()) -> dict:
    """Observed dR2 plus a scan-clustered bootstrap CI (resamples whole scans).

    `dR2` is always the point estimate on the real data; the bootstrap only supplies the interval.
    With n_boot=0 the point estimate still comes back and the CI fields are NaN.
    """
    obs = increment_point(df, scan_key, group_key, hi, lo, estimator, drop)
    rng = np.random.default_rng(seed)
    scan_id = df[scan_key].astype(str).agg("|".join, axis=1)
    scans = np.array(sorted(scan_id.unique()))
    diffs = []
    for _ in range(n_boot):
        pick = rng.choice(scans, size=len(scans), replace=True)
        boot = pd.concat([df[scan_id == s] for s in pick], ignore_index=True)
        b_scan = boot[scan_key].astype(str).agg("|".join, axis=1)
        if b_scan.nunique() < 2:
            continue
        diffs.append(increment_point(boot, scan_key, group_key, hi, lo, estimator, drop))
    diffs = np.array([d for d in diffs if np.isfinite(d)])
    if len(diffs) == 0:
        return {"dR2": obs, "ci_lo": np.nan, "ci_hi": np.nan, "p_ge_0": np.nan}
    return {
        "dR2": obs,
        "ci_lo": float(np.percentile(diffs, 2.5)),
        "ci_hi": float(np.percentile(diffs, 97.5)),
        "p_ge_0": float((diffs <= 0).mean()),  # frac of bootstrap where increment <= 0
    }
