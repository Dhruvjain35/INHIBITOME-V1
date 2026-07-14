"""Grouped out-of-sample validation (docs/02 §4).

Two invariants, enforced here so no caller can violate them:
  1. GROUPING — all rows sharing an EM root_id stay in the same fold (one neuron may map to many ROIs).
  2. PRIMARY SPLIT — leave-one-scan-out: hold out an entire imaging scan at a time.

We report held-out R2 (R2_oos) pooled across folds, plus per-fold values (their spread is part of the
evidence). scan-clustered bootstrap CIs live in models/nested.py.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterator

import numpy as np
import pandas as pd


@dataclass
class Fold:
    name: str
    train_idx: np.ndarray
    test_idx: np.ndarray


def leave_one_scan_out(df: pd.DataFrame, scan_key: list[str], group_key: str) -> Iterator[Fold]:
    """Yield one fold per unique scan. Rows of a group_key are never split across train/test.

    Because coregistration is one-neuron-to-one-scan-region, grouping by root_id and holding out by
    scan are compatible; we still assert no group straddles the split.
    """
    scan_id = df[scan_key].astype(str).agg("|".join, axis=1)
    for scan in sorted(scan_id.unique()):
        test_mask = (scan_id == scan).to_numpy()
        train_groups = set(df.loc[~test_mask, group_key])
        test_groups = set(df.loc[test_mask, group_key])
        leaked = train_groups & test_groups
        if leaked:
            # A neuron appears in both -> force it fully into test to preserve the invariant.
            move = df[group_key].isin(leaked).to_numpy()
            test_mask = test_mask | move
        yield Fold(
            name=scan,
            train_idx=np.where(~test_mask)[0],
            test_idx=np.where(test_mask)[0],
        )


def grouped_random(df: pd.DataFrame, group_key: str, n_splits: int, seed: int) -> Iterator[Fold]:
    """SECONDARY diagnostic only: group-preserving random K-fold. Never the headline evidence."""
    rng = np.random.default_rng(seed)
    groups = df[group_key].to_numpy()
    uniq = np.array(sorted(pd.unique(groups)))
    rng.shuffle(uniq)
    chunks = np.array_split(uniq, n_splits)
    for k, held in enumerate(chunks):
        test_mask = np.isin(groups, held)
        yield Fold(f"randfold{k}", np.where(~test_mask)[0], np.where(test_mask)[0])


def r2_oos(
    df: pd.DataFrame,
    y_col: str,
    fit_predict: Callable[[pd.DataFrame, pd.DataFrame], np.ndarray],
    folds: Iterator[Fold],
) -> dict:
    """Pool held-out predictions across folds; return overall + per-fold R2.

    `fit_predict(train_df, test_df) -> yhat_test` isolates the model from the CV machinery so the
    same validation is reused for every M0..M5 and every null.
    """
    per_fold = {}
    y_true_all, y_pred_all = [], []
    for fold in folds:
        tr, te = df.iloc[fold.train_idx], df.iloc[fold.test_idx]
        yhat = np.asarray(fit_predict(tr, te), float)
        yt = te[y_col].to_numpy(float)
        m = np.isfinite(yhat) & np.isfinite(yt)
        per_fold[fold.name] = _r2(yt[m], yhat[m])
        y_true_all.append(yt[m])
        y_pred_all.append(yhat[m])
    yt = np.concatenate(y_true_all)
    yp = np.concatenate(y_pred_all)
    return {"r2_oos": _r2(yt, yp), "per_fold": per_fold,
            "median_fold": float(np.nanmedian(list(per_fold.values())))}


def _r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return float(1 - np.sum((y_true - y_pred) ** 2) / ss_tot) if ss_tot > 1e-12 else np.nan
