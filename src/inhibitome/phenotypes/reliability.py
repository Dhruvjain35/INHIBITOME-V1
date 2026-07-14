"""R_i — trial-to-trial visual response reliability (the "oracle" score).

Definition (docs/01 §5): for repeated presentations of the same movie clip, oracle reliability is the
mean correlation between each presentation's response and the leave-one-out average of the others.

We compute this OURSELVES from the repeated-movie responses (never trust a released score blind) and,
if a released score exists, use it only as a cross-check.

Reliability must not be confused with mean activity, selectivity, amplitude, or digital-twin fit.
"""
from __future__ import annotations

import numpy as np


def oracle_score(responses: np.ndarray) -> float:
    """Leave-one-out oracle correlation for one neuron.

    Parameters
    ----------
    responses : (n_repeats, n_frames)
        The neuron's response time-course on each repeated presentation of the same clip, resampled
        onto a common frame grid.

    Returns
    -------
    float : mean over repeats of Corr(repeat_k, mean of the other repeats). NaN if < 2 usable repeats.
    """
    r = np.asarray(responses, dtype=float)
    if r.ndim != 2 or r.shape[0] < 2:
        return np.nan
    k = r.shape[0]
    corrs = []
    for i in range(k):
        held = r[i]
        rest = np.delete(r, i, axis=0).mean(axis=0)
        if np.std(held) < 1e-9 or np.std(rest) < 1e-9:
            continue
        corrs.append(np.corrcoef(held, rest)[0, 1])
    return float(np.mean(corrs)) if corrs else np.nan


def oracle_score_multiclip(clip_responses: dict[int, np.ndarray]) -> float:
    """Average oracle score across all oracle clips for one neuron.

    clip_responses: {clip_id: (n_repeats, n_frames)}. Averaging across clips gives a single, more
    stable R_i per neuron.
    """
    scores = [oracle_score(v) for v in clip_responses.values()]
    scores = [s for s in scores if np.isfinite(s)]
    return float(np.mean(scores)) if scores else np.nan


def split_half_reliability(clip_responses: dict[int, np.ndarray], seed: int) -> tuple[float, float]:
    """Test-retest QC: oracle scores computed on two disjoint halves of the repeats.

    Returns (R_half_A, R_half_B). Downstream we require these to be consistent for the neuron to be
    retained (analogous to the M_i split-half gate). Useful as a per-neuron reliability-of-reliability.
    """
    rng = np.random.default_rng(seed)
    a, b = {}, {}
    for cid, resp in clip_responses.items():
        k = resp.shape[0]
        if k < 4:
            continue
        idx = rng.permutation(k)
        half = k // 2
        a[cid], b[cid] = resp[idx[:half]], resp[idx[half:]]
    return oracle_score_multiclip(a), oracle_score_multiclip(b)
