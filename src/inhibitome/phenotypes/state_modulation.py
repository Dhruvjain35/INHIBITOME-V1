"""M_i — locomotion/pupil-linked state modulation of visual responses.

A raw correlation between activity and running is NOT sufficient: running, pupil, stimulus, and eye
position are correlated. We fit a per-neuron STIMULUS-AWARE encoding model (docs/03 Days 4-5):

    y_i(t) = f_i(S_t) + b_L*L_t + b_P*P_t + b_E*E_t + b_LP*(L_t*P_t) + eps

`f_i(S_t)` is the stimulus prediction — the residual after removing it is what behavior explains.
We take the stimulus prediction from a frozen stimulus-only model (or the released digital-twin
prediction) and model the residual with behavior. Evaluation is on CONTIGUOUS held-out time blocks,
because behavior autocorrelates and random-sample CV would leak.

Every M_i must pass a split-half test-retest gate (config: retest_reliability_min) before use.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import Ridge


@dataclass
class StateModulation:
    """Per-neuron state-modulation summary (these become predictor targets in Aim 3)."""
    loco_given_pupil: float      # locomotion effect controlling for pupil
    pupil_given_loco: float      # pupil effect controlling for locomotion
    additive_shift: float        # baseline offset attributable to state
    multiplicative_gain: float   # response-gain change (from the L*stim interaction)
    r2_behavior_oos: float       # held-out R2 of the behavior model on stimulus residual
    uncertainty: float           # bootstrap SD of the primary coefficient


def _blocks(n: int, n_blocks: int) -> list[np.ndarray]:
    """Contiguous time blocks for block-CV (never shuffle time)."""
    edges = np.linspace(0, n, n_blocks + 1).astype(int)
    return [np.arange(edges[i], edges[i + 1]) for i in range(n_blocks)]


def fit_state_modulation(
    y: np.ndarray,            # (T,) neuron activity
    stim_pred: np.ndarray,    # (T,) stimulus-only prediction f_i(S_t)
    loco: np.ndarray,         # (T,)
    pupil: np.ndarray,        # (T,)
    eye_xy: np.ndarray,       # (T, 2)
    *,
    n_blocks: int = 5,
    alpha: float = 1.0,
    seed: int = 0,
) -> StateModulation:
    """Fit the residual behavior model with contiguous block-CV and return the summary."""
    y = np.asarray(y, float)
    resid = y - np.asarray(stim_pred, float)      # remove the stimulus drive
    L, P = _z(loco), _z(pupil)
    E = np.column_stack([_z(eye_xy[:, 0]), _z(eye_xy[:, 1])])
    X = np.column_stack([L, P, E, L * P, L * _z(stim_pred)])  # last col ~ multiplicative gain
    #             0  1  2:4  4(L*P)  5(L*stim)
    T = len(y)

    # Held-out R2 over contiguous blocks.
    blocks = _blocks(T, n_blocks)
    r2s, coefs = [], []
    for i, test in enumerate(blocks):
        train = np.concatenate([b for j, b in enumerate(blocks) if j != i])
        m = Ridge(alpha=alpha).fit(X[train], resid[train])
        pred = m.predict(X[test])
        r2s.append(_r2(resid[test], pred))
        coefs.append(m.coef_)
    coefs = np.array(coefs)
    coef_mean = coefs.mean(axis=0)

    return StateModulation(
        loco_given_pupil=float(coef_mean[0]),
        pupil_given_loco=float(coef_mean[1]),
        additive_shift=float(np.mean(resid)),
        multiplicative_gain=float(coef_mean[5]),
        r2_behavior_oos=float(np.mean(r2s)),
        uncertainty=float(coefs[:, 0].std()),
    )


def test_retest(
    y, stim_pred, loco, pupil, eye_xy, *, seed: int = 0, **kw
) -> float:
    """Split the recording into two contiguous halves, fit M_i on each, return their agreement.

    We correlate the primary coefficient vector across halves (loco, pupil, gain). The phenotype gate
    (config: retest_reliability_min) is applied to this value.
    """
    T = len(y)
    h = T // 2
    a = fit_state_modulation(y[:h], stim_pred[:h], loco[:h], pupil[:h], eye_xy[:h], seed=seed, **kw)
    b = fit_state_modulation(y[h:], stim_pred[h:], loco[h:], pupil[h:], eye_xy[h:], seed=seed, **kw)
    va = np.array([a.loco_given_pupil, a.pupil_given_loco, a.multiplicative_gain])
    vb = np.array([b.loco_given_pupil, b.pupil_given_loco, b.multiplicative_gain])
    if np.std(va) < 1e-9 or np.std(vb) < 1e-9:
        return np.nan
    return float(np.corrcoef(va, vb)[0, 1])


def _z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    s = x.std()
    return (x - x.mean()) / s if s > 1e-9 else x - x.mean()


def _r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return float(1 - ss_res / ss_tot) if ss_tot > 1e-12 else np.nan
