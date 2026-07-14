"""Days 4-5 — build R_i (oracle reliability) and M_i (state modulation) + the phenotype gate.

Requires DANDI 000402 downloaded (docs/01 §4). Iterates scans, computes per-neuron phenotypes, and
enforces the split-half test-retest gate before a neuron's M_i is trusted.
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from inhibitome.config import CFG
from inhibitome.data.functional import list_nwb_files, load_scan
from inhibitome.phenotypes.reliability import oracle_score_multiclip
from inhibitome.phenotypes.state_modulation import fit_state_modulation, test_retest


def main() -> int:
    files = list_nwb_files()  # raises with a helpful dandi-download hint if missing
    retest_min = CFG.gates["retest_reliability_min"]

    rows = []
    for f in files:
        scan = load_scan(f)  # NotImplementedError until NWB schema resolved on Day 4 (see functional.py)
        # stim_pred: from a frozen stimulus-only model or the released digital-twin prediction.
        # For the pilot we accept a per-scan stimulus predictor passed alongside the NWB; see docs/03.
        stim_pred = _stimulus_prediction(scan)
        for u, unit_id in enumerate(scan.unit_ids):
            y = scan.activity[u]
            clip_resp = _repeat_matrix(y, scan)
            R = oracle_score_multiclip(clip_resp)
            M = fit_state_modulation(y, stim_pred[u], scan.locomotion, scan.pupil, scan.eye_xy,
                                     seed=CFG.seed)
            retest = test_retest(y, stim_pred[u], scan.locomotion, scan.pupil, scan.eye_xy,
                                 seed=CFG.seed)
            rows.append({
                "session": scan.session, "scan_idx": scan.scan_idx, "unit_id": int(unit_id),
                "reliability": R,
                "loco_given_pupil": M.loco_given_pupil, "pupil_given_loco": M.pupil_given_loco,
                "multiplicative_gain": M.multiplicative_gain, "additive_shift": M.additive_shift,
                "state_r2_oos": M.r2_behavior_oos, "state_uncertainty": M.uncertainty,
                "retest": retest, "retest_pass": bool(np.isfinite(retest) and retest >= retest_min),
            })

    df = pd.DataFrame(rows)
    out = CFG.path("processed") / "phenotypes.parquet"
    df.to_parquet(out)

    n_pass = int(df["retest_pass"].sum())
    (CFG.path("outputs") / "phenotype_qc.md").write_text(
        f"# Phenotype QC (Days 4-5)\n\n"
        f"- Units with phenotypes: {len(df)}\n"
        f"- Passing M_i split-half gate (>= {retest_min}): **{n_pass}**\n"
        f"- Median reliability R_i: {df['reliability'].median():.3f}\n"
    )
    print(f"Wrote {out}; {n_pass}/{len(df)} units pass the M_i test-retest gate.")
    return 0


def _stimulus_prediction(scan) -> np.ndarray:
    """Placeholder for the stimulus-only prediction f_i(S_t) per unit.

    Day 4-5: supply either the released digital-twin prediction aligned to this scan, or fit a frozen
    stimulus-only encoding model. Shape (n_units, T). Kept explicit so the residual model in
    state_modulation.py is genuinely stimulus-aware (docs/03).
    """
    raise NotImplementedError(
        "Provide f_i(S_t): released digital-twin prediction or a frozen stimulus-only model. "
        "See docs/03 Days 4-5."
    )


def _repeat_matrix(y: np.ndarray, scan) -> dict:
    """Slice a unit's activity into {clip_id: (n_repeats, n_frames)} using scan.oracle_repeats."""
    out = {}
    for cid, runs in scan.oracle_repeats.items():
        length = min(e - s for s, e in runs)
        out[cid] = np.stack([y[s : s + length] for s, e in runs])
    return out


if __name__ == "__main__":
    sys.exit(main())
