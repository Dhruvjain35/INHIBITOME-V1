"""Days 8-9 — nested-model ladder + matched nulls, then the Day-10 decision. See docs/02, docs/03.

Freezes the prereg (writes outputs/prereg_lock.json with the git SHA), runs M0..M5 with
leave-one-scan-out, benchmarks against N2/N3 matched nulls, and writes the decision report.
"""
from __future__ import annotations

import json
import subprocess
import sys

import pandas as pd

from inhibitome.config import CFG
from inhibitome.models.nested import run_ladder
from inhibitome.nulls import run_null


def _git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def _assemble() -> pd.DataFrame:
    """Join master + phenotypes + fingerprints into the modeling frame (one row per ROI)."""
    proc = CFG.path("processed")
    master = pd.read_parquet(proc / "master_neurons.parquet")
    pheno = pd.read_parquet(proc / "phenotypes.parquet")
    fp = pd.read_parquet(proc / "fingerprints.parquet")
    df = (master
          .merge(pheno, on=["session", "scan_idx", "unit_id"], how="inner")
          .merge(fp, on="pt_root_id", how="left"))
    return df


def main() -> int:
    # Freeze the pre-registration.
    lock = {"prereg_git_sha": _git_sha(), "materialization_version": CFG.materialization_version,
            "seed": CFG.seed}
    (CFG.path("outputs") / "prereg_lock.json").write_text(json.dumps(lock, indent=2))

    df = _assemble()
    scan_key, group_key = CFG.validation["scan_key"], CFG.validation["group_key"]
    common = dict(scan_key=scan_key, group_key=group_key, seed=CFG.seed)

    results = {}
    # Primary endpoint: state modulation. Use loco_given_pupil as the headline M_i target.
    for endpoint, y_col, excl in [
        ("state_modulation", "loco_given_pupil", False),
        ("reliability", "reliability", True),
    ]:
        sub = df.dropna(subset=[y_col])
        ladder = run_ladder(sub, y_col, endpoint=endpoint, exclude_reliability=excl,
                            n_boot=CFG.validation["n_bootstrap"], **common)
        nulls = {
            n: run_null(sub, y_col, n, endpoint=endpoint,
                        n_perm=CFG.validation["n_permutation"], exclude_reliability=excl, **common)
            for n in ["N2_layer_morph", "N3_compartment"]
        }
        results[endpoint] = {"r2": ladder.r2, "increments": ladder.increments, "nulls": nulls}

    (CFG.path("outputs") / "blinded_comparison.md").write_text(_render(results))
    _write_decision(results)
    print(open(CFG.path("outputs") / "DAY10_DECISION.md").read())
    return 0


def _render(results: dict) -> str:
    lines = ["# Blinded comparison (Days 8-9)", ""]
    for ep, r in results.items():
        lines += [f"## Endpoint: {ep}", "", "### R2_oos by model", ""]
        lines += [f"- {m}: {v:.4f}" for m, v in r["r2"].items()]
        lines += ["", "### Pre-registered increments (dR2, scan-clustered bootstrap 95% CI)", ""]
        for k, inc in r["increments"].items():
            lines.append(f"- {k}: dR2={inc['dR2']:.4f} "
                        f"[{inc['ci_lo']:.4f}, {inc['ci_hi']:.4f}] p(<=0)={inc['p_ge_0']:.3f}")
        lines += ["", "### Matched nulls", ""]
        for n, nr in r["nulls"].items():
            lines.append(f"- {n}: observed={nr['observed_dR2']:.4f} "
                        f"null_mean={nr['null_mean']:.4f} p={nr['p_value']:.3f} "
                        f"beats={'✅' if nr['beats_null'] else '❌'}")
        lines.append("")
    return "\n".join(lines)


def _write_decision(results: dict) -> None:
    def survives(ep: str) -> bool:
        r = results[ep]
        inc = r["increments"]["M5-M2"]
        beats = all(n["beats_null"] for n in r["nulls"].values())
        return bool(inc["dR2"] > 0 and inc["ci_lo"] > 0 and beats)

    state_ok = survives("state_modulation")
    reliab_ok = survives("reliability")
    if state_ok:
        verdict = "✅ GREENLIGHT full development (primary endpoint survives)."
    elif reliab_ok:
        verdict = "🟡 NARROW claim: reliability endpoint survives, state modulation does not (docs/04)."
    else:
        verdict = "❌ KILL: neither endpoint beats matched nulls out-of-sample (docs/04)."

    (CFG.path("outputs") / "DAY10_DECISION.md").write_text(
        f"# Day 10 Decision\n\n{verdict}\n\n"
        f"- state_modulation survives: {state_ok}\n"
        f"- reliability survives: {reliab_ok}\n\n"
        f"See outputs/blinded_comparison.md for the numbers behind this.\n"
    )


if __name__ == "__main__":
    sys.exit(main())
