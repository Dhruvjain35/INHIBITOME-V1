"""Day-3 sample accounting + the hard DATA GATE (docs/03, config.gates).

Produces outputs/sample_accounting.md and returns a pass/fail dict. The gate decides whether the
pilot continues at all.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from inhibitome.config import CFG


def sample_accounting(master: pd.DataFrame, synapses: pd.DataFrame,
                     out_path: Path | None = None) -> dict:
    g = CFG.gates
    scan_key = CFG.validation["scan_key"]
    scan_id = master[scan_key].astype(str).agg("|".join, axis=1) if set(scan_key).issubset(
        master.columns) else pd.Series(["?"] * len(master))

    n_neurons = master["pt_root_id"].nunique()
    n_scans = scan_id.nunique()
    neurons_per_scan = master.groupby(scan_id)["pt_root_id"].nunique()

    inh = synapses[
        synapses["pre_ei"].astype(str).str.lower().str.startswith("inhibitory")
    ] if "pre_ei" in synapses.columns else synapses.iloc[0:0]
    inh_per_neuron = inh.groupby("post_pt_root_id").size() if not inh.empty else pd.Series(dtype=int)
    # A fingerprint is only as good as the number of typed inhibitory inputs behind it. Measured
    # median is 129/neuron (min 26), so this floor drops the tail rather than most of the cohort.
    min_inh = g["min_typed_inhibitory_per_neuron"]
    n_with_fingerprint = int((inh_per_neuron >= min_inh).sum())

    comp_labeled = _labeled_fraction(synapses, "compartment")
    class_labeled = _labeled_fraction(synapses, "pre_ei")
    frac_typed = (float(master["frac_typed"].median())
                  if "frac_typed" in master.columns else float("nan"))

    checks = {
        "min_coreg_neurons": (n_neurons, g["min_coreg_neurons"], n_neurons >= g["min_coreg_neurons"]),
        "min_scans": (n_scans, g["min_scans"], n_scans >= g["min_scans"]),
        "min_neurons_with_fingerprint": (
            n_with_fingerprint, g["min_neurons_with_fingerprint"],
            n_with_fingerprint >= g["min_neurons_with_fingerprint"]),
        "min_compartment_labeled_fraction": (
            round(comp_labeled, 3), g["min_compartment_labeled_fraction"],
            comp_labeled >= g["min_compartment_labeled_fraction"]),
    }
    passed = all(ok for _, _, ok in checks.values())

    md = _render(n_neurons, n_scans, neurons_per_scan, inh_per_neuron, comp_labeled,
                class_labeled, checks, passed, frac_typed, min_inh)
    out_path = out_path or (CFG.path("outputs") / "sample_accounting.md")
    out_path.write_text(md)
    return {"passed": passed, "checks": checks, "report": str(out_path)}


def _labeled_fraction(df: pd.DataFrame, col: str) -> float:
    if col not in df.columns or len(df) == 0:
        return 0.0
    s = df[col].astype(str).str.lower()
    return float((~s.isin(["", "nan", "none", "unknown"])).mean())


def _render(n_neurons, n_scans, per_scan, inh_per_neuron, comp, cls, checks, passed,
            frac_typed=float("nan"), min_inh=0) -> str:
    def q(s, f):
        return int(s.quantile(f)) if len(s) else 0

    lines = ["# Sample accounting (Day 3)", ""]
    lines += [f"- Materialization: **{CFG.materialization_version}**",
              f"- Unique usable EM neurons: **{n_neurons}**",
              f"- Functional scans: **{n_scans}**",
              f"- Neurons/scan: min {int(per_scan.min()) if len(per_scan) else 0}, "
              f"median {int(per_scan.median()) if len(per_scan) else 0}, "
              f"max {int(per_scan.max()) if len(per_scan) else 0}",
              f"- Typed inhibitory synapses/neuron: min {q(inh_per_neuron, 0)}, "
              f"q25 {q(inh_per_neuron, .25)}, median {q(inh_per_neuron, .5)}, "
              f"q75 {q(inh_per_neuron, .75)}, max {q(inh_per_neuron, 1)} "
              f"(fingerprint floor: {min_inh})",
              f"- Median `frac_typed` (reconstruction completeness): **{frac_typed:.1%}**",
              f"- Synapses with compartment label: **{comp:.1%}**",
              f"- Synapses with presynaptic class: **{cls:.1%}**",
              "",
              "> `frac_typed` is low by construction — ~91% of incoming synapses come from orphan",
              "> axon fragments with no soma in the volume. That is a property of the EM",
              "> reconstruction, not of this cohort. It is why all fractions are normalized over",
              "> typed inputs and why `frac_typed` is an M0 covariate. See docs/01 §3.",
              ""]
    lines += ["## DATA GATE", "", "| check | value | threshold | pass |", "|---|---|---|---|"]
    for name, (val, thr, ok) in checks.items():
        lines.append(f"| {name} | {val} | {thr} | {'✅' if ok else '❌'} |")
    lines += ["", f"### Verdict: {'✅ PROCEED' if passed else '❌ STOP — see docs/04'}"]
    return "\n".join(lines)
