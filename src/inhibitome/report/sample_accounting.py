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

    inh = synapses[synapses.get("pre_ei", pd.Series(dtype=str)).astype(str).str.lower().str.startswith(
        ("inh", "i"))] if "pre_ei" in synapses.columns else synapses.iloc[0:0]
    inh_per_neuron = inh.groupby("post_pt_root_id").size() if not inh.empty else pd.Series(dtype=int)
    n_with_fingerprint = int((inh_per_neuron >= 5).sum())  # "reliable fingerprint" ~ >=5 inh synapses

    comp_labeled = _labeled_fraction(synapses, "compartment")
    class_labeled = _labeled_fraction(synapses, "pre_ei")

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
                class_labeled, checks, passed)
    out_path = out_path or (CFG.path("outputs") / "sample_accounting.md")
    out_path.write_text(md)
    return {"passed": passed, "checks": checks, "report": str(out_path)}


def _labeled_fraction(df: pd.DataFrame, col: str) -> float:
    if col not in df.columns or len(df) == 0:
        return 0.0
    s = df[col].astype(str).str.lower()
    return float((~s.isin(["", "nan", "none", "unknown"])).mean())


def _render(n_neurons, n_scans, per_scan, inh_per_neuron, comp, cls, checks, passed) -> str:
    lines = ["# Sample accounting (Day 3)", ""]
    lines += [f"- Materialization: **{CFG.materialization_version}**",
              f"- Unique usable EM neurons: **{n_neurons}**",
              f"- Functional scans: **{n_scans}**",
              f"- Neurons/scan: min {int(per_scan.min()) if len(per_scan) else 0}, "
              f"median {int(per_scan.median()) if len(per_scan) else 0}, "
              f"max {int(per_scan.max()) if len(per_scan) else 0}",
              f"- Inhibitory synapses/neuron: median "
              f"{int(inh_per_neuron.median()) if len(inh_per_neuron) else 0}",
              f"- Synapses with compartment label: **{comp:.1%}**",
              f"- Synapses with presynaptic class: **{cls:.1%}**", ""]
    lines += ["## DATA GATE", "", "| check | value | threshold | pass |", "|---|---|---|---|"]
    for name, (val, thr, ok) in checks.items():
        lines.append(f"| {name} | {val} | {thr} | {'✅' if ok else '❌'} |")
    lines += ["", f"### Verdict: {'✅ PROCEED' if passed else '❌ STOP — see docs/04'}"]
    return "\n".join(lines)
