"""Days 6-7 — build inhibitory fingerprints Z_i (pilot feature set). See docs/03."""
from __future__ import annotations

import sys

import pandas as pd

from inhibitome.config import CFG
from inhibitome.fingerprints import build_fingerprints


def main() -> int:
    proc = CFG.path("processed")
    synapses = pd.read_parquet(proc / "incoming_synapses.parquet")

    # Optional per-neuron dendrite length for normalization (from skeletons; else None).
    dendrite_length = None
    dl_path = proc / "dendrite_length.parquet"
    if dl_path.exists():
        dl = pd.read_parquet(dl_path).set_index("pt_root_id")["dendrite_length_um"]
        dendrite_length = dl

    fp = build_fingerprints(synapses, dendrite_length=dendrite_length, pilot=True)
    out = proc / "fingerprints.parquet"
    fp.to_parquet(out)

    (CFG.path("outputs") / "fingerprint_qc.md").write_text(
        f"# Fingerprint QC (Days 6-7)\n\n"
        f"- Neurons with a fingerprint: {len(fp)}\n"
        f"- Median inhibitory synapses/neuron: {fp['inh_synapse_count'].median():.0f}\n"
        f"- Median perisomatic fraction: {fp['inh_frac_perisomatic'].median():.3f}\n"
        f"- Median dendritic fraction: {fp['inh_frac_dendritic'].median():.3f}\n"
    )
    print(f"Wrote {out} ({len(fp)} neurons).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
