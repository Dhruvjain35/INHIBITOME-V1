"""Day 3 — sample accounting + DATA GATE. See docs/03 and docs/04."""
from __future__ import annotations

import sys

import pandas as pd

from inhibitome.config import CFG
from inhibitome.report import sample_accounting


def main() -> int:
    proc = CFG.path("processed")
    master = pd.read_parquet(proc / "master_neurons.parquet")
    synapses = pd.read_parquet(proc / "incoming_synapses.parquet")

    result = sample_accounting(master, synapses)
    print(open(result["report"]).read())

    if not result["passed"]:
        print("\n❌ DATA GATE FAILED — see docs/04_KILL_AND_SUCCESS.md before continuing.")
        return 2
    print("\n✅ DATA GATE PASSED — proceed to scripts/03_functional_targets.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
