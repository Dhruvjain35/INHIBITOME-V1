"""Days 1-2 — pin materialization and build the master neuron table. See docs/03."""
from __future__ import annotations

import json
import sys

from inhibitome.config import CFG
from inhibitome.data.cave import Cave
from inhibitome.data.join import build_master


def main() -> int:
    cave = Cave()
    print(f"Pinned datastack={CFG.datastack} version={CFG.materialization_version}")

    provenance = {
        "datastack": CFG.datastack,
        "materialization_version": CFG.materialization_version,
        "tables": CFG.tables,
    }
    try:
        provenance["query_timestamp"] = str(
            cave.client.materialize.get_timestamp(CFG.materialization_version)
        )
    except Exception as e:  # noqa: BLE001
        provenance["query_timestamp_error"] = str(e)

    (CFG.path("outputs") / "data_provenance.json").write_text(json.dumps(provenance, indent=2))

    paths = build_master(cave)
    print("Wrote:")
    for k, v in paths.items():
        print(f"  {k}: {v}")
    print("\nNext: python scripts/02_sample_accounting.py  (Day 3 DATA GATE)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
