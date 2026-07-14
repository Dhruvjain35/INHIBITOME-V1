"""Functional data loader — calcium activity, behavior, and repeated-movie stimuli from DANDI 000402.

The functional plane is NOT in CAVE (docs/01_DATA_ACCESS.md §4). It lives in NWB files on DANDI
dandiset 000402. This module downloads (once) and provides typed accessors that the phenotype code
consumes. NWB internal paths vary by file; the accessors below resolve them defensively and the exact
resolved paths are logged the first time each file is opened.

Download (one time):
    pip install dandi
    dandi download DANDI:000402   # -> data/raw/dandi_000402/
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from inhibitome.config import CFG


@dataclass
class ScanFunctional:
    """One imaging scan's aligned functional data.

    Arrays are time-aligned to `t` (seconds). `activity` is (n_units, n_time) deconvolved/df-f.
    `unit_ids` link to the coregistration table's `unit_id` for this (session, scan_idx).
    """
    session: int
    scan_idx: int
    t: np.ndarray                 # (T,)
    activity: np.ndarray          # (n_units, T)
    unit_ids: np.ndarray          # (n_units,)
    locomotion: np.ndarray        # (T,) treadmill velocity
    pupil: np.ndarray             # (T,) pupil diameter
    eye_xy: np.ndarray            # (T, 2) eye position
    stim_id: np.ndarray           # (T,) stimulus/clip id per frame
    oracle_repeats: dict          # {clip_id: list of (start,end) index ranges over the 10 repeats}


def dandi_root() -> Path:
    return (CFG.path("raw") / "dandi_000402")


def list_nwb_files() -> list[Path]:
    root = dandi_root()
    files = sorted(root.rglob("*.nwb"))
    if not files:
        raise FileNotFoundError(
            f"No NWB files under {root}. Run: `dandi download DANDI:{CFG.raw['data']['functional_source']['dandiset']}`"
        )
    return files


def load_scan(nwb_path: Path) -> ScanFunctional:
    """Open one NWB file and extract the aligned functional bundle.

    Uses pynwb. The exact NWB neurodata paths (behavior TimeSeries names, stimulus tables, the
    plane-segmentation column holding pt_root_id/unit_id) differ across files, so resolution is
    defensive and logged. This is the one place that must be validated against a real file on Day 4
    before trusting downstream phenotypes — see the assertions at the end.
    """
    from pynwb import NWBHDF5IO  # local import; heavy dependency

    with NWBHDF5IO(str(nwb_path), "r") as io:
        nwb = io.read()
        # --- these resolvers are intentionally explicit; adjust to the real schema on Day 4 ---
        raise NotImplementedError(
            "Resolve NWB internal paths against a real DANDI 000402 file on Day 4, then implement:\n"
            "  - activity: ophys ROIResponseSeries (fluorescence or deconvolved)\n"
            "  - unit_ids/pt_root_id: PlaneSegmentation columns\n"
            "  - locomotion/pupil/eye_xy: behavior module TimeSeries\n"
            "  - stim_id + oracle_repeats: stimulus presentation table (find the 10x repeated clips)\n"
            "Return a fully-populated ScanFunctional. Keep the defensive column resolution pattern "
            "used in data/join.py."
        )


def find_oracle_repeats(stim_id: np.ndarray, min_repeats: int = 8) -> dict:
    """Identify clips presented >= min_repeats times (the oracle set) and their frame ranges.

    Pure/testable given a stim_id-per-frame array. Groups contiguous runs of the same clip id into
    presentations, then keeps clips with enough presentations.
    """
    presentations: dict[int, list[tuple[int, int]]] = {}
    if len(stim_id) == 0:
        return {}
    start = 0
    for i in range(1, len(stim_id) + 1):
        if i == len(stim_id) or stim_id[i] != stim_id[start]:
            cid = int(stim_id[start])
            presentations.setdefault(cid, []).append((start, i))
            start = i
    return {cid: runs for cid, runs in presentations.items() if len(runs) >= min_repeats}
