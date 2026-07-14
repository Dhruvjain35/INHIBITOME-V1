"""Unit tests for the pure, network-free logic: oracle reliability, encoding-model summaries,
fingerprint aggregation, grouped/leave-one-scan-out validation, and matched-null invariants.

These are the pieces whose correctness must not depend on live CAVE/DANDI access. Run: `pytest -q`.
"""
import numpy as np
import pandas as pd
import pytest

from inhibitome.fingerprints.build import build_fingerprints
from inhibitome.models.validation import leave_one_scan_out
from inhibitome.nulls.controls import compartment_location_null, matched_permutation_null
from inhibitome.phenotypes.reliability import oracle_score
from inhibitome.phenotypes.state_modulation import fit_state_modulation
from inhibitome.data.functional import find_oracle_repeats


def test_oracle_score_perfect_and_noise():
    base = np.sin(np.linspace(0, 6, 100))
    perfect = np.tile(base, (10, 1))
    assert oracle_score(perfect) > 0.99
    rng = np.random.default_rng(0)
    noise = rng.normal(size=(10, 100))
    assert abs(oracle_score(noise)) < 0.3  # unreliable -> ~0
    assert np.isnan(oracle_score(base[None, :]))  # <2 repeats


def test_find_oracle_repeats_groups_runs():
    stim = np.array([1, 1, 1, 2, 2, 1, 1, 1, 2, 2] * 1)  # clip1 x2 presentations, clip2 x2
    reps = find_oracle_repeats(stim, min_repeats=2)
    assert set(reps) == {1, 2}
    assert len(reps[1]) == 2 and len(reps[2]) == 2


def test_state_modulation_recovers_injected_locomotion_gain():
    rng = np.random.default_rng(1)
    T = 2000
    loco = rng.normal(size=T)
    pupil = rng.normal(size=T)
    eye = rng.normal(size=(T, 2))
    stim = rng.normal(size=T)
    y = stim + 0.8 * loco + 0.1 * pupil + 0.05 * rng.normal(size=T)
    m = fit_state_modulation(y, stim, loco, pupil, eye, seed=1)
    assert m.loco_given_pupil > m.pupil_given_loco  # locomotion is the injected driver
    assert m.r2_behavior_oos > 0.5


def test_build_fingerprints_normalizes_and_labels():
    # 3 inhibitory (soma-heavy) + 2 excitatory synapses onto neuron 100.
    syn = pd.DataFrame({
        "post_pt_root_id": [100] * 5,
        "pre_pt_root_id": [1, 1, 2, 3, 4],
        "pre_ei": ["inhibitory", "inhibitory", "inhibitory", "excitatory", "excitatory"],
        "compartment": ["soma", "soma", "distal_basal", "spine", "spine"],
    })
    fp = build_fingerprints(syn, pilot=True)
    row = fp.iloc[0]
    assert row["inh_synapse_count"] == 3
    assert row["inh_fraction"] == pytest.approx(3 / 5)
    assert row["n_inh_source_neurons"] == 2  # neurons 1 and 2
    assert row["inh_frac_perisomatic"] == pytest.approx(2 / 3)
    assert row["inh_frac_dendritic"] == pytest.approx(1 / 3)


def test_leave_one_scan_out_no_group_leakage():
    df = pd.DataFrame({
        "pt_root_id": [1, 1, 2, 3, 4, 5],
        "session": [1, 1, 1, 1, 2, 2],
        "scan_idx": [0, 0, 0, 0, 0, 0],
    })
    folds = list(leave_one_scan_out(df, ["session", "scan_idx"], "pt_root_id"))
    assert len(folds) == 2
    for f in folds:
        tr = set(df.iloc[f.train_idx]["pt_root_id"])
        te = set(df.iloc[f.test_idx]["pt_root_id"])
        assert tr.isdisjoint(te)  # invariant #1: no neuron in both splits


def test_compartment_null_preserves_total_count():
    df = pd.DataFrame({
        "inh_frac_soma": [0.5], "inh_frac_proximal": [0.2],
        "inh_frac_distal_basal": [0.2], "inh_frac_apical": [0.1],
    })
    shuffled = compartment_location_null(df, seed=0)
    orig_cols = ["inh_frac_soma", "inh_frac_proximal", "inh_frac_distal_basal", "inh_frac_apical"]
    assert shuffled[orig_cols].sum(axis=1).iloc[0] == pytest.approx(1.0)


def test_matched_permutation_keeps_rows_within_strata():
    df = pd.DataFrame({
        "area": ["V1", "V1", "AL", "AL"],
        "layer": ["L23", "L23", "L5", "L5"],
        "mtype": ["a", "a", "b", "b"],
        "inh_frac_soma": [0.1, 0.9, 0.2, 0.8],
    })
    out = matched_permutation_null(df, cols=["inh_frac_soma"], seed=3)
    # Values only permute within the (V1,L23,a) and (AL,L5,b) strata -> multiset preserved per stratum.
    assert sorted(out.iloc[:2]["inh_frac_soma"]) == pytest.approx([0.1, 0.9])
    assert sorted(out.iloc[2:]["inh_frac_soma"]) == pytest.approx([0.2, 0.8])
