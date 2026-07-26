"""End-to-end tests for the model ladder and the null machinery (docs/02 §2-§4).

test_core_logic.py covers the leaf functions; nothing there ever calls `run_ladder` or `run_null`,
which is how three defects reached main. These tests exercise the two entry points that actually
produce the confirmatory numbers, on synthetic data with a known planted effect.
"""
import numpy as np
import pandas as pd
import pytest

from inhibitome.fingerprints.build import build_fingerprints
from inhibitome.models.nested import FEATURE_BLOCKS, increment_point, run_ladder
from inhibitome.nulls.controls import run_null

SCAN_KEY = ["session", "scan_idx"]
GROUP_KEY = "pt_root_id"


def _synthetic(n=400, n_scans=8, effect=0.7, seed=0) -> pd.DataFrame:
    """Neuron table with a planted inh_frac_soma -> target effect, one row per neuron."""
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({
        "pt_root_id": np.arange(n),
        "session": rng.integers(0, 2, n),
        "scan_idx": rng.integers(0, n_scans // 2, n),
        "area": rng.choice(["V1", "AL"], n),
        "layer": rng.choice(["L23", "L4", "L5"], n),
        "mtype": rng.choice(["23P", "4P", "5P"], n),
    })
    for c in ["depth", "imaging_quality", "em_boundary_dist", "dendrite_length",
              "baseline_activity", "tuning", "response_amplitude", "reliability",
              "total_exc_input", "total_inh_input", "inh_source_entropy"]:
        df[c] = rng.normal(size=n)
    for c in ["inh_frac_soma", "inh_frac_proximal", "inh_frac_distal_basal", "inh_frac_apical",
              "inh_dominant_source_frac", "inh_effective_n_classes", "inh_multisynaptic_frac"]:
        df[c] = rng.random(n)
    df["n_inh_source_neurons"] = rng.integers(1, 20, n)
    df["target"] = effect * df["inh_frac_soma"] + 0.3 * df["depth"] + rng.normal(0, 0.3, n)
    return df


def test_ladder_reports_finite_increments_without_bootstrap():
    """n_boot=0 must still yield the observed dR2 — permutation tests depend on it."""
    df = _synthetic()
    res = run_ladder(df, "target", scan_key=SCAN_KEY, group_key=GROUP_KEY,
                     endpoint="state_modulation", n_boot=0, seed=1)
    assert set(res.r2) == {"M0", "M1", "M2", "M3", "M4", "M5"}
    for name, inc in res.increments.items():
        assert np.isfinite(inc["dR2"]), f"{name} dR2 is NaN with n_boot=0"


def test_planted_compartment_effect_is_recovered():
    """A real inh_frac_soma effect must show up as a positive M4-M3 increment."""
    df = _synthetic(effect=1.5)
    res = run_ladder(df, "target", scan_key=SCAN_KEY, group_key=GROUP_KEY,
                     endpoint="state_modulation", n_boot=0, seed=1)
    assert res.increments["M4-M3"]["dR2"] > 0.05


def test_exclude_reliability_does_not_leak_across_calls():
    """The reliability endpoint drops 'reliability'; the next endpoint must get it back."""
    df = _synthetic()
    before = list(FEATURE_BLOCKS["functional"])
    run_ladder(df, "target", scan_key=SCAN_KEY, group_key=GROUP_KEY,
               endpoint="reliability", exclude_reliability=True, n_boot=0, seed=1)
    assert FEATURE_BLOCKS["functional"] == before
    assert "reliability" in FEATURE_BLOCKS["functional"]


def test_exclude_reliability_actually_removes_the_target_column():
    df = _synthetic().assign(__y__=lambda d: d["target"])
    with_rel = increment_point(df, SCAN_KEY, GROUP_KEY, "M2", "M0")
    without = increment_point(df, SCAN_KEY, GROUP_KEY, "M2", "M0",
                              drop=frozenset({"reliability"}))
    assert with_rel != pytest.approx(without), "drop= had no effect on the fitted columns"


def test_bootstrap_ci_brackets_the_point_estimate():
    df = _synthetic(effect=1.5)
    res = run_ladder(df, "target", scan_key=SCAN_KEY, group_key=GROUP_KEY,
                     endpoint="state_modulation", n_boot=25, seed=1)
    inc = res.increments["M4-M3"]
    assert np.isfinite(inc["ci_lo"]) and np.isfinite(inc["ci_hi"])
    assert inc["ci_lo"] <= inc["ci_hi"]


def test_run_null_returns_a_usable_p_value():
    """The whole Day 8-9 verdict hangs on this not being NaN."""
    df = _synthetic(effect=1.5)
    r = run_null(df, "target", "N3_compartment", scan_key=SCAN_KEY, group_key=GROUP_KEY,
                 endpoint="state_modulation", n_perm=20, seed=1)
    assert np.isfinite(r["observed_dR2"])
    assert np.isfinite(r["p_value"])
    assert 0 < r["p_value"] <= 1  # +1-corrected, so never exactly 0


def test_null_destroys_a_planted_effect():
    """Shuffling compartments must cost the model the signal it had."""
    df = _synthetic(effect=1.5)
    r = run_null(df, "target", "N3_compartment", scan_key=SCAN_KEY, group_key=GROUP_KEY,
                 endpoint="state_modulation", hi="M4", lo="M3", n_perm=20, seed=1)
    assert r["observed_dR2"] > r["null_mean"]
    assert r["beats_null"]


def test_pure_noise_target_does_not_beat_the_null():
    """No planted effect -> the fingerprint must not look better than matched chance wiring."""
    df = _synthetic(effect=0.0)
    r = run_null(df, "target", "N3_compartment", scan_key=SCAN_KEY, group_key=GROUP_KEY,
                 endpoint="state_modulation", hi="M4", lo="M3", n_perm=30, seed=2)
    assert not r["beats_null"]


# --- M5 must be a real model, not a silent alias for M4 (docs/02 §2) ----------------------------

def _syn_table() -> pd.DataFrame:
    """Two neurons: 200 gets diverse inhibitory sources, 201 gets a single source."""
    return pd.DataFrame({
        "post_pt_root_id": [200] * 4 + [201] * 4,
        "pre_pt_root_id": [1, 2, 3, 4, 9, 9, 9, 9],
        "pre_ei": ["inhibitory"] * 8,
        "pre_mtype": ["BC", "MC", "BPC", "NGC", "BC", "BC", "BC", "BC"],
        "compartment": ["soma", "soma", "distal_basal", "apical"] * 2,
    })


def test_pilot_fingerprints_include_the_m5_source_block():
    fp = build_fingerprints(_syn_table(), pilot=True)
    for c in ["inh_source_entropy", "inh_dominant_source_frac", "inh_effective_n_classes"]:
        assert c in fp.columns, f"{c} missing -> M5 silently collapses onto M4"


def test_source_diversity_separates_diverse_from_single_source_neurons():
    fp = build_fingerprints(_syn_table(), pilot=True).set_index("pt_root_id")
    assert fp.loc[200, "inh_source_entropy"] > fp.loc[201, "inh_source_entropy"]
    assert fp.loc[200, "inh_effective_n_classes"] == pytest.approx(4.0)
    assert fp.loc[201, "inh_effective_n_classes"] == pytest.approx(1.0)
    assert fp.loc[201, "inh_dominant_source_frac"] == pytest.approx(1.0)


def test_coarse_only_labels_give_nan_not_a_fake_zero():
    """Without pre_mtype, source diversity is undefined -- must not read as 'no signal'."""
    syn = _syn_table().drop(columns=["pre_mtype"])
    fp = build_fingerprints(syn, pilot=True)
    assert fp["inh_source_entropy"].isna().all()
    assert fp["inh_effective_n_classes"].isna().all()
