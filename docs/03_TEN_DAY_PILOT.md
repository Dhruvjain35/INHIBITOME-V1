# 03 — Ten-Day Feasibility Pilot (runbook)

The pilot exists to answer one question fast: **is there enough clean, matched signal to justify full
development?** Each phase ends at a **gate**. If a gate fails, you stop and read `docs/04`. No phase is
"just build it and see."

Each day maps to a numbered script in `scripts/`. Scripts write to `outputs/` and cache to `data/`, and
they are **idempotent** (safe to re-run; cached CAVE queries are reused).

---

## Days 1–2 — Freeze & join · `scripts/01_freeze_and_join.py`

**Goal:** one row per unique EM neuron, with everything attached, at a pinned materialization.

1. Pin materialization **1300** (`config/pilot.yaml`); record the version + query timestamp in
   `outputs/data_provenance.json`.
2. Pull and join:
   - manual coregistration (`coregistration_manual_v4`) → the cohort;
   - cell classes (E/I, m-type) + restrict to **excitatory** post-synaptic neurons;
   - functional-area + soma depth (via `standard-transform`);
   - proofreading status + reconstruction-quality flags;
   - incoming synapses (`synapses_pni_2`, `post_ids=`) with pre-synaptic E/I labels;
   - synapse compartment predictions (`synapse_target_predictions_ssa`).
3. Emit `data/processed/master_neurons.parquet` — one row per `pt_root_id`, with list-valued or
   pre-aggregated synapse columns.

**Deliverable:** the master table + provenance JSON.

---

## Day 3 — Sample accounting · `scripts/02_sample_accounting.py`  ← **DATA GATE**

Print and save (`outputs/sample_accounting.md`):
- unique usable neurons; # scans; neurons per scan; neurons per layer/area;
- incoming inhibitory synapses per neuron (distribution);
- % of synapses with usable **compartment** labels; % with usable **presynaptic class** labels;
- missingness patterns; proofreading coverage.

**GATE — proceed only if ALL hold** (`config/pilot.yaml → gates`):
- ≥ **1,000** high-confidence coregistered excitatory neurons,
- ≥ **10** functional scans,
- **hundreds** of neurons with reliable inhibitory fingerprints,
- usable behavior across multiple scans,
- enough structural variability to support prediction.

If it fails → `docs/04` kill criteria. Do not "borrow" from the auto-coreg table to inflate counts
without re-checking quality; note any expansion explicitly.

---

## Days 4–5 — Functional targets · `scripts/03_functional_targets.py`

**Reliability `R_i`:** reproduce the oracle score from the repeated-movie NWB responses
(`phenotypes/reliability.py`); cross-check against any released score. Keep it distinct from mean
activity, selectivity, amplitude, or digital-twin performance.

**State modulation `M_i`:** fit the stimulus-aware encoding model per neuron
(`phenotypes/state_modulation.py`):
```
y_i(t) = f_i(S_t) + b_L·L_t + b_P·P_t + b_E·E_t + b_LP·(L_t·P_t) + eps
```
Evaluate on **contiguous held-out time blocks** (not shuffled samples — locomotion autocorrelates).
Derive per neuron: locomotion modulation | pupil (and vice-versa), additive shift, multiplicative gain,
stimulus-dependent modulation, and an uncertainty estimate.

**PHENOTYPE GATE:** split recordings into independent blocks A/B; require `Corr(M_i^A, M_i^B) ≥ 0.30`.
Noise-dominated neurons are dropped or uncertainty-weighted. Save `outputs/phenotype_qc.md`.

---

## Days 6–7 — Basic fingerprints · `scripts/04_fingerprints.py`

**Pilot fingerprint only** (motifs excluded — those are full-development):
- total inhibitory input; inhibitory fraction;
- perisomatic/proximal **vs** distal input (validated vs `allen_v1_column_types_slanted_ref` first);
- number of inhibitory source neurons;
- broad source-class composition.

Always normalize by dendritic size / total input — never raw counts alone. Emit
`data/processed/fingerprints.parquet` and `outputs/fingerprint_qc.md` (incl. agreement with the manual
census).

---

## Days 8–9 — First blinded comparison · `scripts/05_blinded_comparison.py`

Freeze `docs/02` (tag `prereg-frozen`, write `outputs/prereg_lock.json`). Then run, with
**leave-one-scan-out** validation, the reduced ladder:
- technical baseline (M0/M1),
- + total inhibition (M3),
- + inhibitory fingerprint (M4/M5),
- **morphology-matched shuffled fingerprint (N2)** and **compartment-location null (N3)**.

Report per-fold `dR2`, CIs (scan-clustered bootstrap), and where the observed increment sits in the null
distributions. Save `outputs/blinded_comparison.md` + `outputs/figures/`.

---

## Day 10 — Hard decision

**Greenlight full development only if ALL hold:**
- the data join is large and clean;
- reliability *or* state modulation is reproducibly measurable;
- inhibitory features add signal **beyond morphology and depth**;
- the improvement appears in a **majority of scan folds**;
- **matched nulls perform worse**;
- no direct overlapping paper has appeared since the freeze.

Write the verdict to `outputs/DAY10_DECISION.md` with the numbers that drove it. Anything short of the
above → narrow the claim (if only reliability survives) or kill (`docs/04`).

---

## Operating notes

- **Cache aggressively.** CAVE queries are slow and rate-limited; every query is cached to `data/cache/`
  keyed by (table, version, filter). Re-runs are cheap.
- **Respect the 500k row cap** on `synapses_pni_2` — chunk by id batches.
- **One mouse.** Nothing here estimates across animals. The statistical unit is the scan/block.
- **Blind yourself.** Do not look at `dR2` before the nulls are wired up; that is the whole point of the
  freeze.
