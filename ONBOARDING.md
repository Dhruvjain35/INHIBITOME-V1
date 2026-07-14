# Onboarding — INHIBITOME-V1

Welcome. This repo is a **research plan + a runnable feasibility pilot** asking one question about the
MICrONS mouse-visual-cortex dataset:

> Does the amount, cellular **source**, and subcellular **placement** of inhibitory synapses onto a
> single excitatory neuron predict how *reliably* it responds to vision (`R_i`) and how strongly
> locomotion/pupil arousal *modulates* that response (`M_i`) — beyond layer, area, depth, morphology,
> and tuning?

Read the honest framing first: this is a **10-day feasibility pilot with hard go/no-go gates**, not an
unconditional build. It's designed to *kill fast* if static inhibitory anatomy adds nothing past
morphology. The governing limitation is that all matched EM+function data come from **one mouse** — every
claim is "a principle within the MICrONS specimen," never a law of cortex.

## Read these, in order (15 min)

1. `README.md` — the elevator version + repo map.
2. `docs/00_PROJECT_PLAN.md` — question, three aims, novelty line.
3. `docs/02_PREREGISTRATION.md` — the frozen `M0→M5` model ladder, seven null models, endpoints. This is
   the scientific contract; don't change it after the `prereg-frozen` tag.
4. `docs/03_TEN_DAY_PILOT.md` — the day-by-day runbook with the exact gates.
5. `docs/01_DATA_ACCESS.md` — how to actually get the data (skim now, live-reference later).

`docs/04` (kill/success criteria) and `docs/05` (ISEF + publication) are reference.

## How the code is laid out

```
config/pilot.yaml          Every reproducibility-critical string is FROZEN here (tables, version, gates,
                           seed). Never hard-code a table name in a module — read it from CFG.
src/inhibitome/
  data/        CAVE wrapper (pinned materialization + query cache) · master join · DANDI functional loader
  phenotypes/  reliability.py (oracle R_i) · state_modulation.py (stimulus-aware M_i + split-half gate)
  fingerprints/ inhibitory amount/location/source/diversity features (always normalized)
  models/      nested.py (M0–M5 ladder + bootstrap ΔR²) · validation.py (grouped leave-one-scan-out)
  nulls/       matched-permutation / compartment / presynaptic-identity controls
  report/      Day-3 sample accounting + the DATA GATE
scripts/00–05              Numbered, runnable pilot steps that mirror the runbook exactly.
tests/                     Pure-logic unit tests (no network). 7/7 passing.
```

## Get running

```bash
make setup            # uv-managed venv (Python 3.10–3.12; the heavy connectomics stack needs <3.13)
make token            # one-time CAVE auth (accept ToS + mint token; see docs/01 §1)
make test             # sanity-check the pure logic offline
make pilot            # runs scripts 01→05 in order, stopping at the first failed gate
```

## The two things that will bite you

1. **Function is NOT in CAVE.** Anatomy (synapses, cell types, compartments, coregistration) comes from
   CAVE `minnie65_public`. Activity + behavior + repeated-movie ("oracle") stimuli come from **DANDI
   dandiset 000402** (NWB). The coregistration table bridges them. Download NWB with `dandi download
   DANDI:000402`.
2. **The pipeline intentionally stops with `NotImplementedError` in `data/functional.py` and
   `scripts/03`.** Resolving the NWB internal schema (which TimeSeries hold locomotion/pupil, where
   `pt_root_id`/`unit_id` live, how the 10×-repeated clips are marked) must be done against a **real
   DANDI file on Day 4** — it can't be guessed from docs. That is the current frontier, not a bug.

## Current status & next step

Scaffold complete; plan frozen in `docs/`; unit tests green. **Nothing has touched live data yet.**
Next real work is **Days 1–2**: `make token`, then `python scripts/01_freeze_and_join.py` to build the
master neuron table, then the **Day-3 DATA GATE** (`scripts/02`) decides whether the pilot proceeds at
all. If a gate fails, go straight to `docs/04_KILL_AND_SUCCESS.md` — a clean kill is a legitimate outcome.

## Ground rules (from the pre-registration)

- Evidence = **held-out** `ΔR²` under **leave-one-scan-out**, beating **layer/morphology-matched null
  shuffles**. Training fit and random splits are not evidence.
- Interpretable models first (ridge/GAM); gradient-boosting is a *control* only; no GNNs/Transformers as
  primary evidence.
- The statistical unit is the **scan/block**, not the synapse — never quote a tiny p-value off millions
  of synapses from one animal.
- Report negative results and all exclusions. Reproducibility bar in `docs/05`.
