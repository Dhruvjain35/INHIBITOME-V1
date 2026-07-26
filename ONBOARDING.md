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
tests/       test_core_logic.py      leaf functions (oracle, encoding model, fingerprints, CV)
             test_ladder_and_nulls.py  run_ladder / run_null end-to-end + the silent-zero traps
```
23/23 passing, no network required.

## Get running

```bash
make setup            # uv-managed venv (Python 3.10–3.12; the heavy connectomics stack needs <3.13)
make token            # one-time CAVE auth (accept ToS + mint token; see docs/01 §1)
make test             # sanity-check the pure logic offline
make pilot            # runs scripts 01→05 in order, stopping at the first failed gate
```

No `make` on Windows? Every target is a one-liner — `uv sync --extra dev`, `uv run pytest -q`,
`uv run python scripts/01_freeze_and_join.py`, and so on. `make token` is interactive and needs a
real terminal; the token lands in `~/.cloudvolume/secrets/cave-secret.json`.

The Days 1–2 pull transfers ~47M synapse rows and takes a while. It caches per batch, so killing it
and re-running resumes rather than restarting. `make clean` drops the cache and forces a refetch.

## The things that will bite you

1. **Function is NOT in CAVE.** Anatomy (synapses, cell types, compartments, coregistration) comes from
   CAVE `minnie65_public`. Activity + behavior + repeated-movie ("oracle") stimuli come from **DANDI
   dandiset 000402** (NWB). The coregistration table bridges them. Download NWB with `dandi download
   DANDI:000402`.
2. **The pipeline intentionally stops with `NotImplementedError` in `data/functional.py` and
   `scripts/03`.** Resolving the NWB internal schema (which TimeSeries hold locomotion/pupil, where
   `pt_root_id`/`unit_id` live, how the 10×-repeated clips are marked) must be done against a **real
   DANDI file on Day 4** — it can't be guessed from docs. That is the current frontier, not a bug.
3. **91.4% of a neuron's incoming synapses come from orphan axon fragments with no soma.** They can
   never be typed, so every fraction here is normalized over *typed* inputs, and `frac_typed` rides
   along in the M0 technical block. If you ever see an inhibitory fraction near 0.05, you are looking
   at the wrong denominator — the typed one gives ~0.6. See `fingerprints/build.py`.
4. **`pt_root_id` cannot be filtered server-side** on the cell-type / coregistration tables — it lives
   on the referenced nucleus table, and CAVE 500s with `KeyError: 'pt_root_id'`. `Cave.query_table`
   filters client-side for exactly this reason.
5. **A missing feature column does not raise.** `_columns_for` skips columns that aren't present, so an
   absent feature block silently collapses one rung of the ladder onto the one below and reads as a
   clean negative result. Three separate bugs of this shape have already been fixed here; if a model
   comparison comes back suspiciously close to zero, check the columns actually exist first.

## Current status & next step

**Days 1–2 are done and the Day-3 DATA GATE PASSES.** 27 unit tests green. Config verified against
the live API — all 16 tables exist at materialization **1822**.

| gate check | value | threshold |
|---|---|---|
| coregistered excitatory neurons | **15,282** | ≥ 1,000 |
| functional scans | **16** | ≥ 10 |
| neurons with a fingerprint | **15,223** | ≥ 300 |
| compartment labels | **98.8%** | ≥ 50% |

Master table: 19,004 ROI rows / 15,282 neurons (V1 13,059 · RL 4,710 · AL 1,211 · LM 24). Pulled
4,925,938 typed synapses out of 52,347,211 total. Typed inhibitory input per neuron: median 181,
q25 120, q75 260, max 4,502. Note **LM has only 24 ROIs** — too few for area-stratified claims.

### To pick this up

1. `uv sync --extra dev`, then mint a CAVE token (`scripts/00_setup_cave_token.py`, needs a real
   terminal — it prompts).
2. **`data/` is gitignored, so the pull does not come with the repo.** Re-running
   `scripts/01_freeze_and_join.py` takes ~4h against CAVE and is resumable per batch. Copying a
   collaborator's `data/cache/` directory across is far faster and is exactly equivalent —
   everything in it is keyed by materialization version.
3. `scripts/02_sample_accounting.py` should reproduce the table above.

### Next real work: Days 4–5 (`scripts/03`)

Needs `dandi download DANDI:000402` and the NWB internal schema resolved against a real file —
`data/functional.py` raises `NotImplementedError` until then. That is the last genuine unknown in
the pipeline.

### Open decisions, deliberately not made

- **`synapse_target_predictions_ssa_v2` is flagged by its owners as "Table in development…
  Citation forthcoming, reach out if wanting to use for publication."** It is the basis of M4. The
  v1 table (`synapse_target_predictions_ssa`, 204.3M rows) is the published one behind the 2025
  census and also exists at 1822, so switching costs nothing on the version pin. For a
  publication target, v1 as primary and v2 as a robustness check is the safer shape.
- **`test_retest` correlates a 3-element coefficient vector per neuron.** A correlation over three
  points is close to meaningless; the docs/02 §6 gate (`Corr ≥ 0.30`) more plausibly means
  correlating `M_i` *across neurons* between time blocks.
- **`scripts/05` uses `loco_given_pupil` as the headline `M_i`**, but docs/02 §1 defines the primary
  endpoint as state-dependent visual *gain* — that is `multiplicative_gain`. The prereg and the code
  currently disagree about the primary endpoint.

All three are cheap to change now and awkward after the `prereg-frozen` tag (which does not yet
exist — nothing is actually frozen).

## Ground rules (from the pre-registration)

- Evidence = **held-out** `ΔR²` under **leave-one-scan-out**, beating **layer/morphology-matched null
  shuffles**. Training fit and random splits are not evidence.
- Interpretable models first (ridge/GAM); gradient-boosting is a *control* only; no GNNs/Transformers as
  primary evidence.
- The statistical unit is the **scan/block**, not the synapse — never quote a tiny p-value off millions
  of synapses from one animal.
- Report negative results and all exclusions. Reproducibility bar in `docs/05`.
