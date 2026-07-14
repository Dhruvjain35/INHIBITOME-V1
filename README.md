# INHIBITOME-V1

**How inhibitory synaptic architecture shapes stable and state-dependent visual computation.**

> **One-sentence question.** Does the amount, cellular source, and subcellular placement of
> inhibitory synapses onto an individual visual-cortex neuron predict both how *reliably* it
> responds to visual stimuli and how strongly locomotion / pupil-linked arousal *modulates* that
> response — over and above the neuron's layer, area, depth, morphology, and tuning?

This repository is the full research plan **and** the executable pilot for that question, built on the
[MICrONS](https://www.microns-explorer.org/) cubic-millimeter mouse visual-cortex dataset (same-animal
electron microscopy + two-photon calcium imaging + behavior).

---

## The honest framing

This is a **10-day feasibility pilot, not an unconditional full build**. The upside is legitimately
ISEF- and specialist-publication-worthy; the central risk is that *static* inhibitory anatomy adds
little predictive signal once you control for layer, morphology, depth, and visual tuning. The whole
design is built to find that out fast and kill cleanly if the signal isn't there.

| Dimension | Rating |
|---|---|
| Scientific importance | 9/10 |
| Dataset quality | 9/10 |
| Prob. of a *complete, defensible* study | 8/10 |
| Prob. inhibitory wiring adds a *strong independent* effect | 6.5–7/10 |
| ISEF ceiling with a decisive result | 8.5–9/10 |
| Largest limitation | **One mouse** |

The single-animal design means every claim is *"a principle observed within the MICrONS specimen,"*
never a universal law of mouse cortex. That framing is load-bearing and appears in every output.

---

## Repository map

```
docs/                         The plan. Read these in order.
  00_PROJECT_PLAN.md          Full scientific framework, aims, novelty positioning.
  01_DATA_ACCESS.md           How to get every MICrONS table/asset used here.
  02_PREREGISTRATION.md       Frozen hypotheses, endpoints, model hierarchy, null models.
  03_TEN_DAY_PILOT.md         Day-by-day runbook with hard go/no-go gates.
  04_KILL_AND_SUCCESS.md      What counts as success; what triggers a kill.
  05_ISEF_AND_PUBLICATION.md  Category, narrative, figures, venues.
config/
  pilot.yaml                  FROZEN parameters: materialization version, table names, gates.
src/inhibitome/               The pipeline package (data -> phenotypes -> fingerprints -> models -> nulls).
scripts/                      Numbered, runnable pilot steps mirroring the 10-day runbook.
notebooks/                    Exploration only; nothing load-bearing lives here.
outputs/                      Generated tables, figures, and the sample-accounting report (gitignored).
data/                         Local data cache (gitignored; fully reproducible from the pipeline).
```

## Quick start

```bash
# 1. Environment (uv is fast and reproducible; conda also fine)
uv sync --extra dev

# 2. Configure MICrONS / CAVE access (see docs/01_DATA_ACCESS.md)
#    Sets up the CAVE auth token once.
uv run python scripts/00_setup_cave_token.py

# 3. Run the pilot, one gate at a time (see docs/03_TEN_DAY_PILOT.md)
uv run python scripts/01_freeze_and_join.py       # Days 1-2: pin data, build master table
uv run python scripts/02_sample_accounting.py     # Day 3:   DATA GATE - proceed only if it passes
uv run python scripts/03_functional_targets.py    # Days 4-5: reliability + state-modulation targets
uv run python scripts/04_fingerprints.py          # Days 6-7: inhibitory fingerprints
uv run python scripts/05_blinded_comparison.py    # Days 8-9: nested models + matched nulls
#                                                    Day 10:  HARD DECISION (see docs/04)
```

## The core statistical claim (what we actually test)

For each functionally-coregistered excitatory neuron `i`, we test whether an **inhibitory fingerprint**
`Zi = [amount, source, location, diversity, motifs]` improves *out-of-sample* prediction of:

- **Ri** — trial-to-trial visual **response reliability** (oracle score on repeated movies), and
- **Mi** — locomotion/pupil-linked **state modulation** of visual responses,

**after** a nested baseline (technical -> cellular -> functional -> total-input). The headline number is
the held-out increment `dR2`, evaluated with **leave-one-scan-out** validation and benchmarked against
**layer/morphology-matched shuffled fingerprints**. A model that fits training data but does not beat
matched nulls out-of-sample is **not evidence**.

See `docs/02_PREREGISTRATION.md` for the frozen hierarchy `M0...M5`, the seven null models, and the
confirmatory endpoints.

## Status

Pilot scaffold. Data-access strings in `config/pilot.yaml` are pinned to a specific MICrONS
materialization for reproducibility — see `docs/01_DATA_ACCESS.md` before changing them.
