# 00 — Project Plan

*How inhibitory synaptic architecture shapes stable and state-dependent visual computation.*

---

## 1. The question, precisely

Does the **amount**, **cellular source**, and **subcellular placement** of inhibitory synapses onto an
individual visual-cortex neuron predict both:

1. how **reliably** it responds to repeated visual stimuli (a *stable* phenotype), and
2. how strongly locomotion / pupil-linked arousal **modulates** its visual response (a *flexible* phenotype),

**after controlling for** cortical area, layer, depth, imaging session, excitatory morphological type,
dendritic size, baseline activity, visual tuning, total synaptic input, and imaging quality?

The arrows we test are *prediction and statistical association*, not proven causation.

---

## 2. Why the question matters

Cortex must do two things that pull against each other:

- **Reliability** — the same visual event should produce a usable, repeatable representation.
- **Flexibility** — behavioral state (locomotion, arousal, attention) should be able to reshape how that
  same input is processed.

Inhibition is the leading candidate mechanism for negotiating that trade-off: it gates whether a neuron
responds, sets response gain, shapes trial-to-trial variability, controls dendritic integration, and
carries much of the state-dependent modulation of cortical activity.

But **"inhibition" is not one number.** Inhibitory interneurons differ systematically in *which*
excitatory cells they target and *where* on those cells they synapse — perisomatic vs. proximal-dendritic
vs. distal-basal vs. apical. The MICrONS inhibitory census demonstrated that this targeting is highly
specific and compartment-organized. The **open question** is whether that anatomical specificity, read
off a single static reconstruction, predicts what the neuron actually *does* during vision and during
changes in behavioral state.

---

## 3. Why MICrONS makes this possible

MICrONS recorded two-photon calcium activity from ~75,000 excitatory neurons across primary and higher
visual cortex in an awake mouse viewing natural and synthetic stimuli, **while also recording treadmill
locomotion, eye position, and pupil diameter.** The *same* tissue was then reconstructed by electron
microscopy — >200,000 cells and ~0.5 billion synapses.

The decisive property: **anatomy and function come from the same neurons in the same animal.** This
sidesteps the false-correspondence problem that sinks cross-animal designs (e.g. matching a functional
cell in one zebrafish to an anatomical cell in another).

Assets we lean on (exact table names in `docs/01_DATA_ACCESS.md`, pinned in `config/pilot.yaml`):

- A **manually-verified functional coregistration** (thousands of EM root IDs matched to functional
  ROIs by human curators) — this is our **primary study cohort**, not the smaller census.
- **Repeated movie clips** presented multiple times per session, enabling an **oracle** measure of
  trial-to-trial response reliability.
- **Dataset-wide cell-type predictions** (exc/inh, predicted morphological type) plus a 2025 release of
  **synapse-target / compartment predictions** for ~200M synapses.
- The **manually-curated ~1,352-cell inhibitory census** (+70k inhibitory synapses) — used as a
  **high-quality validation subset**, not the only cohort.

---

## 4. The variables

For each functionally-coregistered excitatory neuron `i`:

**Structural inhibitory fingerprint** — *how the neuron is inhibited.*
```
Z_i = [ I_amount, I_source, I_location, I_diversity, I_motifs ]
```

**Stable functional phenotype** — *trial-to-trial visual response reliability.*
```
R_i
```

**Flexible functional phenotype** — *locomotion/pupil-linked modulation of visual responses.*
```
M_i
```

We test whether `Z_i -> R_i` and/or `Z_i -> M_i` survive the nested controls in §2 of
`docs/02_PREREGISTRATION.md`.

---

## 5. Hypotheses

**Primary (state modulation).** The composition and subcellular placement of inhibitory input improves
out-of-sample prediction of a neuron's behavioral-state modulation beyond the full baseline. The
headline statistic is the held-out increment:
```
dR2_state = R2(baseline + inhibitory fingerprint) - R2(baseline)
```

**Secondary (reliability).** Inhibitory wiring explains variation in how consistently neurons respond
across repeated presentations. This is a **separate, pre-registered endpoint** — not a fallback invented
after seeing the state-modulation result.

**Mechanistic (directional).** Different inhibitory features map to different functional roles:
- perisomatic / proximal inhibition → overall reliability / global gain;
- distal / apical inhibition → stimulus-specific state modulation;
- greater source diversity → more flexible modulation across states.

These are hypotheses to test, not assumptions.

---

## 6. Novelty: what is already known, and what is not

Already shown by MICrONS-adjacent work (so **not** claimable as novel here):
- inhibitory wiring is highly specific and compartment-organized;
- excitatory connectivity relates to functional similarity;
- functional-model representations can predict anatomical properties;
- arousal changes inferred functional circuitry.

A 2026 preprint analyzed *state-dependent directed functional networks* across >57,000 MICrONS neurons
and reported that synapse count predicted inferred functional connection *strength*, with weaker
structure–function coupling during high arousal. This raises novelty pressure but does **not** perform
our analysis — it studies pairwise *inferred functional* connections and broad structural coupling, not
the single-neuron inhibitory-input fingerprint.

**Our defensible novelty:**
> Testing whether the *complete inhibitory-input fingerprint* of an individual excitatory neuron —
> including presynaptic inhibitory identity **and** postsynaptic compartment targeting — predicts that
> neuron's stable reliability and state-dependent visual gain.

The project must stay centered on this neuron-level, inhibition-specific, compartment-specific question.
Any drift toward "structure relates to function" in general is a novelty failure.

---

## 7. The three aims

### Aim 1 — Trustworthy functional phenotypes
Build `R_i` (oracle reliability on repeated movies, reproduced independently from raw responses for QC)
and `M_i` (state modulation from a **stimulus-aware** encoding model that separates locomotion, pupil,
and eye position, with additive and multiplicative terms). **Every `M_i` must pass a test–retest
reproducibility check** (split-half correlation across time blocks); noise-dominated estimates are
dropped or uncertainty-weighted. Details: `src/inhibitome/phenotypes/`.

### Aim 2 — Inhibitory fingerprints
Compute `Z_i` from synapses onto each neuron:
- **amount** (counts, fraction, per-unit-dendrite normalization, multisynaptic connections);
- **source** (input from proximal-, distal-, interneuron-targeting classes; broad *and* fine labels kept
  separate — a low-confidence 20-class label is not automatically better than a reliable 4-class one);
- **location** (soma / proximal / distal-basal / apical, validated first against the manual column cells
  and published compartment rules);
- **diversity** (entropy, dominant-class fraction, effective # of classes, source concentration,
  compartment entropy);
- **motifs** (shared/convergent/disinhibitory structure) — **secondary; excluded from the pilot.**

Raw counts alone are never used un-normalized (bigger dendrites collect more synapses). Details:
`src/inhibitome/fingerprints/`.

### Aim 3 — What the inhibitory connectome *adds*
Nested, out-of-sample model comparison `M0 -> M5` (technical -> cellular -> functional -> total-input ->
inhibitory-compartment -> full-fingerprint), each evaluated with grouped **leave-one-scan-out**
validation and benchmarked against **matched shuffled-fingerprint nulls**. Details:
`docs/02_PREREGISTRATION.md`, `src/inhibitome/models/`, `src/inhibitome/nulls/`.

---

## 8. Modeling discipline

- Start simple: hierarchical / regularized linear regression, then GAMs.
- Gradient-boosted trees only as a **nonlinear control** — never the primary evidence.
- **Never** start with GNNs / Transformers / deep nets: a flexible model can exploit imaging-depth and
  reconstruction artifacts and destroy the biological interpretation.
- Group all observations from one EM root ID into the same split (one neuron can map to several
  functional ROIs).
- The strongest statistical unit is the **scan / local cortical block**, not the individual synapse.

---

## 9. What success looks like (summary — full text in `docs/04`)

- **Minimum:** a rigorously validated inhibitory-fingerprint dataset + at least one endpoint (reliability
  *or* state modulation) that beats matched nulls out-of-sample across held-out scans.
- **Strong (ISEF-level):** source & placement of inhibition predict *distinct* computations —
  perisomatic → stable reliability, dendritic → state-dependent gain.
- **Exceptional:** static inhibitory wiring separates neurons into identifiable control regimes across
  layers and areas.

## 10. The limitation that governs everything

All matched EM+function data come from **one mouse**. This is also a stated limitation of the flagship
MICrONS analyses. No result here is a universal law of cortex; every result is a *principle observed
within the MICrONS specimen*, observational and un-perturbed. Unmeasured confounds we cannot rule out
include neuromodulator receptor expression, long-range inputs outside the EM volume, interneuron firing,
synaptic strength, intrinsic membrane properties, and incomplete axonal reconstruction.

See `docs/02_PREREGISTRATION.md` next.
