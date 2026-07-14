# 05 — ISEF Positioning & Publication

## ISEF category

**Computational Biology and Bioinformatics**, sub-area **Computational Neuroscience** — the work uses
computational and mathematical methods to study information processing by neural structures.

## The story to tell (judges, not reviewers)

Do **not** lead with dataset size or ML architecture. Lead with the biology:

1. Every cortical neuron receives a *different pattern* of inhibition — different amounts, from different
   interneurons, onto different parts of the cell.
2. MICrONS is the one place we can *see* those individual inhibitory synapses **and** watch the very same
   neuron compute during vision and changing behavioral state.
3. We test whether that static inhibitory wiring predicts two things a neuron does: respond *reliably*
   vs. be *modulated* by arousal/locomotion.
4. We use **matched rewiring controls** — shuffle where the synapses land, or which interneuron they came
   from — to prove it's the *location* and *identity* that matter, not just "more synapses."
5. We say exactly **where static anatomy succeeds and where it fails.**

That narrative is mechanistic, honest about limits, and visually powerful.

## Figures that carry the poster

- **F1 — The idea.** One excitatory neuron, its inhibitory synapses colored by compartment
  (soma/proximal/distal/apical), beside its calcium trace across rest vs. running. (`skeleton_plot` +
  activity.)
- **F2 — Two phenotypes.** Distributions of reliability `R_i` and state modulation `M_i` across the
  cohort; show they're distinct.
- **F3 — The headline.** Nested-model `dR2` ladder (M0→M5) for each endpoint, with leave-one-scan-out
  error bars, and the matched-null distributions overlaid. This is the whole argument in one panel.
- **F4 — Dissociation (if it holds).** Perisomatic inhibition → reliability; dendritic inhibition →
  state gain. A 2×2 that makes the mechanistic claim land.
- **F5 — Honesty panel.** Where it fails: per-scan spread, the one-mouse caveat, effects that vanish
  under a null. Judges reward this.

## Judging-defense one-liners (rehearse these)

- *"Isn't this just structure-relates-to-function, already known?"* → No — we test the *single-neuron
  inhibitory-input fingerprint* (source **and** compartment), not connectivity in general, and we beat
  matched nulls that hold cell type fixed.
- *"It's one mouse."* → Correct, and we say so everywhere. Every claim is a principle *within the MICrONS
  specimen*; it's the only animal on Earth with matched EM + function at this scale. We control for scan,
  depth, and optics precisely because we can't average across animals.
- *"Did the model just overfit?"* → Held-out leave-one-scan-out `dR2`, not training fit; interpretable
  linear/GAM models, not a black box; gradient-boosting only as a control.

## Publication

Realistic targets (specialist venues): *eNeuro*, *Network Neuroscience*, *Journal of Neuroscience*,
*PLOS Computational Biology*, *Cerebral Cortex*, or a neuroinformatics venue.

Why **not** a top general journal, stated plainly:
- only one mouse has matched EM + function;
- the study is observational, no experimental perturbation;
- external replication of the exact structural relationship is unavailable.

A good specialist publication is realistic **if** the effect is strong, the controls are adversarial, and
the full analysis is released reproducibly (this repo, pinned materialization, seeded splits, all nulls,
all null results reported).

## Reproducibility bar (do this regardless of outcome)

- Pin materialization 1300; record provenance JSON.
- Seed every split/bootstrap (`config → random_seed`).
- Commit every `outputs/*.md` report, including negative results.
- Tag `prereg-frozen` before unblinding.
- A one-command rerun: `make pilot`.
