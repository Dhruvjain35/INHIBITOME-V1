# 02 — Pre-registration (freeze before the final experiment)

This document is **frozen before** the Day 8–9 blinded comparison. Its purpose is to prevent fishing:
the endpoints, the model hierarchy, and the null models are all fixed in advance. Anything discovered
after unblinding is labeled **exploratory** and can never rescue a failed confirmatory endpoint.

Freeze protocol: tag the commit that finalizes this file `prereg-frozen`, record the git SHA in
`outputs/prereg_lock.json` (written by `scripts/05`), and do not edit the hierarchy or nulls afterward.

---

## 1. Confirmatory endpoints

**Primary endpoint.** Out-of-sample improvement in predicting **state-dependent visual gain** `M_i`
when adding the inhibitory fingerprint to the cellular + functional baseline:
```
dR2_state = R2_oos(M2 + inhibitory fingerprint) - R2_oos(M2)
```

**Secondary endpoint.** Out-of-sample improvement in predicting **repeated-stimulus response
reliability** `R_i`:
```
dR2_reliability = R2_oos(M2' + inhibitory fingerprint) - R2_oos(M2')
```
where `M2'` is the functional baseline **with reliability itself removed** from the predictors.

**Exploratory (labeled, never confirmatory):** stimulus-specific modulation; pupil-vs-locomotion
dissociation; disinhibitory motifs; visual-area interactions; nonlinear structural interactions; fine
inhibitory connectivity groups; branch-level inhibition.

**Decision rule.** If state modulation fails but the pre-registered reliability endpoint succeeds
strongly (beats matched nulls out-of-sample across scans), the project continues with a **narrower
claim**. If both fail, stop — do not add targets until something is significant.

---

## 2. The nested model hierarchy (M0 -> M5)

Each model adds a block of predictors. The target `Y_i` is either `M_i` (primary) or `R_i` (secondary).
Every model is fit and evaluated identically (§4). We report `R2_oos` for each, and the *increments*.

| Model | Adds | Purpose |
|---|---|---|
| **M0** technical baseline | scan/session, cortical depth, imaging quality, EM boundary distance | soak up pure acquisition/reconstruction artifacts |
| **M1** cellular baseline | + area, layer, morphological type, dendritic size | control for cell identity |
| **M2** functional baseline | + baseline activity, visual tuning, response amplitude, reliability* | control for what the neuron already does |
| **M3** total-input control | + total excitatory input, total inhibitory input | control for connectivity *amount* |
| **M4** inhibitory compartment | + perisomatic / proximal / distal / apical inhibitory fractions | *does placement matter?* |
| **M5** full fingerprint | + presynaptic inhibitory classes, source diversity, connection redundancy | *does source identity matter?* |

\* For the **reliability** endpoint, reliability is obviously excluded from M2 (call it `M2'`).

**Primary tests:**
- Does **M5** generalize better than **M2** and **M3**? (the fingerprint adds signal beyond function and beyond total connectivity)
- Does **M4** beat **M3**? (compartment placement matters beyond amount)
- Does **M5** beat **M4**? (source identity matters beyond placement)

A model that improves *training* fit but not *held-out* `R2_oos` is **not** evidence.

---

## 3. Required adversarial null models

Each null re-runs the full M5 evaluation with one aspect destroyed. The fingerprint is credible only if
it beats **every** applicable null.

| # | Null | Preserves | Destroys | Question it answers |
|---|---|---|---|---|
| N1 | **Excitatory-input control** | build the *same* fingerprint from excitatory inputs | inhibition-specificity | is the signal specific to inhibition, or just "connectivity"? |
| N2 | **Layer/morphology-matched permutation** | area, layer, depth, morph class, dendritic size | which neuron got which fingerprint | is it more than cell-type-goes-with-cell-type? |
| N3 | **Compartment-location null** | per-neuron inhibitory synapse count | which compartment each synapse hits | does *location* matter beyond amount? |
| N4 | **Presynaptic-identity null** | counts + compartments; shuffle source class among matched inputs | which interneuron class provides input | does *source identity* matter? |
| N5 | **Degree-preserving graph null** | in/out degrees; rewire locally | specific wiring | is it the wiring, or just degree? |
| N6 | **Imaging-depth control** | match/regress out depth + signal quality | depth-correlated artifacts | is it biology or optics? |
| N7 | **Behavioral time-shift null** | traces, shifted by large offsets | true temporal alignment of L/P to activity | are `M_i` estimates real, not chance structure? |

Nulls N2–N5 are permutation tests: recompute the observed statistic under ≥1,000 matched shuffles, report
the empirical p-value and where the observed increment sits in the null distribution.

---

## 4. Validation & inference (identical for every model and null)

**Grouping.** All rows sharing an EM `root_id` stay in the same fold (one neuron → possibly many
functional ROIs). No leakage across the anatomical unit.

**Primary split — leave-one-scan-out.** Train on all-but-one imaging scan, test on the held-out scan;
rotate over all scans. Random neuron splits are **not** the primary evidence (they leak scan-specific
imaging/behavior). We *also* report random-grouped splits only as a secondary diagnostic.

**Stress splits (secondary).** Hold out entire layers / morphological classes / cortical areas / local
spatial blocks. Generalizing to a morphological class unseen in training is the strongest possible
result.

**Inference under one-mouse dependence.** Neurons are not independent biological replicates. Use:
scan-clustered bootstrap; spatial-block bootstrap; hierarchical random effects (scan as grouping);
permutation tests preserving area/layer/morphology; FDR correction across structural features;
confidence intervals on *every* increment. **Never** report a tiny p-value from millions of synapses
while ignoring the single-animal design. Report effect per scan, per area/layer, variation across folds,
all exclusions, and all null results.

---

## 5. Model classes (in order; stop escalating once signal is clear)

1. hierarchical linear regression (scan random effect)
2. regularized regression (ridge/elastic-net) with nested CV for the penalty
3. generalized additive models (nonlinearity, still interpretable)
4. gradient-boosted trees — **nonlinear control only**, never the headline

Explicitly **not** used as primary evidence: GNNs, Transformers, deep MLPs.

---

## 6. Frozen thresholds (proposed gates, conservative)

- Data gate (Day 3): ≥1,000 high-confidence coregistered neurons; ≥10 functional scans; hundreds of
  neurons with reliable inhibitory fingerprints; usable behavior across multiple scans.
- Phenotype gate (Day 5): `M_i` split-half test–retest `Corr ≥ 0.3` for the neurons retained (others
  down-weighted or dropped); oracle reliability reproduced within tolerance of the released score.
- Signal gate (Day 10): `dR2` positive and CI-excludes-zero in a **majority of scan folds**, and the
  observed increment beats matched nulls N2 + N3 (minimum) at empirical p < 0.05 after FDR.

These are proposed thresholds, not guaranteed counts. If the data undershoot, that is itself a finding —
report it (see `docs/04`).
