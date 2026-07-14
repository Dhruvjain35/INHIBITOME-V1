# 04 — Kill & Success Criteria

The honest core of the project. Written **before** results so they can't be rationalized after.

---

## Success tiers

**Minimum successful research result.** A rigorously validated inhibitory-fingerprint dataset **plus**
at least one of:
- inhibitory architecture predicts **response reliability** beyond morphology and total input, or
- inhibitory architecture predicts **state modulation** beyond morphology, tuning, and reliability.

The effect must **reproduce across held-out scans** and **beat matched nulls**. This alone is a
defensible, releasable study.

**Strong (ISEF-level) result.** Cellular **source** and subcellular **placement** of inhibition predict
*distinct* computations: perisomatic inhibition → stable reliability; dendritic/apical inhibitory
organization → behavioral-state-dependent visual gain.

**Exceptional result.** Static inhibitory wiring separates neurons into structurally identifiable
**control regimes** (stable / amplified / suppressed / stimulus-dependent) across layers and areas.

---

## What does NOT count as success

These are too weak to support ISEF-winning claims and must be reported as such:
- inhibitory synapse *count* merely correlates with activity;
- a black-box model gets slightly higher `R2` (no interpretable, null-beating increment);
- one cortical layer shows an *uncorrected* association;
- effects vanish after controlling for morphology;
- random neuron splits work but **leave-one-scan-out fails**;
- only automated **low-confidence** labels produce the effect;
- the true fingerprint is **no better than matched shuffled fingerprints**.

---

## Kill criteria (stop the project)

Kill if **any** hold:
- fewer than **1,000** trustworthy neurons survive the join;
- compartment labels missing/unreliable for **most** inputs;
- state-modulation estimates are **not reproducible** (fail the split-half gate);
- inhibitory features only re-encode layer/morphology (fail N2);
- effects disappear under **leave-one-scan-out**;
- matched shuffled fingerprints perform **equally well**;
- results are driven by **one** scan/layer/area;
- **optical/depth quality** explains the findings (fail N6);
- automated annotations disagree badly with the **manual census**;
- only a **large neural network** produces an effect.

**Partial-kill / narrow rule.** If **state modulation fails** but the pre-registered **reliability**
endpoint succeeds strongly (null-beating, cross-scan), continue with the **narrower claim**. If **both**
fail, stop — do not bolt on new targets to chase significance.

---

## Why a clean kill is a good outcome

A negative result that is *rigorously established* — "static inhibitory anatomy does **not** add
predictive signal beyond morphology and tuning for these phenotypes in the MICrONS specimen" — is a real,
reportable, ISEF-defensible finding. It is far better than a fragile positive that dies under a reviewer's
first null model. The 10-day structure is designed so a kill costs days, not months.
