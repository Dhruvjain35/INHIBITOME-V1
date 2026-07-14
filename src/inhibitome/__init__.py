"""INHIBITOME-V1 — inhibitory synaptic architecture vs. stable/state-dependent visual computation.

Package layout mirrors the three aims (see docs/00_PROJECT_PLAN.md):

    data/         Aim 0  — pin MICrONS materialization, join to one row per EM neuron.
    phenotypes/   Aim 1  — reliability R_i (oracle) and state modulation M_i.
    fingerprints/ Aim 2  — inhibitory input fingerprints Z_i (amount/source/location/diversity).
    models/       Aim 3  — nested model hierarchy M0..M5 + grouped leave-one-scan-out validation.
    nulls/                — the seven adversarial matched-null controls.
    report/               — sample accounting / QC reports.

Nothing here estimates across animals: this is one mouse. The statistical unit is the scan/block.
"""

__version__ = "0.1.0"
