# 01 — Data Access

Everything here was verified against the live MICrONS / CAVE docs (July 2026). Exact strings are
mirrored in `config/pilot.yaml`; **change them there, not in code.** Source URLs are at the bottom.

---

## 0. Two data planes (this trips people up)

MICrONS data lives in **two separate systems**, and this project needs both:

| Plane | What | How you get it |
|---|---|---|
| **CAVE** | EM anatomy: synapses, cell types, compartments, coregistration, proofreading | `caveclient` (Python), queried by materialization version |
| **DANDI (NWB)** | Function: calcium activity, treadmill/locomotion, pupil, eye position, repeated-movie stimuli | `dandi` download of dandiset **000402** |

The coregistration table (CAVE) is the bridge: it maps an EM `pt_root_id` to a functional
`(session, scan_idx, unit_id)`, which is how you find that neuron's activity in the NWB files.

---

## 1. CAVE setup (one time per machine)

Package: **`caveclient`** (v8.2.1 verified on PyPI, 2026-07-10). A token **is required even for public
data**.

```python
from caveclient import CAVEclient

# First time only — mint and save a token:
client = CAVEclient()                    # no datastack yet
client.auth.setup_token(make_new=True)   # opens a browser to log in via global.daf-apis.com
client.auth.save_token(token="PASTE_TOKEN_HERE")   # add overwrite=True to replace

# You must also accept the Terms of Service once (Google-associated email):
#   https://global.daf-apis.com/sticky_auth/api/v1/tos/2/accept
```

Token is stored at `~/.cloudvolume/secrets/cave-secret.json` (key: `token`). `scripts/00_setup_cave_token.py`
walks through this and verifies access. **Never commit the token** (`.gitignore` covers the secrets dir
and `.env`).

Then, for all real work:
```python
client = CAVEclient("minnie65_public")
client.version = 1300      # pin the materialization (config/pilot.yaml)
```

## 2. Materialization versions

- List: `client.materialize.get_versions()` · timestamp: `client.materialize.get_timestamp(v)`.
- **Live public versions (Mar 2026):** 117, 943, **1300**, 1507, 1621, 1718.
- **Expired (static download only, not queryable):** 343, 661, 795, 1078, 1181, 1412.
- Materializations are removed ~2 years after release, so **pin one now** and record it. We pin **1300**.
- ⚠️ Some table names contain `_v661` even though the v661 *materialization* is expired — the table name
  and the materialization version are independent. Query those tables at version 1300 fine.

## 3. Tables used (all in `config/pilot.yaml`)

**Connectivity**
- `synapses_pni_2` — ~337M synapses. Columns: `pre_pt_root_id`, `post_pt_root_id`, `*_pt_position`, `size`.
  Query onto a neuron: `client.materialize.synapse_query(post_ids=root_id)`. **500,000-row cap per
  query** — chunk by id batches for hub neurons.
- `synapse_target_predictions_ssa` — ~204M spine/shaft/soma compartment predictions (2025 release). This
  is what turns "N inhibitory synapses" into "N inhibitory synapses *on the soma / proximal / distal /
  apical* compartment."

**Cell identity**
- `baylor_log_reg_cell_type_coarse_v1` — excitatory vs inhibitory.
- `aibs_metamodel_celltypes_v661` — broadest-coverage soma/nucleus classifier (Elabbady et al. 2025).
- `aibs_metamodel_mtypes_v661_v2` — data-driven morphological types (census).
- `allen_column_mtypes_v1` — within-V1-column m-types.
- `allen_v1_column_types_slanted_ref` — **expert-curated column cells (~2,204)**; the manual validation
  set for compartment rules and inhibitory classes.

**Soma / area / QC**
- `nucleus_detection_v0` (~144K), `nucleus_ref_neuron_svm`, `nucleus_functional_area_assignment`
  (V1/AL/RL/LM), `proofreading_status_and_strategy`.

**Coregistration (the bridge)**
- `coregistration_manual_v4` — **PRIMARY cohort: 15,352 root IDs / 19,181 ROIs**, human-verified.
- `coregistration_auto_phase3_fwd_apl_vess_combined_v2` — automated, ~83K ROIs (secondary / expansion).
- Linking columns: nucleus id, `pt_root_id`, `session`, `scan_idx`, `unit_id`, plus quality metrics
  `residual` / `score`. Query: `client.materialize.query_table("coregistration_manual_v4")`.
- `digital_twin_properties_bcm_coreg_v4` — orientation/direction tuning etc. on the manual coreg cells
  (same schema), a convenient source of `visual tuning` baseline features.

## 4. Functional data (DANDI 000402, NWB)

- Dandiset **000402**, "MICrONS two-photon functional imaging." Download:
  ```bash
  pip install dandi
  dandi download DANDI:000402            # into data/raw/dandi_000402/
  ```
- Each NWB has calcium activity, behavior (treadmill rotation → locomotion velocity; eye video → pupil
  diameter + eye position), and stimulus timing. EM root IDs are attached as plane-segmentation columns
  derived from CAVE, so you can align to `pt_root_id`.
- Alternative source: the original DataJoint DB / SWDB Code Ocean xarray asset (heavier setup; NWB is
  the recommended path here).

## 5. The oracle (reliability) score

- Stimulus: six ~1-minute movies presented **10× per session** ("oracle" repeats).
- **Oracle score** = mean signal correlation between each presentation's response and the average of the
  other nine → a per-cell reliability metric (higher = more reliable). This is our `R_i` target.
- ⚠️ **[verify at ingest]** The definition and that it's a *derived property* are confirmed; a single
  CAVE column literally named `oracle_score` was **not** confirmed. Plan: **compute it ourselves from
  the NWB repeated-movie responses** (`src/inhibitome/phenotypes/reliability.py`), and if a released
  score is present, use it only as a cross-check. Never trust it blind.

## 6. Python libraries

| Package | Purpose | Version note |
|---|---|---|
| `caveclient` | CAVE query / materialization / auth | **8.2.1 verified** |
| `standard-transform` | voxel/nm → cortical-depth/layer-aligned coords | pin via `pip index versions` |
| `meshparty` | mesh↔skeleton, annotation handling | ~1.16.x (confirm) |
| `nglui` | build Neuroglancer states for figures/QC | confirm |
| `cloud-volume` | read meshes/skeletons/imagery | confirm |
| `pcg-skel` | skeletonize from the chunkedgraph | confirm |
| `skeleton_plot` | 2D compartment-colored skeleton figures | confirm |

Only `caveclient`'s version is firmly verified; pin the rest at `uv sync` time and record in `uv.lock`.

## 7. Key papers (cite these)

- **Flagship:** MICrONS Consortium, "Functional connectomics spanning multiple areas of mouse visual
  cortex," *Nature* 640:435–447 (2025). doi:10.1038/s41586-025-08790-w.
- **Inhibitory census:** Schneider-Mizell et al., "Inhibitory specificity from a connectomic census of
  mouse visual cortex," *Nature* 640:448–458 (2025). doi:10.1038/s41586-024-07780-8. (Correction: *Nature*
  642, E9, 2025.)
- **Wiring rule / digital twin:** "Functional connectomics reveals a general wiring rule in mouse visual
  cortex," *Nature* (2025). doi:10.1038/s41586-025-08840-3.
- **Cell types:** Elabbady et al., "Perisomatic ultrastructure efficiently classifies cells in mouse
  cortex," 2025 (basis for `aibs_metamodel_celltypes_v661`).

## 8. Flagged / unverified (confirm at ingest, do not hard-code trust)

- Exact token file path (`~/.cloudvolume/secrets/cave-secret.json` — very likely, confirm on first run).
- The literal `field`/`unit_id` addressing and any dedicated `oracle_score` column.
- Current PyPI versions of standard-transform / meshparty / nglui / cloud-volume / skeleton_plot.
- A public digital-twin *model-weights* download URL (not needed for the pilot — we use the released
  tuning-property tables, not the model itself).
- **A specific 2026 state-dependent functional-network preprint could NOT be verified.** Search bioRxiv
  directly before citing it in the novelty section; the plan does not depend on it existing.

## Sources

CAVEclient PyPI · MICrONS CAVEclient setup tutorial · CAVE authentication guide · materialization-version
tutorial · synapse-query tutorial · annotation-tables tutorial · SWDB functional-data book · DANDI 000402 ·
Nature flagship (PMC11981939) · inhibitory census · wiring-rule paper · awesome-connectomics.
(Full URLs in the research log; the tutorial root is https://tutorial.microns-explorer.org/.)
