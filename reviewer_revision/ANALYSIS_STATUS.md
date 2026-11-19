# NF-kB manuscript reviewer-response validation — working audit (19 Nov 2026)

> **Status:** NEW post-review analysis. This branch is NOT the historic model-development
> record and must not be merged into a claimed original-code archive without explicit
> provenance labels and execution verification.
> Original Phase-6/7 probabilities, original Phase-8 2,000-compound screen, and the
> final 12 candidates MUST remain unchanged.

## Verified source inputs

The review audit used these files from the author's research archive:
- NFkB_Phase5_Multitask_Matrix_751.csv (751 compounds / 497 scaffold groups)
- NFkB_Phase6_Deployment_OOF.csv (751 compound rows; 3 tasks)
- NFkB_Phase7D_Transferred_HTTr_Biology_751.csv (751 compounds / 27 features)
- NFkB_Phase7D_HTTr_kNN_OOF.csv (1,221 labeled prediction records)
- NFkB_Phase7D_HTTr_kNN_Teacher_Ablation.csv
- NFkB_Phase7D_HTTr_kNN_ScaffoldBootstrap_Delta_CI.csv
- NFkB_Phase7D_HTTr_kNN_Reliability.csv (EPA nearest-neighbor similarity controls)
- NFkB_Phase8B_All_2000_Predictions.csv
- NFkB_Phase8C_Prioritized_All_2000.csv
- NFkB_Phase8D_Final_12_Leads.csv
- NFkB_Phase8D_Orthogonal_Validation.xlsx
- NFkB_Phase5_Model_Benchmarks.csv

Data files are not automatically available in this GitHub branch. Run scripts locally
using the original archived datasets and preserve their input checksums.

## Completed: direct reconstruction of original Phase-7 AUPRC

| Task | Structure AP | Structure + NF-kB HTTr AP | Direct AP difference | Reported bootstrap-mean difference |
|---|---:|---:|---:|---:|
| Bcell489023 | 0.589164 | 0.617438 | +0.028274 | +0.027642 |
| TNFa489022 | 0.514893 | 0.460260 | -0.054633 | -0.052279 |
| TNFa2337 | 0.764125 | 0.764584 | +0.000459 | -0.001338 |

The bootstrap mean is **not** the same estimand as the arithmetic
difference between AUPRCs displayed in Table 3. Explain both quantities
explicitly. Table 2's TNFa489022 AUPRC 0.505588 refers to the separately
Platt-calibrated deployment model; Table 3's 0.514893 is the raw XGB comparator.

## Completed: original Phase-8 screening integrity and gate sensitivity

The source file contains 2,000 unique IDs. Calculated original counts exactly
equal archived flags: 199 strict-confidence, 193 strict+any-HC-active, and
172 tie-robust strict; also 20 primary, 20 reserve and final 12 IDs.

NEW retrospective sensitivity: multiply each frozen structure AD threshold by
0.9, 1.0 or 1.1, and adjust each entropy cutoff by -0.02, 0 or +0.02.
All frozen probabilities, the EPA all-cell HTTr AD gate and tie-sensitive
annotations were fixed. These are threshold perturbations, NOT newly
estimated 1st/10th percentile training thresholds.

| AD factor | Entropy shift | All three structural AD | Strict | Strict + active | Tie-robust strict | Original 12 still strict | Four strict B-cell selected leads still strict |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.9 | -0.02 | 643 | 218 | 215 | 182 | 5 | 1 |
| 0.9 | 0.00 | 643 | 309 | 300 | 255 | 9 | 4 |
| 0.9 | +0.02 | 643 | 456 | 437 | 373 | 11 | 4 |
| 1.0 | -0.02 | 394 | 154 | 151 | 135 | 5 | 1 |
| 1.0 | 0.00 | 394 | 199 | 193 | 172 | 9 | 4 |
| 1.0 | +0.02 | 394 | 282 | 270 | 234 | 11 | 4 |
| 1.1 | -0.02 | 241 | 101 | 98 | 88 | 4 | 1 |
| 1.1 | 0.00 | 241 | 126 | 122 | 110 | 5 | 1 |
| 1.1 | +0.02 | 241 | 169 | 162 | 144 | 7 | 1 |

**Interpretation:** Original exact counts are reproducible; substantial count
sensitivity remains. Do NOT claim threshold invariance. Candidate panel was
not re-ranked or changed. Strict B-cell lead eligibility in particular
is sensitive to tightening the AD threshold.

## Completed: scaffold versus nearest-neighbor audit (archived OOF values)

The reviewer-cited CID pairs are both in the 751-compound matrix and have
different Murcko scaffold strings:
- 2211281 vs 2279440
- 24856284 vs 24978754

Archived Phase-6 OOF maximum ECFP4 similarity columns:
| Task | N | Median max similarity | OOF rows with max similarity >=0.85 | >=0.90 |
|---|---:|---:|---:|---:|
| Bcell489023 | 469 | 0.396552 | 9 | 2 |
| TNFa489022 | 469 | 0.400000 | 7 | 2 |
| TNFa2337 | 283 | 0.593220 | 4 | 4 |

These numbers are derived from the recorded OOF similarity columns; the
original source script that computed them remains unverified. An additional
connected-component ECFP4 similarity-blocked evaluation is supplied as a
NEW sensitivity script, not as proof the original CV was similarity-blocked.

## Pre-existing 40-to-12 lead decision ledger

The original NFkB_Phase8D_Orthogonal_Validation.xlsx documents
7 automatic primary advances + 3 phenotype-preserving primary additions
+ 2 promotions from reserve (P440AM9226 and P440AM61488).
Readiness score = 0.50 biological robustness + 0.35 developability
+ 0.15 CLUE structural-prior percentile (components on 0..100 scale).

Hard-risk rule: both Lipinski and Veber fail; OR ESOL logS < -5.5;
OR >=2 NIH alerts; OR aldehyde/epoxide/quinone/thiol motif.
Structural alerts are triage warnings, not measured toxicology.
CLUE nearest-neighbor target annotations are hypotheses, not target validation.

## Code in this review branch

- run_prospective_sensitivity.py: repeats the Phase-8B nine-grid threshold
  check with strict baseline count assertions; uses original data only.
- run_nested_validation.py: NEW, selection-aware 5 outer x 3 inner grouped
  CV for three tasks. Compares RDKit2D-only, RDKit2D + EPA-neighborhood
  descriptors, RDKit2D + 27 NF-kB HTTr variables. Independently selects
  family, small parameter set, raw/Platt and feature arm inside outer train.
  Supports --groups ecfp4_085 similarity-connected components; --controls
  adds EPA response permutation and non-NFkB transfer arms once available.
- build_httr_controls.py: reconstructs EPA response shuffle BEFORE
  Tanimoto kNN averaging and a matched nine-feature non-NFkB response
  control from the **original broad EPA response matrix**.

### Inputs for scripts

Put the named source files in a local folder (not necessarily GitHub),
then run e.g.:

python reviewer_revision/run_prospective_sensitivity.py --phase8b "data/NFkB_Phase8B_All_2000_Predictions.csv" --phase8c "data/NFkB_Phase8C_Prioritized_All_2000.csv" --final12 "data/NFkB_Phase8D_Final_12_Leads.csv" --out "reviewer_outputs/gate_sensitivity.csv"

python reviewer_revision/run_nested_validation.py --matrix "data/NFkB_Phase5_Multitask_Matrix_751.csv" --httr "data/NFkB_Phase7D_Transferred_HTTr_Biology_751.csv" --reliability "data/NFkB_Phase7D_HTTr_kNN_Reliability.csv" --out "reviewer_outputs/scaffold"

python reviewer_revision/run_nested_validation.py --matrix "data/NFkB_Phase5_Multitask_Matrix_751.csv" --httr "data/NFkB_Phase7D_Transferred_HTTr_Biology_751.csv" --reliability "data/NFkB_Phase7D_HTTr_kNN_Reliability.csv" --groups ecfp4_085 --out "reviewer_outputs/similarity085"

To use **two crucial EPA controls**, supply the original broad 4,142-row
chemical/cell HTTr response source, its DTXSID-to-SMILES mapping, and a
JSON list of nine *exact* non-NFkB response columns. Then generate
control features with build_httr_controls.py and pass --controls to
run_nested_validation.py. Never substitute a shuffle of the already
transferred 751-row features: that is not the requested EPA-response control.

## NOT YET COMPLETED / no fabricated results

- New nested outer-CV training and ECFP4-connected blocked CV have not run:
  the current model-side Python execution service failed. No new nested
  AP values, calibration decisions or model-selection frequencies can be claimed.
- The parent broad EPA response matrix and nine exact non-NFkB summaries
  were not supplied with current extracted files; EPA shuffle and
  non-NFkB controls therefore remain unexecuted.
- Docking Glide release, grid coordinates/dimensions, receptor/ligand
  prep and pose count require original Schrödinger output/job metadata.
- CyberAIDD version/generation settings require original design logs.
- The public default branch still lacks the historical executable model
  and fitted HTTr B-cell deployment model. This draft review branch
  contains NEW scripts, **not a substitute** for original provenance.
- The author licence placeholder and the wrong journal name in original
  README remain to be corrected by authenticated submission and licence
  ownership details before a release is finalized.

**Do not merge/submit these scripts as "original source code".
Review locally, run in the intended environment, inspect outputs,
then document changes and whether analysis was prospective or post hoc.**


## Recovered local original EPA source path (from research workspace inventory)

The author's previous local inventory identifies the broad source table:
E:\College\Research\Pharmacodia Libararies\Chemoinformatics\1\Git Files\chemical_httr_assay_aggregated.feather

It also identifies the original NF-kB-only subset:
E:\College\Research\Pharmacodia Libararies\Chemoinformatics\1\Additional_file_4_Reproducibility_archive\data\NFkB_Phase7C_HTTr_Selected_Features.csv

The second is not sufficient for Reviewer 2's non-NFkB comparator.
Use extract_epa_response_controls.py --feather "<original broad table>"
--output-dir "<local output>" to inspect its true column names without
uploading the ~394 MB source file. Supply a nine-name JSON list selected
on biological rationale BEFORE looking at assay-specific performance.
Then export and retransfer from the small subset.

### Interpreting EPA-derived suffixes

A related 2026 HTTr/Tox21 methods paper (Cattebeke et al.,
doi:10.1186/s13321-026-01231-4) describes _max as the median of the two
highest-magnitude concentration responses; _dose_at_max as the dose level
of the greatest absolute response; and _AUC_neg as trapezoidal AUC for
negative responses (sign convention inverted for some interpretations).
Confirm equivalence with this project's actual source feature generation.
The nine targeted NF-kB _max summaries were a compact biology-driven
selection; do NOT assert a measured advantage over _AUC_neg without
the corresponding sensitivity analysis and actual source values.
