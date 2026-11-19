# Draft supplementary tables for reviewer revision
**Working draft — 19 Nov 2026. Not part of original deposited supplementary material.**

All numerical entries below were read from the original archived Phase-5/7/8
tables or recalculated from archived, frozen predictions. None of the values
in this draft is a claim about the UNRUN nested model-selection evaluation.

## Table S10. Phase-5 model-family benchmarks (original archived results)

| Task | Architecture/model | AUROC | AUPRC | MCC | Brier |
|---|---|---:|---:|---:|---:|
| B-cell AID489023 | Single-task LR | .655159 | .544171 | .203219 | .272037 |
| B-cell AID489023 | Single-task RBF-SVM | .695454 | .589164 | .286241 | .216242 |
| B-cell AID489023 | Single-task RF | .676287 | .561234 | .202147 | .219703 |
| B-cell AID489023 | Single-task XGB | .645029 | .527112 | .174199 | .251575 |
| B-cell AID489023 | Masked shared-trunk MT-MLP | .664497 | .546794 | .289600 | .301095 |
| TNFα AID489022 | Single-task LR | .681990 | .448866 | .262077 | .252346 |
| TNFα AID489022 | Single-task RBF-SVM | .703039 | .484670 | .174018 | .194487 |
| TNFα AID489022 | Single-task RF | .699061 | .506122 | .203048 | .192403 |
| TNFα AID489022 | Single-task XGB | .711817 | .512191 | .209530 | .205558 |
| TNFα AID489022 | Masked shared-trunk MT-MLP | .699693 | .458614 | .237105 | .262110 |
| TNFα AID2337 | Single-task LR | .662547 | .693439 | .273228 | .287146 |
| TNFα AID2337 | Single-task RBF-SVM | .764356 | .764125 | .474428 | .192691 |
| TNFα AID2337 | Single-task RF | .733900 | .742339 | .401348 | .209954 |
| TNFα AID2337 | Single-task XGB | .730361 | .749968 | .353521 | .221973 |
| TNFα AID2337 | Masked shared-trunk MT-MLP | .699348 | .736600 | .279309 | .272548 |

*Source:* original NFkB_Phase5_Model_Benchmarks.csv. These are development
benchmark values; the winning model's performance is not proven
selection-adjusted. The masked multi-task MLP was evaluated and should be
mentioned in Methods and table captions.

## Table S11. Complete historical Phase-7 HTTr ablation metrics

| Task | Input | AUROC | AUPRC | MCC | Balanced accuracy | Brier |
|---|---|---:|---:|---:|---:|---:|
| B-cell AID489023 | Structure-only | .695454 | .589164 | .286241 | .632286 | .216242 |
| B-cell AID489023 | Structure + NFkB HTTr | .711432 | .617438 | .288417 | .629881 | .211108 |
| TNFα AID489022 | Structure-only **raw XGB** | .711480 | .514893 | .228802 | .612644 | .204427 |
| TNFα AID489022 | Structure + NFkB HTTr | .668856 | .460260 | .174329 | .580923 | .217117 |
| TNFα AID2337 | Structure-only | .764356 | .764125 | .474428 | .733040 | .192691 |
| TNFα AID2337 | Structure + NFkB HTTr | .759251 | .764584 | .489142 | .738626 | .195267 |

*Source:* original NFkB_Phase7D_HTTr_kNN_Teacher_Ablation.csv;
Phase-7 1,221-row OOF predictions were independently checked for AUPRC.
These compare already-chosen model families; not selection-adjusted.

**Correct Table 3 footnotes:** The *observed* paired AUPRC differences are
+0.028274 (B-cell), -0.054633 (TNFα AID489022), +0.000459 (TNFα AID2337).
The reported +0.027642, -0.052279, -0.001338 are respective **means of
2,000 scaffold bootstrap AUPRC differences** and do not always equal
direct plug-in differences, nor even their sign for near-zero effects.
The interval for TNFα2337 includes zero. Final AID489022 deployment
AUPRC .505588 in Table 2 uses a **different calibrated predictor** than
the .514893 raw comparator here. Do not conflate them.

## Table S12. Frozen prospective sensitivity — new retrospective analysis

Each model's *fixed* structural AD cut-off is multiplied by the indicated
factor; each *fixed* entropy cut-off is shifted additively (capped at 1).
EPA all-cell HTTr AD and tie flags are held constant. No rescoring.

| AD factor | Entropy shift | All-three AD | Strict confidence | Strict + at least one active | Tie-robust strict | Final12 in strict | Strict B-cell four in strict |
|---|---:|---:|---:|---:|---:|---:|---:|
| .90 | -.02 | 643 | 218 | 215 | 182 | 5 | 1 |
| .90 | .00 | 643 | 309 | 300 | 255 | 9 | 4 |
| .90 | +.02 | 643 | 456 | 437 | 373 | 11 | 4 |
| 1.00 | -.02 | 394 | 154 | 151 | 135 | 5 | 1 |
| 1.00 | .00 | 394 | 199 | 193 | 172 | 9 | 4 |
| 1.00 | +.02 | 394 | 282 | 270 | 234 | 11 | 4 |
| 1.10 | -.02 | 241 | 101 | 98 | 88 | 4 | 1 |
| 1.10 | .00 | 241 | 126 | 122 | 110 | 5 | 1 |
| 1.10 | +.02 | 241 | 169 | 162 | 144 | 7 | 1 |

*Source:* recalculation using frozen NFkB_Phase8B_All_2000_Predictions.csv
with the original Phase-8C selection records and Phase-8D final 12 IDs.
Original row-level flags reproduced EXACTLY (199,193,172).
Sensitivity scenarios alter *retrospective eligibility*, not the frozen
prospective list or historical 12-lead selection. Counts are sensitive to
modest changes; biological and mechanistic conclusions remain provisional.

## Table S13. Phase-8D original 40-to-12 decision ledger

Original multiobjective scoring and orthogonal-screening outputs are in
NFkB_Phase8D_Orthogonal_Validation.xlsx.

| Final rank | ID | Original selection pathway | Readiness score (0-100) |
|---|---|---|---:|
| 1 | P440AM19590 | Primary automatic advance | 85.52 |
| 2 | P440AM34458 | Primary automatic advance | 84.32 |
| 3 | P440AM31921 | Primary automatic advance | 81.44 |
| 4 | P440AM10109 | Primary automatic advance | 80.59 |
| 5 | P440AM9226 | Reserve promotion | 75.29 |
| 6 | P440AM34824 | Primary automatic advance | 75.28 |
| 7 | P440AM61488 | Reserve promotion | 73.47 |
| 8 | P440AM32581 | Primary automatic advance | 67.08 |
| 9 | P440AM33770 | Primary automatic advance | 61.01 |
| 10 | P440AM34830 | Primary phenotype preservation | 56.44 |
| 11 | P440AM61633 | Primary phenotype preservation | 54.51 |
| 12 | P440AM872 | Primary phenotype preservation | 54.20 |

Readiness = .50 (biology robustness) + .35 (developability) +
.15 (CLUE nearest reference percentile), with 0–100 component scales.
Hard risk: both Lipinski and Veber rules fail OR ESOL logS < -5.5
OR >=2 NIH alerts OR aldehyde/epoxide/quinone/thiol motif.
Three compounds were designated for phenotype preservation rather than
strictly descending readiness rank; the 40-to-12 decision is not a
single scalar-threshold operation. CLUE similarity is a *structural
hypothesis prior*, NOT validated target engagement.

## Table S14. Archived Phase-6 fold-neighbor ECFP4 similarity audit

| Assay | Annotated OOF compounds | Median recorded max-neighbor ECFP4 Tanimoto | Number >=.85 | Number >=.90 |
|---|---:|---:|---:|---:|
| B-cell AID489023 | 469 | .396552 | 9 | 2 |
| TNFα AID489022 | 469 | .400000 | 7 | 2 |
| TNFα AID2337 | 283 | .593220 | 4 | 4 |

Both reviewer-specified pairs (CID 2211281 vs 2279440, and CID
24856284 vs 24978754) are in the curated 751-compound matrix and
have *different* Murcko scaffold strings. Thus, scaffold grouping alone
does not block all near-duplicate fingerprint similarity. These counts
come from the archived OOF nearest-similarity fields and cannot replace
the proposed NEW similarity-connected-blocked CV performance test.

## Status and boundaries

Tables S10/S11/S13/S14 are recovered original-development evidence
or direct description of recorded OOF fields. Table S12 is a NEW
post-review, deterministic sensitivity computation. None substitutes
for selection-adjusted nested metrics, EPA response-shuffle controls or
non-NFkB response controls. These remain separate planned computations
in reviewer_revision/run_nested_validation.py and build_httr_controls.py.
