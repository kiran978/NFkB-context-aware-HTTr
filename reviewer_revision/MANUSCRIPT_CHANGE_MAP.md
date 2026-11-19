# JCIM revision: manuscript change map for all remaining reviewer comments
**Working draft; freeze original model and 12 leads.** Do not submit
until new nested and EPA-control outputs have been independently checked.

## Major methods changes

1. **§3.5 — Model benchmarking and validation.**
   Cite Supplementary Table S10 (15 original candidate rows including the
   masked shared-trunk multi-task MLP), and the existing gated-ablation
   results (previous proposed Table S8). Disclose that the model family
   originally used for deployment was selected using scaffold-aware
   development comparisons, and that these comparisons alone were
   not independent selection-adjusted performance estimates.
   When the NEW nested analysis is run, specify common whole-matrix
   outer scaffold folds; inner folds; candidate list; calibration and
   HTTr decisions; held-out test metrics. Explicitly distinguish the
   historical method from the post-review nested sensitivity.

2. **§3.6 — Calibration.**
   Describe separately raw XGB (Phase7 ablation comparator AUPRC .514893)
   and Platt-calibrated deployed XGB (Phase6/Table2 AUPRC .505588).
   Both values can be valid under different model configurations;
   their AUPRC/Brier performance should not be mixed in one row.
   State exact cross-fitting procedure only if original code or new
   post-review code verifies it.

3. **§3.4 & §3.8 — Chemical split claims.**
   Replace "leakage-safe" with "Murcko-scaffold-grouped" or
   "scaffold-aware" in Figure 1 and text. Explain that separate
   Murcko scaffold strings do NOT prohibit high ECFP4 Tanimoto
   cross-fold analogues. Include source OOF max-similarity
   distributions (draft Table S14) and NEW ECFP4 connected-component
   blocked sensitivity (when actually run).
   Confirm reviewer-cited CID pairs have differing Murcko scaffolds.

4. **§3.11 — HTTr features: biological and mathematical meaning.**
   State the nine selected features by exact name and cell line:
   NFkB, TNFa, RELA, NFKB1, NFKB2,
   HALLMARK_TNFA_SIGNALING_VIA_NFKB,
   GILMORE_CORE_NFKB_PATHWAY,
   BIOCARTA_NFKB_PATHWAY and BIOCARTA_RELA_PATHWAY,
   each with _max suffix. State that the nine were a compact,
   biology-driven PRIOR rather than an exhaustive screen of all TLR
   and NFkB signatures. Do not say nine are optimal or unique.
   Related HTTr methods describe _max as median of two highest
   magnitude responses, _dose_at_max as the normalized dose location
   of the maximum, and _AUC_neg as negative-profile trapezoidal AUC
   (doi:10.1186/s13321-026-01231-4); verify source-version
   equivalence before describing the exact current dataset.
   Justify focusing on _max as a dimension-controlled peak-response
   design choice, not as experimentally proven superior to _AUC_neg.

5. **§3.12 — HTTr attribution controls.**
   The original 27 features are smoothed STRUCTURE-DEPENDENT
   estimates of EPA transcriptomic summaries, not measured B-cell
   transcriptomes. Test:
   (a) RDKit2D-only,
   (b) RDKit2D plus label-independent EPA ECFP4 neighborhood summaries,
   (c) RDKit2D plus original 27 NFkB responses,
   (d) same kNN averaging after within-cell-line EPA RESPONSE-VECTOR
       permutation across chemicals,
   (e) nine non-NFkB pathway response summaries per cell line with
       same kNN/weighting and predictor dimensionality.
   ALL model and feature-arm choices must occur inside outer training.
   Do not claim NFkB biological specificity unless (c) outperforms
   (b), (d) and (e) on held-out outer folds.

6. **§3.14 / §4.5 — Threshold sensitivity.**
   Add new retrospective threshold sensitivity (Table S12).
   Each original structure AD threshold multiplied .9,1,1.1;
   each entropy cutoff shifted -.02,0,+.02; HTTr all-cell AD and
   original probabilities fixed. Exact original counts verified
   199 strict, 193 strict-active, 172 tie-robust.
   New strict populations range 101..456. Do not say funnel or
   all four B-cell leads are insensitive: tightened +10% AD
   retains only 1/4 strict B-cell leads under the original
   entropy cutoff. No prospectively selected compounds are changed.

7. **§3.15–3.16 / §4.6 — Explain final 40-to-12.**
   Original Phase8D worksheet specifies 41 assessed including RA
   comparator; original 40 ranked candidate panel was 20 primary
   and 20 reserve. Seven primary automatically advance; three
   primary are retained to preserve rare phenotypes (P440AM34830,
   P440AM61633, P440AM872); reserve P440AM9226/P440AM61488
   are promoted. Score formula:
   0.50*BiologyRobustness + 0.35*Developability
   + 0.15*CLUE_NFkB_StructuralPrior_Percentile
   (all 0..100). Hard risk: simultaneous Lipinski+Veber
   failure OR ESOL logS<-5.5 OR >=2 NIH alerts OR aldehyde,
   epoxide, quinone or thiol alerts.
   Cite Table S13 and the underlying 40-candidate decision ledger.
   Emphasize score is a ranking aid; phenotype preservation is
   an explicit non-score-based exception. Structural prior is
   not target engagement.

## Existing-result corrections (no new model fitting)

8. **Table 3 and text, §4.4.** Attach full historical
   AUROC/AUPRC/MCC/balanced accuracy/Brier rows (Table S11).
   Show both observed AUPRC differences
   +.028274/-.054633/+.000459 and bootstrap-mean deltas
   +.027642/-.052279/-.001338. Their differences arise
   because bootstrapping a nonlinear metric does not, in
   general, reproduce the direct plug-in difference.
   Do not call the small AID2337 Δ a sign/data error; interval
   includes zero and there is no reliable benefit.

9. **Cytotoxicity, §3.2/4.1/Discussion.**
   AID489021 has only three positives among 469. It is an
   observed retrospective artifact annotation, not a trained
   predictive head or a prospective safety exclusion. This
   addresses Reviewer 1 Comment 1 and Reviewer 2 Minor 3 once.

10. **Docking, §3.17 and Supplementary methods.**
    Report Schrödinger release, LigPrep/Epik and protein prep
    settings, Glide SP/XP precision, receptor grid centers
    and box dimensions, ligand protonation/stereochemistry,
    retained poses and native ligand redocking metrics.
    Original 8V4X MALT1 RMSD=1.9 A and PKCβ 2I0E=2.4 A;
    the latter is NOT strict pose recovery.
    Source registry identifies three mutations in 2I0E
    construct; describe whether retained as deposited or
    reverted only after inspecting the actual job records.
    DO NOT invent software settings from current computer install.

11. **CyberAIDD/library, §3.13.**
    Cite actual CyberAIDD software/paper, specific version
    and generation rules from the original design record.
    Explain rosmarinic acid was used as a chemically motivated
    seed/reference comparator, not a validated parent
    inhibitor; preserve actual screening of 2000 IDs.
    Do not confuse the original library-generation tool
    with RDKit used for later structure QC.

12. **Reproducibility + Figure 1 + journal naming.**
    Figure 1 panel4 should describe prospective screening;
    correct panel1 typo to AID2337 / 283 compounds / 157
    actives. Change "leakage-safe" caption to scaffold-aware.
    Public GitHub main currently lacks historical model
    training code, the B-cell HTTr-augmented fitted model
    and claimed prospectives, and LICENSE contains author
    placeholder. Fix deposited artifact and licensing
    FIRST, or narrow claims in Data availability.
    README references old Journal of Cheminformatics and
    must match actual target journal, Journal of
    Chemical Information and Modeling (JCIM).
    The new reviewer branch does not replace missing
    historical model provenance.

## Conclusion claim to withhold until controls

Replace "a statistically supported NFkB-specific HTTr gain"
with "a preliminary positive association under the original
scaffold-aware OOF comparison, with the independent
selection-adjusted magnitude and source of incremental
information evaluated in additional post-review controls."

If the new controls show no biological specificity, reframe
as chemical-neighborhood-dependent transfer/auxiliary
response information and adjust title, abstract and conclusion
accordingly. Do not hide adverse control results.

## Open numerical slots (never invent)

- NEW nested per-task AUPRC/AUROC/MCC/Brier, fold-level
  selected model/HTTr retention and uncertainty.
- NEW matched chemical-neighborhood, EPA response-shuffle,
  non-NFkB response control metrics and paired CIs.
- NEW similarity-connected-component-blocked performance.
- Original Glide job settings, CyberAIDD version,
  author-confirmed license and final public archive links.
