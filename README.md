# Reproducibility archive — context-specific NF-kB cheminformatics study

This archive contains the source code, curated machine-readable inputs, frozen out-of-fold predictions, final fitted model, HTTr transfer inputs, prospective predictions and summary results required to reproduce the quantitative conclusions reported in the Journal of Cheminformatics submission.

## Directory map
- `source_code/03_phase5_context_aware_nfkB.py` — original Phase-5 context-aware data/model analysis script.
- `source_code/04_phase6_context_model.py` — original Phase-6 gated-model, calibration, uncertainty and applicability-domain script.
- `source_code/05_reproduce_httr_transfer.py` — reconstructs the 27 transferred HTTr features from the compact EPA feature table supplied here.
- `source_code/06_reproduce_httr_statistics.py` — paired scaffold-bootstrap comparison of structure-only vs HTTr-augmented OOF predictions.
- `source_code/07_reproduce_prospective_summary.py` — reproduces the prospective-screen funnel counts.
- `data/NFkB_Phase5_Multitask_Matrix_751.csv` — curated 751-compound multi-context matrix.
- `data/NFkB_Phase7C_HTTr_Selected_Features.csv` — compact 4,142-row EPA HTTr feature subset used for transfer/feasibility analyses.
- `data/SMILES_data.csv`, `data/chemical_dtxsid_map.csv` — EPA chemical identity/structure mappings used for the HTTr source set.
- `results/` — frozen Phase-6/7/8 OOF predictions, confidence intervals and prospective outputs.
- `models/NFkB_Phase6_Final_Context_Model.joblib` — frozen fitted structure-only context model bundle used before the Phase-7 B-cell augmentation layer was added.

## Public source data
The underlying HTTr experiments are public through GEO: GSE272548 (MCF-7), GSE274318 (U-2 OS) and GSE284321 (HepaRG). The EPA high-throughput transcriptomics pipeline is openly available from USEPA/CompTox-httrpl. PubChem assay data are public under AID 489023, 489022, 489021, 2337 and 1183838.

The compact `NFkB_Phase7C_HTTr_Selected_Features.csv` supplied in this archive contains the exact selected HTTr summaries used in this study, so reproduction of the reported transfer and statistical comparisons does not require downloading the full ~394 MB project intermediate Feather matrix.

## Environment
The analysis environment used for final audit is listed in `requirements.txt`. Random seeds and fold logic are documented in the source scripts; the main scaffold-CV seed was 42.

## Reproduce core numerical checks
1. Run `04_phase6_context_model.py` against the curated matrix using its documented input paths, or inspect the frozen OOF and bootstrap outputs in `results/`.
2. Recreate the 27 transfer features with `05_reproduce_httr_transfer.py` if desired.
3. Run `06_reproduce_httr_statistics.py` to reproduce the paired scaffold-bootstrap HTTr delta analysis.
4. Run `07_reproduce_prospective_summary.py` on the frozen Phase-8 files to reproduce the reported screening funnel.

## Important interpretation boundary
The study deliberately rejected a universal transcriptomic teacher. HTTr augmentation was retained only for Bcell489023 after paired scaffold-bootstrap evidence; the TNF models remained structure-only. The target/docking analyses are hypothesis-generating and are not required to reproduce the primary cheminformatics conclusions.

## License
Source code in this supplementary archive is intended to be released under the MIT License. Replace the `[AUTHORS TO INSERT]` placeholder in `LICENSE` before public deposition/submission.
