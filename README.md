# Reliability-aware context-specific NF-κB chemoinformatics — reproducibility repository

Supporting materials for the manuscript **“Reliability-aware context-specific chemoinformatics with selective transcriptomic transfer for prospective NF-κB inhibitor prioritization”** (Journal of Molecular Graphics and Modelling submission; see article for current venue and status).

## Repository inventory (verified)

### Analysis / reconstruction utilities
- `source_code/05_reproduce_httr_transfer.py` — archived helper to reconstruct 27 nearest-neighbour HTTr-derived features from the source EPA feature subset; this is a *reconstruction utility*, not the fitted deployment model.
- `source_code/06_reproduce_httr_statistics.py` — archived paired scaffold-bootstrap HTTr ablation utility using saved out-of-fold (OOF) predictions.
- `source_code/07_reproduce_prospective_summary.py` — archived Phase-8 prospective-screen funnel summary utility using saved prospective CSVs.

### Results available in this repository
- `results/NFkB_Phase5_Model_Benchmarks.csv` — complete 15-row Phase-5 benchmark across model families and assay contexts.
- `results/NFkB_Phase8D_Final_12_Leads.csv` — final 12-compound prospective lead audit with predictions and selection metadata.

Additional saved evaluation tables:
- `results/NFkB_Phase6_Deployment_Metrics.csv` — Phase-6 locked structure-model metrics.
- `results/NFkB_Phase6_Deployment_ScaffoldBootstrap_CI.csv` — scaffold-bootstrap confidence limits.
- `results/NFkB_Phase7D_HTTr_kNN_Teacher_Ablation.csv` — matched Phase-7 HTTr ablation metrics.
- `results/NFkB_Phase7D_HTTr_kNN_ScaffoldBootstrap_Delta_CI.csv` — Phase-7 paired scaffold-bootstrap differences.

### Other repository files
- `requirements.txt` — package versions for the archived final-audit environment, not necessarily all historical modeling runs.
- `DATA_SOURCES.md` — public assay and transcriptomic source locations.

## Supplementary article archive

The study's `Additional_file_4_Reproducibility_archive` package includes the full curated assay and EPA input files, Phase-6 and Phase-7 saved OOF predictions, Phase-8 all-2,000-compound results and a fitted **Phase-6 structure-only** context-model bundle. The Phase-6 fitted bundle should not be represented as the deployed Phase-7 HTTr-augmented B-cell model.

The original archive also contains documentation-only outlines named `03_phase5_context_aware_nfkB.py` and `04_phase6_context_model.py`; these are **not executable model-training scripts**. They are not listed here as reproducible source code. The three published utilities above reproduce selected processing/statistics from saved data but do not independently establish an exact end-to-end rerun of the original fitting and screening stages.

## Scientific interpretation

Transferred HTTr features were an exploratory auxiliary representation rather than directly measured responses in the assay compounds. Later similarity-blocked and pathway-control analyses qualified the original development-stage B-cell improvement. See the revised manuscript and corresponding supplementary sensitivity analyses for the full interpretation.

## Data provenance and licensing

PubChem assay identifiers: 489023, 489022, 489021, 2337 and 1183838. EPA HTTr source records: GEO GSE272548 (MCF-7), GSE274318 (U-2 OS) and GSE284321 (HepaRG). Data reuse remains subject to the respective upstream source terms.

A repository `LICENSE` file is present. Its copyright-holder line requires confirmation from the contributing rights holders before it can be considered final; that legal attribution should not be guessed from authorship alone.
