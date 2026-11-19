#!/usr/bin/env python3
"""Reviewer-response analysis: independent outer-fold evaluation of all selections.

NEW ANALYSIS (2026-11): not the original Phase-6/7 training script.
Never overwrite archived OOF probabilities, fitted models, or prospective outputs.

Runs matched structure-only, chemical-neighborhood, and NF-kB HTTr arms,
plus optional EPA response-permuted and non-NF-kB HTTr controls.
Each outer-training partition independently chooses:
- model family / small hyperparameter candidate set via inner scaffold-aware OOF AP;
- raw vs cross-fitted Platt probabilities via inner-training Brier score;
- among the 3 feature arms, the arm with best inner OOF AP.
Reports held-out OUTER-fold metrics for each arm and for the complete selector.
This measures selection-aware performance; no outer-test labels are used to choose.

Usage:
 python reviewer_revision/run_nested_validation.py --matrix data/NFkB_Phase5_Multitask_Matrix_751.csv \
   --httr data/NFkB_Phase7D_Transferred_HTTr_Biology_751.csv \
   --reliability data/NFkB_Phase7D_HTTr_kNN_Reliability.csv --out reviewer_outputs/scaffold
 python reviewer_revision/run_nested_validation.py --matrix ... --httr ... \
   --reliability ... --groups ecfp4_085 --out reviewer_outputs/similarity085

EPA response permutations and NON-NF-kB biological controls require ORIGINAL EPA
chemical-response tables and are included only when --controls is supplied.
"""
import argparse, json, hashlib, sys
from pathlib import Path
from collections import Counter
import numpy as np, pandas as pd
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn import __version__ as SK_VERSION
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import average_precision_score, roc_auc_score, brier_score_loss, matthews_corrcoef
try:
    from xgboost import XGBClassifier, __version__ as XGB_VERSION
except ImportError:
    XGBClassifier = None
    XGB_VERSION = "not installed"

TASKS = [
    ("Bcell489023", "TaskMask_Bcell489023", "Bcell489023_Label"),
    ("TNFa489022", "TaskMask_TNFa489022", "TNFa489022_Label"),
    ("TNFa2337", "TaskMask_TNFa2337", "TNFa2337_Label"),
]
CANDIDATES = ["LR_C1", "SVM_C1", "SVM_C2", "RF_leaf2"] + (
    ["XGB_depth3"] if XGBClassifier is not None else []
)

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def descriptor_matrix(smiles):
    funs = [f for _, f in Descriptors.descList]
    out = np.empty((len(smiles), len(funs)), dtype=float)
    mols = []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            raise ValueError("Unparseable SMILES at row %d" % i)
        mols.append(mol)
        for j, fn in enumerate(funs):
            try:
                out[i, j] = float(fn(mol))
            except Exception:
                out[i, j] = np.nan
    out[~np.isfinite(out)] = np.nan
    return out, mols

def read_reliability(path, compound_ids):
    if path is None:
        return None
    d = pd.read_csv(path)
    required = {"CID", "Cell_Line", "Max_Tanimoto", "Mean_Top5_Tanimoto"}
    if not required.issubset(d):
        raise ValueError("EPA neighborhood file missing columns: %s" % (required-set(d)))
    p = d.pivot_table(index="CID", columns="Cell_Line",
                      values=["Max_Tanimoto", "Mean_Top5_Tanimoto"],
                      aggfunc="mean")
    p.columns = ["EPA_" + str(metric) + "_" + str(cell) for metric, cell in p.columns]
    p = p.reindex(pd.Index(compound_ids.astype(int)))
    return p.to_numpy(dtype=float)

def ecfp4_components(mols, cutoff=0.85):
    """Connected components: any pair with similarity >= cutoff stays together."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(m) for m in mols]
    parent = list(range(len(fps)))
    def root(i):
        while i != parent[i]:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        for j, s in enumerate(sims):
            if s >= cutoff:
                a, b = root(i), root(j)
                if a != b:
                    parent[b] = a
    return np.array(["COMP_%d" % root(i) for i in range(len(fps))])

def make_estimator(name, y_train, seed):
    prep = [("imp", SimpleImputer(strategy="median")),
            ("var", VarianceThreshold(0))]
    if name.startswith("LR"):
        model = LogisticRegression(C=1.0, class_weight="balanced", max_iter=2000, random_state=seed)
        return Pipeline(prep + [("scale", StandardScaler()), ("model", model)])
    if name.startswith("SVM"):
        c = 1.0 if name == "SVM_C1" else 2.0
        model = SVC(C=c, gamma=0.005, kernel="rbf", class_weight="balanced",
                    probability=True, random_state=seed)
        return Pipeline(prep + [("scale", StandardScaler()), ("model", model)])
    if name == "RF_leaf2":
        model = RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                       class_weight="balanced_subsample",
                                       n_jobs=1, random_state=seed)
        return Pipeline(prep + [("model", model)])
    if name == "XGB_depth3":
        y_train = np.asarray(y_train, dtype=int)
        npos = max(int(y_train.sum()), 1)
        spw = float((len(y_train) - y_train.sum()) / npos)
        model = XGBClassifier(n_estimators=250, max_depth=3, learning_rate=0.05,
                              min_child_weight=2, subsample=0.9, colsample_bytree=0.8,
                              reg_lambda=3, scale_pos_weight=spw,
                              tree_method="hist", eval_metric="logloss",
                              n_jobs=1, random_state=seed)
        return Pipeline(prep + [("model", model)])
    raise ValueError(name)

def crossfit_probs(x, y, gidx, folds, model_name, seed):
    """OOF probabilities within ONE outer-training partition, no outer-test rows."""
    p = np.full(len(gidx), np.nan, dtype=float)
    local_lookup = {int(g): i for i, g in enumerate(gidx)}
    for train_global, valid_global in folds:
        tr = [local_lookup[int(t)] for t in train_global if int(t) in local_lookup]
        va = [local_lookup[int(t)] for t in valid_global if int(t) in local_lookup]
        if not tr or not va or len(np.unique(y[tr])) < 2:
            raise ValueError("Empty or one-class inner task fold; revise grouping deterministically")
        est = make_estimator(model_name, y[tr], seed)
        est.fit(x[tr], y[tr])
        p[va] = est.predict_proba(x[va])[:, 1]
    if not np.isfinite(p).all():
        raise ValueError("Incomplete inner OOF predictions")
    return p

def platt_inputs(p):
    p = np.clip(np.asarray(p, dtype=float), 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p)).reshape(-1, 1)

def choose_calibration(p_inner, y_inner, idx_global, inner_folds, seed):
    """Calibration comparison uses only INNER data; outer holdout remains untouched."""
    y_inner = np.asarray(y_inner, dtype=int)
    lookup = {int(g): i for i, g in enumerate(idx_global)}
    p_cross = np.full(len(p_inner), np.nan)
    for train_global, valid_global in inner_folds:
        tr = [lookup[int(t)] for t in train_global if int(t) in lookup]
        va = [lookup[int(t)] for t in valid_global if int(t) in lookup]
        if len(np.unique(y_inner[tr])) < 2:
            raise ValueError("One-class Platt fit partition")
        cal = LogisticRegression(C=1.0, max_iter=2000, random_state=seed)
        cal.fit(platt_inputs(p_inner[tr]), y_inner[tr])
        p_cross[va] = cal.predict_proba(platt_inputs(p_inner[va]))[:, 1]
    raw_brier = brier_score_loss(y_inner, p_inner)
    calibrated_brier = brier_score_loss(y_inner, p_cross)
    choice = "platt" if calibrated_brier < raw_brier - 0.0001 else "raw"
    cal_final = LogisticRegression(C=1.0, max_iter=2000, random_state=seed)
    cal_final.fit(platt_inputs(p_inner), y_inner)
    return choice, cal_final, raw_brier, calibrated_brier

def metrics(df):
    if not len(df):
        return {}
    y = df["Label"].to_numpy(dtype=int)
    p = df["Probability"].to_numpy(dtype=float)
    return dict(N=len(y), Active=int(y.sum()),
                AUPRC=average_precision_score(y, p),
                AUROC=roc_auc_score(y, p) if len(np.unique(y))==2 else np.nan,
                MCC=matthews_corrcoef(y, (p >= 0.5).astype(int)),
                Brier=brier_score_loss(y, p))

def run(a):
    matrix = pd.read_csv(a.matrix).reset_index(drop=True)
    required = ["InChIKey", "CID", "Standardized_SMILES", "Murcko_Scaffold"]
    for x in required:
        if x not in matrix:
            raise ValueError("Matrix missing column " + x)
    if matrix.InChIKey.duplicated().any():
        raise ValueError("Input matrix must be unique by standardized InChIKey")
    httr = pd.read_csv(a.httr)
    cols = [c for c in httr.columns if c.endswith("_kNN")]
    if len(cols) != 27:
        raise ValueError("Expected exactly 27 transferred HTTr columns; found %d" % len(cols))
    merged = matrix[["InChIKey"]].merge(httr[["InChIKey"] + cols],
                                        on="InChIKey", how="left", validate="one_to_one")
    if len(merged) != len(matrix):
        raise ValueError("Missing matched HTTr rows")
    x_struct, mols = descriptor_matrix(matrix.Standardized_SMILES)
    x_httr = merged[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    x_neigh = read_reliability(a.reliability, matrix.CID)
    arms = {"Structure": x_struct, "NFkB_HTTr": np.column_stack((x_struct, x_httr))}
    if x_neigh is not None:
        arms["EPA_Neighborhood"] = np.column_stack((x_struct, x_neigh))
    if a.controls:
        ctrl = pd.read_csv(a.controls)
        shuffled_cols = sorted(c for c in ctrl if c.startswith("SHUF_") and c.endswith("_kNN"))
        other_cols = sorted(c for c in ctrl if c.startswith("NONNF_") and c.endswith("_kNN"))
        if len(shuffled_cols) != 27 or len(other_cols) != 27:
            raise ValueError("Control file must contain 27 SHUF and 27 NONNF columns")
        aligned = matrix[["InChIKey"]].merge(ctrl[["InChIKey"] + shuffled_cols + other_cols],
                                             on="InChIKey", how="left", validate="one_to_one")
        xs = aligned[shuffled_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        xn = aligned[other_cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        arms["EPA_ResponsePermuted"] = np.column_stack((x_struct, xs))
        arms["EPA_NonNFkB"] = np.column_stack((x_struct, xn))
    if a.groups == "scaffold":
        groups = matrix.Murcko_Scaffold.fillna("").astype(str).to_numpy()
        groups = np.where(groups == "", matrix.InChIKey.astype(str), groups)
    else:
        groups = ecfp4_components(mols, cutoff=0.85)
    bmask = pd.to_numeric(matrix.TaskMask_Bcell489023, errors="coerce").fillna(0).eq(1).to_numpy()
    yb = pd.to_numeric(matrix.Bcell489023_Label, errors="coerce").fillna(0).astype(int).astype(str)
    yt = pd.to_numeric(matrix.TNFa2337_Label, errors="coerce").fillna(0).astype(int).astype(str)
    global_strata = np.where(bmask, "B_" + yb, "T_" + yt)
    outer_cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=a.seed)
    outer = list(outer_cv.split(np.zeros(len(matrix)), global_strata, groups))
    all_predictions, fold_audit = [], []
    for fold_idx, (global_train, global_test) in enumerate(outer, start=1):
        inner_cv = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=a.seed + fold_idx)
        inner = [(global_train[tr], global_train[val])
                 for tr, val in inner_cv.split(np.zeros(len(global_train)),
                                                np.asarray(global_strata)[global_train],
                                                groups[global_train])]
        for task, maskcol, labelcol in TASKS:
            m = pd.to_numeric(matrix[maskcol], errors="coerce").fillna(0).eq(1).to_numpy()
            global_t = global_train[m[global_train]]
            global_v = global_test[m[global_test]]
            if len(global_v) == 0:
                raise ValueError("Missing task test fold: %s fold %d" % (task, fold_idx))
            y_t = pd.to_numeric(matrix[labelcol], errors="coerce").iloc[global_t].to_numpy(int)
            y_v = pd.to_numeric(matrix[labelcol], errors="coerce").iloc[global_v].to_numpy(int)
            if len(np.unique(y_t)) < 2 or len(np.unique(y_v)) < 2:
                raise ValueError("One-class task fold: %s fold %d" % (task, fold_idx))
            choices = {}
            for arm_name, xall in arms.items():
                x_t = xall[global_t]
                best_name, best_ap, best_p = None, -1.0, None
                for model_name in CANDIDATES:
                    p = crossfit_probs(x_t, y_t, global_t, inner, model_name, a.seed)
                    ap = average_precision_score(y_t, p)
                    if ap > best_ap:
                        best_name, best_ap, best_p = model_name, float(ap), p
                decision, calibrator, rb, cb = choose_calibration(
                    best_p, y_t, global_t, inner, a.seed)
                final_model = make_estimator(best_name, y_t, a.seed)
                final_model.fit(x_t, y_t)
                p_raw = final_model.predict_proba(xall[global_v])[:, 1]
                p_test = calibrator.predict_proba(platt_inputs(p_raw))[:, 1] \
                    if decision == "platt" else p_raw
                choices[arm_name] = (best_ap, best_name, decision, p_test)
                for ind, prob in zip(global_v, p_test):
                    all_predictions.append(dict(Task=task, OuterFold=fold_idx,
                        Mode=a.groups, Arm=arm_name, CID=int(matrix.CID.iloc[ind]),
                        InChIKey=matrix.InChIKey.iloc[ind],
                        Murcko_Scaffold=matrix.Murcko_Scaffold.iloc[ind],
                        Label=int(pd.to_numeric(matrix[labelcol].iloc[ind])),
                        Probability=float(prob)))
                fold_audit.append(dict(Task=task, OuterFold=fold_idx, Mode=a.groups,
                    Arm=arm_name, SelectedModel=best_name,
                    InnerAUPRC=best_ap, Calibration=decision,
                    InnerBrierRaw=rb, InnerBrierPlattCrossfit=cb,
                    OuterTrainN=len(global_t), OuterTestN=len(global_v)))
            # The WHOLE arm-selection procedure is evaluated solely on outer-test.
            winning_arm = max(choices, key=lambda n: choices[n][0])
            for ind, prob in zip(global_v, choices[winning_arm][3]):
                all_predictions.append(dict(Task=task, OuterFold=fold_idx,
                    Mode=a.groups, Arm="ADAPTIVE_SELECTOR",
                    CID=int(matrix.CID.iloc[ind]),
                    InChIKey=matrix.InChIKey.iloc[ind],
                    Murcko_Scaffold=matrix.Murcko_Scaffold.iloc[ind],
                    Label=int(pd.to_numeric(matrix[labelcol].iloc[ind])),
                    Probability=float(prob)))
            fold_audit.append(dict(Task=task, OuterFold=fold_idx, Mode=a.groups,
                Arm="ADAPTIVE_SELECTOR", SelectedModel=choices[winning_arm][1],
                InnerAUPRC=choices[winning_arm][0],
                Calibration=choices[winning_arm][2], SelectedFeatureArm=winning_arm,
                OuterTrainN=len(global_t), OuterTestN=len(global_v)))
        print("Completed outer fold %d / 5" % fold_idx, flush=True)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    pred = pd.DataFrame(all_predictions)
    audit = pd.DataFrame(fold_audit)
    pred.to_csv(out/"outer_predictions.csv", index=False)
    audit.to_csv(out/"fold_selection_audit.csv", index=False)
    met = []
    for (task, arm), g in pred.groupby(["Task", "Arm"]):
        met.append(dict(Task=task, Arm=arm, Mode=a.groups, **metrics(g)))
    pd.DataFrame(met).to_csv(out/"nested_performance.csv", index=False)
    manifest = dict(analysis="NEW_NESTED_SENSITIVITY_NOT_ORIGINAL_PIPELINE",
       grouping=a.groups, outer_folds=5, inner_folds=3, seed=a.seed,
       candidates=CANDIDATES, arms=list(arms),
       matrix_sha256=sha256(a.matrix), httr_sha256=sha256(a.httr),
       neighborhood_sha256=sha256(a.reliability) if a.reliability else None,
       rdkit_version=rdBase.rdkitVersion, sklearn_version=SK_VERSION,
       xgboost_version=XGB_VERSION,
       controls_sha256=sha256(a.controls) if a.controls else None,
       missing_controls=[] if a.controls else ["EPA_RESPONSE_PERMUTED_BEFORE_NEIGHBOR_AVERAGING",
                         "NON_NFKB_EPA_RESPONSE_SUMMARIES"])
    (out/"analysis_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n",encoding="utf-8")
    print(pd.DataFrame(met).to_string(index=False))

def cli():
    p=argparse.ArgumentParser()
    p.add_argument("--matrix",required=True)
    p.add_argument("--httr",required=True)
    p.add_argument("--reliability")
    p.add_argument("--controls", help="Optional EPA permuted/non-NFkB control features CSV")
    p.add_argument("--out",required=True)
    p.add_argument("--seed",type=int,default=42)
    p.add_argument("--groups",choices=["scaffold","ecfp4_085"],default="scaffold")
    return p.parse_args()

if __name__ == "__main__":
    run(cli())
