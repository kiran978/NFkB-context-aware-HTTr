#!/usr/bin/env python3
"""Paired Murcko-scaffold bootstrap of archived Phase-7 HTTr OOF predictions.

Uses the published Phase-7 OOF columns CID, Task, Label,
Structure_Probability and StructurePlusHTTr_Probability. Scaffolds are
joined from the curated 751-compound matrix by CID. This is an independent
reconstruction utility, not the original fitted model or its exact RNG trace.

Example:
  python source_code/06_reproduce_httr_statistics.py \
    --oof results/NFkB_Phase7D_HTTr_kNN_OOF.csv \
    --nfkb data/NFkB_Phase5_Multitask_Matrix_751.csv \
    --out reproduced_phase7_delta.csv
"""
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

def boot_delta(y, p0, p1, groups, n=2000, seed=42):
    rng = np.random.default_rng(seed)
    unique_groups = np.unique(groups)
    indices = {g: np.flatnonzero(groups == g) for g in unique_groups}
    deltas = []
    for _ in range(n):
        sampled = rng.choice(unique_groups, size=len(unique_groups), replace=True)
        ix = np.concatenate([indices[g] for g in sampled])
        yy = y[ix]
        if len(np.unique(yy)) != 2:
            continue
        deltas.append(average_precision_score(yy, p1[ix]) - average_precision_score(yy, p0[ix]))
    if not deltas:
        raise ValueError('No valid bootstrap replicates')
    values = np.asarray(deltas)
    return float(values.mean()), float(np.quantile(values, .025)), float(np.quantile(values, .975)), float((values > 0).mean()), int(len(values))

def main(a):
    oof = pd.read_csv(a.oof)
    matrix = pd.read_csv(a.nfkb)
    required = {'CID', 'Task', 'Label', 'Structure_Probability', 'StructurePlusHTTr_Probability'}
    if missing := required.difference(oof.columns):
        raise ValueError(f'Phase-7 OOF columns missing: {sorted(missing)}')
    mapping = matrix[['CID', 'Murcko_Scaffold']].drop_duplicates('CID')
    if len(mapping) != len(matrix):
        raise ValueError('Duplicate CIDs in curated matrix; cannot uniquely map scaffolds')
    oof = oof.merge(mapping, on='CID', how='left', validate='many_to_one')
    if oof.Murcko_Scaffold.isna().any():
        raise ValueError('One or more OOF CIDs lack a Murcko scaffold in the curated matrix')
    rows = []
    for task in sorted(oof.Task.unique()):
        d = oof[oof.Task == task].dropna(subset=['Label', 'Structure_Probability', 'StructurePlusHTTr_Probability', 'Murcko_Scaffold'])
        y = d.Label.to_numpy(dtype=int)
        p0 = d.Structure_Probability.to_numpy(dtype=float)
        p1 = d.StructurePlusHTTr_Probability.to_numpy(dtype=float)
        delta = average_precision_score(y, p1) - average_precision_score(y, p0)
        mean, lo, hi, p, n = boot_delta(y, p0, p1, d.Murcko_Scaffold.astype(str).to_numpy(), a.replicates, a.seed)
        rows.append([task, len(d), delta, mean, lo, hi, p, n])
    pd.DataFrame(rows, columns=['Task', 'N', 'Observed_Delta_AUPRC', 'Bootstrap_Mean_Delta_AUPRC', 'CI95_Lower', 'CI95_Upper', 'P_Delta_gt0', 'Replicates']).to_csv(a.out, index=False)

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--oof', required=True)
    p.add_argument('--nfkb', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--replicates', type=int, default=2000)
    p.add_argument('--seed', type=int, default=42)
    main(p.parse_args())
