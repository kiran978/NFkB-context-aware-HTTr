#!/usr/bin/env python3
"""Reproduce paired HTTr ablation statistics from archived OOF predictions.
The supplied Phase-7 OOF file and scaffold mapping are the frozen predictions
used in the paper. Bootstrap resampling is performed at Murcko-scaffold level.
"""
import argparse, numpy as np, pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

def boot_delta(y,p0,p1,groups,n=2000,seed=42):
    rng=np.random.default_rng(seed); ug=np.unique(groups); vals=[]
    by={g:np.flatnonzero(groups==g) for g in ug}
    for _ in range(n):
        gs=rng.choice(ug,size=len(ug),replace=True)
        ix=np.concatenate([by[g] for g in gs])
        yy=y[ix]
        if len(np.unique(yy))<2: continue
        vals.append(average_precision_score(yy,p1[ix])-average_precision_score(yy,p0[ix]))
    vals=np.array(vals)
    return np.mean(vals),np.quantile(vals,.025),np.quantile(vals,.975),np.mean(vals>0),len(vals)

def main(a):
    oof=pd.read_csv(a.oof); nf=pd.read_csv(a.nfkb)[['InChIKey','Murcko_Scaffold']]
    oof=oof.merge(nf,on='InChIKey',how='left')
    results=[]
    for task in sorted(oof.Task.unique()):
        d=oof[oof.Task==task].dropna(subset=['Label','Structure_only_Probability','Structure_plus_HTTr_kNN_Probability','Murcko_Scaffold'])
        est,lo,hi,p,n=boot_delta(d.Label.to_numpy(int),d.Structure_only_Probability.to_numpy(float),d.Structure_plus_HTTr_kNN_Probability.to_numpy(float),d.Murcko_Scaffold.astype(str).to_numpy(),a.replicates,a.seed)
        results.append([task,est,lo,hi,p,n])
    pd.DataFrame(results,columns=['Task','Delta_AUPRC','CI95_Lower','CI95_Upper','P_Delta_gt0','Replicates']).to_csv(a.out,index=False)
if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--oof',required=True); p.add_argument('--nfkb',required=True); p.add_argument('--out',required=True); p.add_argument('--replicates',type=int,default=2000); p.add_argument('--seed',type=int,default=42); main(p.parse_args())
