#!/usr/bin/env python3
"""Reconstruct the 27 NF-kB-focused HTTr transfer features from the compact
EPA HTTr feature subset supplied in data/NFkB_Phase7C_HTTr_Selected_Features.csv.

The source table contains the 4,142 chemical-cell-context rows and selected
response summaries. Structures are mapped by DTXSID using data/SMILES_data.csv.
For every NF-kB query molecule, ECFP4/Tanimoto similarity to EPA chemicals is
computed separately by cell line. The five most similar EPA structures are
combined with similarity^2 weights. Nine response features are transferred per
cell line, giving 27 auxiliary features.

Usage:
  python source_code/05_reproduce_httr_transfer.py \
    --nfkb data/NFkB_Phase5_Multitask_Matrix_751.csv \
    --httr data/NFkB_Phase7C_HTTr_Selected_Features.csv \
    --smiles data/SMILES_data.csv --out reproduced_httr_751.csv

This is a reconstruction utility, not the original frozen deployment model.
"""
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

FEATURES=[
'NFkB_max','TNFa_max','RELA_max','NFKB1_max','NFKB2_max',
'HALLMARK_TNFA_SIGNALING_VIA_NFKB_max','GILMORE_CORE_NFKB_PATHWAY_max',
'BIOCARTA_NFKB_PATHWAY_max','BIOCARTA_RELA_PATHWAY_max']
CELL_LINES=['MCF7','U2OS','HepRG']

def fp(smiles):
    m=Chem.MolFromSmiles(str(smiles))
    return None if m is None else AllChem.GetMorganGenerator(radius=2,fpSize=2048).GetFingerprint(m)

def weighted_knn(qfp, ref_fps, vals, k=5):
    sims=np.array([DataStructs.TanimotoSimilarity(qfp,x) for x in ref_fps])
    idx=np.argsort(-sims)[:k]
    w=sims[idx]**2
    if w.sum()==0: return np.nan, float(sims.max())
    return float(np.average(vals[idx],weights=w)), float(sims.max())

def main(a):
    nf=pd.read_csv(a.nfkb)
    ht=pd.read_csv(a.httr)
    sm=pd.read_csv(a.smiles).rename(columns={'dsstox_substance_id':'outcome_id'})[['outcome_id','SMILES']]
    ht=ht.merge(sm,on='outcome_id',how='left')
    out=nf[['CID','InChIKey','Standardized_SMILES']].copy()
    for cell in CELL_LINES:
        sub=ht[ht.cell_type.astype(str).str.upper()==cell.upper()].copy()
        agg=sub.groupby(['outcome_id','SMILES'],dropna=False)[FEATURES].mean().reset_index()
        agg=agg[agg.SMILES.notna()].copy()
        agg['_fp']=[fp(s) for s in agg.SMILES]
        agg=agg[agg['_fp'].notna()].copy()
        ref_fps=list(agg['_fp'])
        for feat in FEATURES:
            vals=agg[feat].astype(float).to_numpy()
            transferred=[]; maxsim=[]
            for s in nf.Standardized_SMILES:
                q=fp(s)
                if q is None: transferred.append(np.nan); maxsim.append(np.nan); continue
                v,ms=weighted_knn(q,ref_fps,vals,5)
                transferred.append(v); maxsim.append(ms)
            out[f'{cell}_{feat}_kNN']=transferred
            if feat==FEATURES[0]: out[f'{cell}_MaxEPA_Tanimoto']=maxsim
    out.to_csv(a.out,index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--nfkb',required=True); p.add_argument('--httr',required=True); p.add_argument('--smiles',required=True); p.add_argument('--out',required=True); main(p.parse_args())
