#!/usr/bin/env python3
"""Build two biologically interpretable HTTr negative controls from RAW EPA responses.

1. SHUF_*: shuffle complete 9-feature EPA response vectors among EPA chemicals
   WITHIN cell line before ECFP4 k=5 similarity^2 neighbor averaging.
2. NONNF_*: transfer nine separately specified non-NFkB EPA response summaries
   by the IDENTICAL chemical-neighbor procedure.

CANNOT be reconstructed from only the final 27-feature NFkB matrix.
Requires original broad EPA response summaries and original EPA SMILES mapping.
Choose and document nine exact non-NFkB columns in a JSON list in advance.
Non-NFkB response summaries are NOT necessarily biologically inert; this is a
specificity comparator, not a biological null.

Output matches 751 NF-kB compounds by InChIKey. Does not use activity labels.
"""
import argparse, json, hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

NF = [
 "NFkB_max", "TNFa_max", "RELA_max", "NFKB1_max", "NFKB2_max",
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB_max",
 "GILMORE_CORE_NFKB_PATHWAY_max",
 "BIOCARTA_NFKB_PATHWAY_max", "BIOCARTA_RELA_PATHWAY_max"
]
CELLS = ["MCF7", "U2OS", "HepRG"]

def fingerprint(smiles, gen):
    mol = Chem.MolFromSmiles(str(smiles))
    return None if mol is None else gen.GetFingerprint(mol)

def weighted_k5(q, fps, values):
    sims = np.array(DataStructs.BulkTanimotoSimilarity(q, fps), dtype=float)
    ix = np.argsort(-sims, kind="stable")[:5]
    weights = np.square(sims[ix])
    v = values[ix, :]
    finite = np.isfinite(v)
    w = weights[:, None] * finite
    denom = w.sum(axis=0)
    result = np.full(v.shape[1], np.nan)
    good = denom > 0
    result[good] = (np.nan_to_num(v) * w).sum(axis=0)[good] / denom[good]
    return result, float(sims.max()), float(sims[ix].mean())

def run(args):
    nonnf = json.loads(Path(args.non_nfkb_columns_json).read_text())
    if not isinstance(nonnf, list) or len(nonnf) != 9 or len(set(nonnf)) != 9:
        raise ValueError("Provide a JSON list of EXACTLY nine distinct non-NFkB EPA columns")
    matrix = pd.read_csv(args.matrix)
    epa = pd.read_csv(args.epa_response)
    mp = pd.read_csv(args.epa_smiles)
    if "dsstox_substance_id" in mp:
        mp = mp.rename(columns={"dsstox_substance_id":"outcome_id"})
    required = ["outcome_id", "SMILES"]
    if not set(required).issubset(mp):
        raise ValueError("EPA mapping must contain outcome_id and SMILES")
    needed = ["outcome_id","cell_type"] + NF + nonnf
    missing = set(needed) - set(epa)
    if missing:
        raise ValueError("Original EPA table missing response columns: " + repr(sorted(missing)))
    if epa.empty:
        raise ValueError("Empty EPA response table")
    mp = mp[["outcome_id","SMILES"]].dropna().drop_duplicates("outcome_id")
    epa = epa.merge(mp,on="outcome_id",how="inner",validate="many_to_one")
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
    qfps=[fingerprint(s,gen) for s in matrix.Standardized_SMILES]
    if any(x is None for x in qfps):
        raise ValueError("Invalid NF-kB SMILES encountered")
    output=matrix[["CID","InChIKey"]].copy()
    manifest={"seed":args.seed,"permutation":"EPA compound response vectors permuted as intact rows within cell line",
              "NFkB_features":NF,"non_NFkB_features":nonnf,"k":5,
              "similarity":"Morgan radius 2/2048 Tanimoto","weight":"similarity squared",
              "cell_lines":CELLS,"n_queries":len(matrix),"epa_unique_by_cell":{}}
    rng=np.random.default_rng(args.seed)
    for cell in CELLS:
        sub=epa[epa.cell_type.astype(str).str.upper()==cell.upper()].copy()
        agg=sub.groupby(["outcome_id","SMILES"],dropna=False)[NF+nonnf].mean().reset_index()
        agg=agg[agg.SMILES.notna()].copy()
        fps=[fingerprint(s,gen) for s in agg.SMILES]
        keep=[i for i,v in enumerate(fps) if v is not None]
        fps=[fps[i] for i in keep]
        agg=agg.iloc[keep].reset_index(drop=True)
        if len(fps)<5:
            raise ValueError("Fewer than 5 valid EPA chemicals in " + cell)
        manifest["epa_unique_by_cell"][cell]=len(fps)
        nfv=agg[NF].apply(pd.to_numeric,errors="coerce").to_numpy(float)
        non=agg[nonnf].apply(pd.to_numeric,errors="coerce").to_numpy(float)
        shuffled=nfv[rng.permutation(len(nfv)),:]
        sh, nn = [], []
        for q in qfps:
            a,m,mean=weighted_k5(q,fps,shuffled)
            b,_,_=weighted_k5(q,fps,non)
            sh.append(a);nn.append(b)
        sh=np.asarray(sh)
        nn=np.asarray(nn)
        for i,name in enumerate(NF):
            output["SHUF_"+cell+"_"+name+"_kNN"]=sh[:,i]
        for i,name in enumerate(nonnf):
            output["NONNF_"+cell+"_"+name+"_kNN"]=nn[:,i]
        print(cell,":",len(fps),"EPA chemicals processed",flush=True)
    out=Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    output.to_csv(out,index=False)
    (out.parent/"epa_control_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print("Control features exported:",out)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--matrix",required=True)
    p.add_argument("--epa-response",required=True)
    p.add_argument("--epa-smiles",required=True)
    p.add_argument("--non-nfkb-columns-json",required=True)
    p.add_argument("--out",required=True)
    p.add_argument("--seed",type=int,default=42)
    run(p.parse_args())
