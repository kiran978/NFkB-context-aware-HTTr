#!/usr/bin/env python3
"""Selectively extract EPA response controls from the ORIGINAL ~394 MB Feather.

Reads Arrow schema first and only desired response columns, avoiding a full
17,029-column pandas allocation. Does NOT read NF-kB activity labels.
First, run --scan to list candidate non-NFkB pathways. THEN choose and
document exactly nine exact names in a JSON array BEFORE model fitting.
Second, run --extract with the nine-column JSON to export compact CSV.

Use the already-existing local chemical_httr_assay_aggregated.feather.
Feather V2/Arrow IPC is expected. Requires: pip install pyarrow pandas.
"""
import argparse, json, re
from pathlib import Path
import pyarrow as pa
import pyarrow.ipc as ipc
import pyarrow.feather as feather

NF = [
 "NFkB_max","TNFa_max","RELA_max","NFKB1_max","NFKB2_max",
 "HALLMARK_TNFA_SIGNALING_VIA_NFKB_max","GILMORE_CORE_NFKB_PATHWAY_max",
 "BIOCARTA_NFKB_PATHWAY_max","BIOCARTA_RELA_PATHWAY_max"
]
PATTERNS = {
 "JAK_STAT":r"JAK|STAT",
 "MAPK":r"MAPK|ERK|P38",
 "PI3K_AKT":r"PI3K|AKT|MTOR",
 "JUN_FOS_AP1":r"JUN|FOS|AP1",
 "NFE2L2_NRF2":r"NFE2L2|NRF2",
}
def run(args):
    with pa.memory_map(str(args.feather),"r") as f:
        try:
            names = list(ipc.open_file(f).schema.names)
        except Exception as e:
            raise RuntimeError("Unable to read Feather V2 schema: " + str(e))
    candidates = {}
    for pathway, pattern in PATTERNS.items():
        candidates[pathway] = [x for x in names
                  if x.endswith("_max") and re.search(pattern,x,re.IGNORECASE)
                  and x not in NF]
    print("Total columns:",len(names))
    print("NF-kB source columns present:",sum(c in names for c in NF),"/",len(NF))
    for label, matches in candidates.items():
        print("\n["+label+"]",len(matches))
        for x in matches[:100]:
            print("  "+x)
    Path(args.output_dir).mkdir(parents=True,exist_ok=True)
    (Path(args.output_dir)/"non_nfkb_candidate_names.json").write_text(
      json.dumps(candidates,indent=2)+"\n",encoding="utf-8")
    if not args.select:
        print("\nScan complete. Choose EXACTLY 9 appropriate columns in a JSON list;")
        print("then rerun with --select to extract a small, scientifically documented matrix.")
        return
    nfn=json.loads(Path(args.select).read_text(encoding="utf-8"))
    if not isinstance(nfn,list) or len(nfn)!=9 or len(set(nfn))!=9:
        raise ValueError("--select must specify exactly 9 distinct feature names")
    select=["outcome_id","cell_type"]+NF+nfn
    absent=sorted(set(select)-set(names))
    if absent:
        raise ValueError("Selected columns absent from original EPA Feather: "+repr(absent))
    t=feather.read_table(str(args.feather),columns=select)
    df=t.to_pandas()
    path=Path(args.output_dir)/"EPA_NFkB_and_NonNFkB_Controls_Compact.csv"
    df.to_csv(path,index=False)
    print("Saved:",path,"rows=",len(df),"columns=",len(df.columns))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--feather",required=True)
    p.add_argument("--output-dir",required=True)
    p.add_argument("--select",help="JSON list of exactly nine named non-NFkB response columns")
    run(p.parse_args())
