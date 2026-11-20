#!/usr/bin/env python3
"""Reproduce the key prospective-screen counts reported in the manuscript."""
import argparse,pandas as pd

def yn(s): return s.astype(str).eq('Yes')
def main(a):
    b=pd.read_csv(a.phase8b); c=pd.read_csv(a.phase8c); d=pd.read_csv(a.final12)
    strict=yn(b.Strict_All_Context_Confidence)
    anyactive=yn(b.Bcell_HighConfidence_Active)|yn(b.TNFa489022_HighConfidence_Active)|yn(b.TNFa2337_HighConfidence_Active)
    out={
      'Prospective_library':len(b),
      'All3_structural_AD':int(b.Full_3Task_Structural_AD.astype(str).eq('IN').sum()),
      'HTTr_allcell_AD':int(b.HTTr_AllCell_AD.astype(str).eq('IN').sum()),
      'Strict_all_context_confidence':int(strict.sum()),
      'Strict_with_any_HC_active':int((strict&anyactive).sum()),
      'Tie_robust_strict':int(yn(b.TieRobust_Strict_Confidence).sum()),
      'Development_clean':int(c.Dev_Chem_Gate.astype(str).str.lower().eq('true').sum()),
      'Primary_panel':int(c.Phase8C_Panel_Status.eq('PRIMARY_VALIDATION').sum()),
      'Reserve_panel':int(c.Phase8C_Panel_Status.eq('RESERVE').sum()),
      'Final_panel':len(d)}
    pd.DataFrame(out.items(),columns=['Metric','N']).to_csv(a.out,index=False)
    print(pd.DataFrame(out.items(),columns=['Metric','N']).to_string(index=False))
if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--phase8b',required=True); p.add_argument('--phase8c',required=True); p.add_argument('--final12',required=True); p.add_argument('--out',required=True); main(p.parse_args())
