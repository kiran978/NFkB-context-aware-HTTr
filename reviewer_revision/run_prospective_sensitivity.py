#!/usr/bin/env python3
"""Retrospective sensitivity of original, frozen Phase-8B NF-kB screening gates.

NEW post-review analysis, not an alteration to the original 2000-compound
prospective screen. Does NOT retrain, rescore, resort or reselect compounds.
Original AD similarity thresholds are multiplied by 0.9, 1.0 or 1.1.
Original entropy thresholds are shifted by -0.02, 0 or +0.02 (clipped [0,1]).
EPA all-cell HTTr AD membership and tie classifications remain FIXED.

Usage:
 python reviewer_revision/run_prospective_sensitivity.py \
  --phase8b data/NFkB_Phase8B_All_2000_Predictions.csv \
  --phase8c data/NFkB_Phase8C_Prioritized_All_2000.csv \
  --final12 data/NFkB_Phase8D_Final_12_Leads.csv \
  --out reviewer_outputs/gate_sensitivity.csv
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

TASKS = ["Bcell489023", "TNFa489022", "TNFa2337"]
BCCELL_STRICT_FINAL = {"P440AM19590", "P440AM34824", "P440AM34830", "P440AM61633"}

def yes(s):
    return s.astype(str).str.strip().eq("Yes")

def run(a):
    df = pd.read_csv(a.phase8b)
    if len(df) != 2000 or df.ID.nunique() != 2000:
        raise ValueError("Expected 2000 unique prospective compounds")
    ostrict = yes(df.Strict_All_Context_Confidence)
    oactive = (yes(df.Bcell_HighConfidence_Active) |
               yes(df.TNFa489022_HighConfidence_Active) |
               yes(df.TNFa2337_HighConfidence_Active))
    otie = yes(df.TieRobust_Strict_Confidence)
    if (int(ostrict.sum()), int((ostrict & oactive).sum()), int(otie.sum())) != (199,193,172):
        raise ValueError("Original frozen counts do not reproduce exactly; abort.")
    primary, reserve, final = set(), set(), set()
    if a.phase8c:
        c = pd.read_csv(a.phase8c)
        primary = set(c.loc[c.Phase8C_Panel_Status.eq("PRIMARY_VALIDATION"), "ID"])
        reserve = set(c.loc[c.Phase8C_Panel_Status.eq("RESERVE"), "ID"])
        if (len(primary), len(reserve)) != (20, 20):
            raise ValueError("Unexpected primary/reserve sizes")
    if a.final12:
        f = pd.read_csv(a.final12)
        final = set(f.ID)
        if len(final) != 12:
            raise ValueError("Expected final 12 unique IDs")
    all_httr_ad = df.HTTr_AllCell_AD.eq("IN")
    not_tie = df.HTTr_AnyTieSensitive.eq("No")
    rows = []
    for admult in [0.9, 1.0, 1.1]:
        ad_per_task = {}
        for t in TASKS:
            ad_per_task[t] = (pd.to_numeric(df[t + "_MaxTrain_ECFP4"]) >=
                              pd.to_numeric(df[t + "_AD_Threshold"])*admult)
        all_ad = np.logical_and.reduce(list(ad_per_task.values()))
        for entropy_shift in [-0.02, 0, 0.02]:
            ent_per_task = {}
            for t in TASKS:
                thr = (pd.to_numeric(df[t + "_Entropy80_Threshold"]) +
                       entropy_shift).clip(0, 1)
                ent_per_task[t] = pd.to_numeric(df[t + "_Entropy"]) <= thr
            all_entropy = np.logical_and.reduce(list(ent_per_task.values()))
            strict = all_ad & all_entropy & all_httr_ad
            any_active = np.logical_or.reduce([
                pd.to_numeric(df[t+"_Probability"]) >= 0.5 for t in TASKS])
            active = strict & any_active
            ties = strict & not_tie
            selected_ids = set(df.loc[strict, "ID"])
            active_ids = set(df.loc[active, "ID"])
            rows.append({
                "AD_threshold_multiplier":admult,
                "Entropy_threshold_shift":entropy_shift,
                "All3_structural_AD_N":int(sum(all_ad)),
                "All3_low_entropy_N":int(sum(all_entropy)),
                "Strict_confidence_N":int(strict.sum()),
                "Strict_plus_active_N":int(active.sum()),
                "Strict_tie_robust_N":int(ties.sum()),
                "Original_strict_retained_N":len(selected_ids & set(df.loc[ostrict,"ID"])),
                "Original_primary_in_strict_N":len(selected_ids & primary),
                "Original_reserve_in_strict_N":len(selected_ids & reserve),
                "Original_final12_in_strict_N":len(selected_ids & final),
                "Original_final12_in_strict_active_N":len(active_ids & final),
                "Four_strict_Bcell_final_in_strict_N":len(selected_ids & BCCELL_STRICT_FINAL),
                "Four_strict_Bcell_final_in_strict_active_N":len(active_ids & BCCELL_STRICT_FINAL),
            })
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    result = pd.DataFrame(rows)
    result.to_csv(out, index=False)
    baseline = result[(result.AD_threshold_multiplier.eq(1.0)) &
                      (result.Entropy_threshold_shift.eq(0))].iloc[0]
    if int(baseline.Strict_confidence_N) != 199 or int(baseline.Strict_plus_active_N) != 193 \
       or int(baseline.Strict_tie_robust_N) != 172:
        raise ValueError("Retrospective baseline differs from archived Phase-8B flags")
    print(result.to_string(index=False))
    print("\nOriginal selection and probabilities have not been changed.")
    print("The 9 comparisons are multiplicative/additive cutoff sensitivity,")
    print("not reestimated training-distribution percentile thresholds.")

if __name__ == "__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--phase8b",required=True)
    p.add_argument("--phase8c")
    p.add_argument("--final12")
    p.add_argument("--out",required=True)
    run(p.parse_args())
