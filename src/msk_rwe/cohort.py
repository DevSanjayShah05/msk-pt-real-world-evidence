"""Build the analytic cohort from MEPS public use files (see docs/ANALYSIS_PLAN.md)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "meps"
PROC = ROOT / "data" / "processed"
PANEL = 26
OPIOID_CLASSES = {60, 191}  # Multum: narcotic analgesics, narcotic analgesic combinations
TC_FIELDS = ["TC1", "TC1S1", "TC1S1_1", "TC1S1_2", "TC2", "TC2S1", "TC2S1_1", "TC2S1_2", "TC3", "TC3S1", "TC3S1_1"]

MSK_GROUPS = {  # 3-character ICD-10 ranges (MEPS releases condition codes at 3 characters)
    "msk_oa": [f"M{n}" for n in range(15, 20)],
    "msk_joint": ["M25"],
    "msk_spine": [f"M{n}" for n in range(43, 55)],
    "msk_soft_tissue": [f"M{n}" for n in range(60, 80)],
}

LONG_VARS = ["DUPERSID", "PANEL", "LONGWT", "VARSTR", "VARPSU", "YEARIND", "AGEY1X", "SEX", "RACETHX", "POVCATY1",
             "INSCOVY1", "REGIONY1", "EDUCYR", "EMPST2", "TOTEXPY1", "TOTEXPY2", "ERTOTY1", "ERTOTY2", "IPDISY1",
             "IPDISY2", "RXEXPY2", "DVTEXPY1", "DVTEXPY2", "VPCS2", "VPCS4", "VMCS2", "VMCS4", "ADPAIN2", "ADPAIN4",
             "ADGENH2", "ADSMOK2", "DIABDXY1_M18", "CHDDXY1", "ARTHDXY1"]


def _neg_to_nan(s: pd.Series) -> pd.Series:
    s = pd.to_numeric(s, errors="coerce")
    return s.where(s >= 0)


def _opioid_fills(path: Path, panel: int) -> pd.Series:
    cols = ["DUPERSID", "PANEL"] + TC_FIELDS
    rx = pd.read_stata(path, columns=cols, convert_categoricals=False)
    rx = rx[rx.PANEL == panel]
    is_op = rx[TC_FIELDS].isin(OPIOID_CLASSES).any(axis=1)
    return rx[is_op].groupby("DUPERSID").size()


def build() -> pd.DataFrame:
    lg = pd.read_stata(RAW / "h244.dta", columns=LONG_VARS, convert_categoricals=False)
    lg = lg[(lg.PANEL == PANEL) & (lg.LONGWT > 0)].copy()

    cond = pd.read_stata(RAW / "h231.dta", columns=["DUPERSID", "PANEL", "ICD10CDX", "CCSR1X"], convert_categoricals=False)
    cond = cond[cond.PANEL == PANEL]
    for g, codes in MSK_GROUPS.items():
        lg[g] = lg.DUPERSID.isin(cond.loc[cond.ICD10CDX.isin(codes), "DUPERSID"]).astype(int)
    lg["msk_pain"] = lg[list(MSK_GROUPS)].max(axis=1)
    non_msk = cond[~cond.ICD10CDX.str.startswith("M") & cond.CCSR1X.str.match(r"^[A-Z]{3}\d{3}$", na=False)]
    lg["n_other_conditions"] = lg.DUPERSID.map(non_msk.groupby("DUPERSID").CCSR1X.nunique()).fillna(0)

    visits = []
    for f, xp in (("h229g.dta", "OBXP21X"), ("H229F.dta", "OPXP21X")):
        v = pd.read_stata(RAW / f, columns=["DUPERSID", "PANEL", "MEDPTYPE_M18", xp], convert_categoricals=False)
        v = v[(v.PANEL == PANEL) & (v.MEDPTYPE_M18 == 10)].rename(columns={xp: "xp"})
        visits.append(v)
    pt = pd.concat(visits)
    lg["pt_visits_2021"] = lg.DUPERSID.map(pt.groupby("DUPERSID").size()).fillna(0)
    lg["pt_spend_2021"] = lg.DUPERSID.map(pt.groupby("DUPERSID").xp.sum()).fillna(0)
    lg["pt"] = (lg.pt_visits_2021 >= 1).astype(int)

    lg["opioid_fills_2021"] = lg.DUPERSID.map(_opioid_fills(RAW / "h229a.dta", PANEL)).fillna(0)
    lg["opioid_fills_2022"] = lg.DUPERSID.map(_opioid_fills(RAW / "h239a.dta", PANEL)).fillna(0)

    for c in ["VPCS2", "VPCS4", "VMCS2", "VMCS4", "ADPAIN2", "ADPAIN4", "ADGENH2", "ADSMOK2", "EDUCYR", "EMPST2",
              "POVCATY1", "INSCOVY1", "REGIONY1", "DIABDXY1_M18", "CHDDXY1", "ARTHDXY1"]:
        lg[c] = _neg_to_nan(lg[c])
    lg["baseline_spend_ex_pt"] = (lg.TOTEXPY1 - lg.pt_spend_2021).clip(lower=0)

    df = lg[(lg.msk_pain == 1) & lg.AGEY1X.between(18, 84)].copy()
    # outcomes
    df["y_total_spend"] = df.TOTEXPY2
    df["y_any_opioid"] = (df.opioid_fills_2022 > 0).astype(int)
    df["y_opioid_fills"] = df.opioid_fills_2022
    df["y_any_er"] = (df.ERTOTY2 > 0).astype(int)
    df["y_any_inpatient"] = (df.IPDISY2 > 0).astype(int)
    df["y_pcs"] = df.VPCS4
    df["y_pain_interference"] = df.ADPAIN4
    df["nc_dental_spend"] = df.DVTEXPY2
    PROC.mkdir(parents=True, exist_ok=True)
    df.to_parquet(PROC / "cohort.parquet", index=False)
    lg.to_parquet(PROC / "panel26_adults_all.parquet", index=False)
    return df


if __name__ == "__main__":
    d = build()
    print(f"cohort n={len(d)}  PT users={int(d.pt.sum())}  weighted PT share="
          f"{np.average(d.pt, weights=d.LONGWT):.3f}")
    print(d[list(MSK_GROUPS)].mean().round(3).to_dict())
