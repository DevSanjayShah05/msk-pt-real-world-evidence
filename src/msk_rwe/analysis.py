"""Run the pre-specified analysis (docs/ANALYSIS_PLAN.md) and write results/."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .causal import (aipw, design_matrix, evalue, propensity, psu_bootstrap_indices, smd, taylor_mean,
                     weighted_diff, weights)

ROOT = Path(__file__).resolve().parents[2]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
OUTCOMES = {  # name -> (column, binary?, label)
    "y_total_spend": ("y_total_spend", False, "2022 total health care spending ($)"),
    "y_any_opioid": ("y_any_opioid", True, "Any opioid fill in 2022"),
    "y_opioid_fills": ("y_opioid_fills", False, "Opioid fills in 2022"),
    "y_any_er": ("y_any_er", True, "Any ER visit in 2022"),
    "y_any_inpatient": ("y_any_inpatient", True, "Any inpatient stay in 2022"),
    "y_pcs": ("y_pcs", False, "VR-12 physical component score, 2022 (higher = better)"),
    "y_pain_interference": ("y_pain_interference", False, "Pain interference with work, 2022 (1-5, lower = better)"),
    "nc_dental_spend": ("nc_dental_spend", False, "NEGATIVE CONTROL: 2022 dental spending ($)"),
}
PRIMARY = ["y_total_spend", "y_any_opioid"]


def estimate(df: pd.DataFrame) -> dict:
    X = design_matrix(df)
    t = df.pt.to_numpy()
    sw = df.LONGWT.to_numpy(float)
    e = propensity(X, t)
    out = {}
    for kind in ("overlap", "ipw_ate"):
        w = sw * weights(e, t, kind)
        for name, (col, binary, _) in OUTCOMES.items():
            diff, m1, m0 = weighted_diff(df[col].to_numpy(float), t, w)
            out[(name, kind)] = {"diff": diff, "treated_mean": m1, "control_mean": m0}
    for name, (col, binary, _) in OUTCOMES.items():
        out[(name, "aipw")] = {"diff": aipw(df[col].to_numpy(float), t, X, e, sw, binary)}
    return out, e, X


def run(n_boot: int = 500, seed: int = 2026) -> dict:
    RES.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(PROC / "cohort.parquet").reset_index(drop=True)
    allad = pd.read_parquet(PROC / "panel26_adults_all.parquet")
    allad = allad[allad.AGEY1X.between(18, 84)]

    # Descriptives (survey-weighted, Taylor SE)
    desc = {}
    for label, y, d in [("msk_pain_prevalence_adults", allad.msk_pain, allad), ("pt_use_among_msk", df.pt, df),
                        ("any_opioid_2021_among_msk", (df.opioid_fills_2021 > 0).astype(int), df)]:
        m, se = taylor_mean(y.to_numpy(float), d.LONGWT.to_numpy(float), d.VARSTR.to_numpy(), d.VARPSU.to_numpy())
        desc[label] = {"estimate": m, "ci95": [m - 1.96 * se, m + 1.96 * se], "n": int(len(d))}
    pop = allad.LONGWT.sum()
    desc["us_adults_represented_millions"] = float(pop / 1e6)
    desc["us_adults_with_msk_pain_millions"] = float(df.LONGWT.sum() / 1e6)

    point, e, X = estimate(df)
    t = df.pt.to_numpy()
    sw = df.LONGWT.to_numpy(float)
    bal = pd.DataFrame({"unweighted": smd(X, t, np.ones(len(t))), "survey_weighted": smd(X, t, sw),
                        "overlap": smd(X, t, sw * weights(e, t, "overlap")),
                        "ipw_ate": smd(X, t, sw * weights(e, t, "ipw_ate"))})
    bal.to_csv(RES / "balance_smd.csv")
    ps = pd.DataFrame({"pt": t, "ps": e})
    ps.to_csv(RES / "propensity_scores.csv", index=False)

    rng = np.random.default_rng(seed)
    boots = {k: [] for k in point}
    for b in range(n_boot):
        idx = psu_bootstrap_indices(df, rng)
        dfb = df.loc[idx].reset_index(drop=True)
        if dfb.pt.nunique() < 2:
            continue
        pb, _, _ = estimate(dfb)
        for k, v in pb.items():
            boots[k].append(v["diff"])

    rows = []
    for (name, kind), v in point.items():
        col, binary, label = OUTCOMES[name]
        bs = np.array(boots[(name, kind)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        row = {"outcome": name, "label": label, "estimator": kind, "binary": binary, "diff": v["diff"],
               "ci_low": lo, "ci_high": hi, "treated_mean": v.get("treated_mean"), "control_mean": v.get("control_mean"),
               "n_nonmissing": int(df[col].notna().sum()), "ci_excludes_zero": bool(lo > 0 or hi < 0)}
        if v.get("control_mean"):
            ratio = v["treated_mean"] / v["control_mean"] if v["control_mean"] > 0 else np.nan
            row["ratio"] = ratio
            row["e_value"] = evalue(ratio) if name in PRIMARY + ["nc_dental_spend"] else None
        rows.append(row)
    res = pd.DataFrame(rows)
    res.to_csv(RES / "effects.csv", index=False)

    summary = {
        "n_cohort": int(len(df)), "n_pt": int(t.sum()), "n_bootstrap": len(boots[("y_total_spend", "overlap")]),
        "descriptives": desc,
        "ps_range": {"pt": [float(e[t == 1].min()), float(e[t == 1].max())],
                     "no_pt": [float(e[t == 0].min()), float(e[t == 0].max())]},
        "ipw_trimmed_n": int(((e < 0.05) | (e > 0.95)).sum()),
        "max_abs_smd": {c: float(bal[c].abs().max()) for c in bal.columns},
        "n_covariates_smd_over_0.1": {c: int((bal[c].abs() > 0.1).sum()) for c in bal.columns},
        "effects": res.to_dict("records"),
    }
    (RES / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return summary


if __name__ == "__main__":
    s = run()
    print(json.dumps({k: v for k, v in s.items() if k != "effects"}, indent=2, default=float))
    r = pd.DataFrame(s["effects"])
    print(r[["outcome", "estimator", "diff", "ci_low", "ci_high", "treated_mean", "control_mean", "e_value"]]
          .round(3).to_string(index=False))
