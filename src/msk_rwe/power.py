"""Experiment design for an MSK program, using MEPS year-1 -> year-2 correlations (pre-specified, see plan).

CUPED (Deng et al., 2013) adjusts the outcome with its pre-period value: Y_adj = Y - theta (X - mean X), which
cuts variance by a factor of (1 - rho^2). Required sample size scales the same way.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
PROC, RES = ROOT / "data" / "processed", ROOT / "results"

PAIRS = {  # outcome -> (year-2 column, year-1 column, effect sizes to power for, unit)
    "total_spend": ("TOTEXPY2", "TOTEXPY1", [500, 1000, 2000], "$ per member per year"),
    "total_spend_winsorized_p99": ("TOTEXPY2", "TOTEXPY1", [500, 1000, 2000], "$ (winsorized at p99)"),
    "pcs": ("VPCS4", "VPCS2", [1, 2, 3], "VR-12 PCS points"),
    "pain_interference": ("ADPAIN4", "ADPAIN2", [0.1, 0.2, 0.3], "points on the 1-5 scale"),
    "any_opioid": ("y_any_opioid", "opioid_fills_2021", [0.02, 0.03, 0.05], "absolute risk"),
}


def n_per_arm(sd: float, delta: float, alpha: float = 0.05, power: float = 0.8) -> int:
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return int(np.ceil(2 * (z * sd / delta) ** 2))


def wstats(x, y, w):
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    cov = np.average((x - mx) * (y - my), weights=w)
    vx, vy = np.average((x - mx) ** 2, weights=w), np.average((y - my) ** 2, weights=w)
    return np.sqrt(vy), cov / np.sqrt(vx * vy)


def run() -> dict:
    df = pd.read_parquet(PROC / "cohort.parquet")
    out = {}
    for name, (ycol, xcol, deltas, unit) in PAIRS.items():
        d = df[[ycol, xcol, "LONGWT"]].dropna()
        y, x = d[ycol].to_numpy(float), d[xcol].to_numpy(float)
        if name == "any_opioid":
            x = (x > 0).astype(float)
        if name.endswith("p99"):
            cap = np.quantile(np.r_[y, x], 0.99)
            y, x = np.minimum(y, cap), np.minimum(x, cap)
        sd, rho = wstats(x, y, d.LONGWT.to_numpy(float))
        reduction = rho ** 2
        out[name] = {"unit": unit, "n": int(len(d)), "sd": float(sd), "rho_pre_post": float(rho),
                     "variance_reduction": float(reduction),
                     "n_per_arm": {str(dl): {"unadjusted": n_per_arm(sd, dl),
                                             "cuped": n_per_arm(sd * np.sqrt(1 - reduction), dl)} for dl in deltas}}
    RES.mkdir(parents=True, exist_ok=True)
    (RES / "power_cuped.json").write_text(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
