"""Propensity-score estimators with survey weights, balance diagnostics, design-based bootstrap, E-values."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, LogisticRegression

NUMERIC = ["AGEY1X", "EDUCYR", "baseline_spend_log", "opioid_fills_2021_any", "n_other_conditions",
           "VPCS2", "VMCS2", "ADPAIN2", "ADGENH2", "msk_oa", "msk_joint", "msk_spine", "msk_soft_tissue",
           "diabetes", "chd", "smoker"]
CATEGORICAL = ["SEX", "RACETHX", "POVCATY1", "INSCOVY1", "REGIONY1", "EMPST2"]


def design_matrix(df: pd.DataFrame) -> pd.DataFrame:
    d = pd.DataFrame(index=df.index)
    d["AGEY1X"] = df.AGEY1X
    d["EDUCYR"] = df.EDUCYR
    d["baseline_spend_log"] = np.log1p(df.baseline_spend_ex_pt)
    d["opioid_fills_2021_any"] = (df.opioid_fills_2021 > 0).astype(float)
    d["n_other_conditions"] = df.n_other_conditions
    for c in ["VPCS2", "VMCS2", "ADPAIN2", "ADGENH2", "msk_oa", "msk_joint", "msk_spine", "msk_soft_tissue"]:
        d[c] = df[c]
    d["diabetes"] = (df.DIABDXY1_M18 == 1).astype(float).where(df.DIABDXY1_M18.notna())
    d["chd"] = (df.CHDDXY1 == 1).astype(float).where(df.CHDDXY1.notna())
    d["smoker"] = (df.ADSMOK2 == 1).astype(float).where(df.ADSMOK2.notna())
    for c in NUMERIC:
        if d[c].isna().any():
            d[f"{c}_missing"] = d[c].isna().astype(float)
            d[c] = d[c].fillna(d[c].median())
    for c in CATEGORICAL:
        dummies = pd.get_dummies(df[c].fillna(-9).astype(int), prefix=c, drop_first=True, dtype=float)
        d = pd.concat([d, dummies], axis=1)
    return d


def propensity(X: pd.DataFrame, t: np.ndarray) -> np.ndarray:
    Z = (X - X.mean()) / X.std(ddof=0).replace(0, 1)
    m = LogisticRegression(C=1.0, max_iter=5000)
    m.fit(Z, t)
    return m.predict_proba(Z)[:, 1]


def weights(e: np.ndarray, t: np.ndarray, kind: str) -> np.ndarray:
    if kind == "overlap":
        return np.where(t == 1, 1 - e, e)
    if kind == "ipw_ate":
        w = np.where(t == 1, 1 / e, 1 / (1 - e))
        return np.where((e >= 0.05) & (e <= 0.95), w, 0.0)
    raise ValueError(kind)


def smd(X: pd.DataFrame, t: np.ndarray, w: np.ndarray) -> pd.Series:
    out = {}
    for c in X.columns:
        x = X[c].to_numpy(float)
        m1 = np.average(x[t == 1], weights=w[t == 1])
        m0 = np.average(x[t == 0], weights=w[t == 0])
        v1, v0 = x[t == 1].var(), x[t == 0].var()  # unweighted pooled SD in the denominator (standard practice)
        out[c] = (m1 - m0) / np.sqrt((v1 + v0) / 2) if (v1 + v0) > 0 else 0.0
    return pd.Series(out)


def weighted_diff(y: np.ndarray, t: np.ndarray, w: np.ndarray) -> tuple[float, float, float]:
    ok = ~np.isnan(y) & (w > 0)
    m1 = np.average(y[ok & (t == 1)], weights=w[ok & (t == 1)])
    m0 = np.average(y[ok & (t == 0)], weights=w[ok & (t == 0)])
    return m1 - m0, m1, m0


def aipw(y: np.ndarray, t: np.ndarray, X: pd.DataFrame, e: np.ndarray, sw: np.ndarray, binary: bool) -> float:
    """Doubly robust ATE (survey-weighted), PS trimmed to [0.05, 0.95]."""
    ok = ~np.isnan(y) & (e >= 0.05) & (e <= 0.95)
    y, t, Xo, e, sw = y[ok], t[ok], X[ok].to_numpy(float), e[ok], sw[ok]
    mu = {}
    for arm in (0, 1):
        idx = t == arm
        if binary:
            m = LogisticRegression(C=1.0, max_iter=5000).fit(Xo[idx], y[idx].astype(int), sample_weight=sw[idx])
            mu[arm] = m.predict_proba(Xo)[:, 1]
        else:
            m = LinearRegression().fit(Xo[idx], y[idx], sample_weight=sw[idx])
            mu[arm] = m.predict(Xo)
    psi = mu[1] - mu[0] + t * (y - mu[1]) / e - (1 - t) * (y - mu[0]) / (1 - e)
    return float(np.average(psi, weights=sw))


def evalue(rr: float) -> float:
    """VanderWeele & Ding E-value for a risk ratio (or ratio of means for a positive outcome)."""
    if not np.isfinite(rr) or rr <= 0:
        return float("nan")
    rr = 1 / rr if rr < 1 else rr
    return float(rr + np.sqrt(rr * (rr - 1)))


def psu_bootstrap_indices(df: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Resample PSUs with replacement within each variance stratum (n_h - 1 rescaling not applied)."""
    idx = []
    for _, g in df.groupby("VARSTR"):
        psus = g.VARPSU.unique()
        pick = rng.choice(psus, size=len(psus), replace=True)
        for p in pick:
            idx.extend(g.index[g.VARPSU == p].tolist())
    return np.array(idx)


def taylor_mean(y: np.ndarray, w: np.ndarray, strata: np.ndarray, psu: np.ndarray) -> tuple[float, float]:
    """Survey-weighted mean with Taylor-linearized SE (with-replacement PSU approximation)."""
    ok = ~np.isnan(y)
    y, w, strata, psu = y[ok], w[ok], strata[ok], psu[ok]
    mean = np.sum(w * y) / np.sum(w)
    z = w * (y - mean) / np.sum(w)
    var = 0.0
    frame = pd.DataFrame({"s": strata, "p": psu, "z": z}).groupby(["s", "p"]).z.sum().reset_index()
    for _, g in frame.groupby("s"):
        n = len(g)
        if n > 1:
            var += n / (n - 1) * np.sum((g.z - g.z.mean()) ** 2)
    return float(mean), float(np.sqrt(var))
