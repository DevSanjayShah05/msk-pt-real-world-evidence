import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from msk_rwe.causal import evalue, smd, taylor_mean, weighted_diff, weights
from msk_rwe.cohort import MSK_GROUPS

ROOT = Path(__file__).resolve().parents[1]


def test_evalue_matches_published_formula():
    # VanderWeele & Ding (2017): RR = 3.9 -> E-value 7.26
    assert evalue(3.9) == pytest.approx(7.26, abs=0.01)
    assert evalue(1 / 3.9) == pytest.approx(7.26, abs=0.01)
    assert evalue(1.0) == pytest.approx(1.0)


def test_overlap_weights_are_bounded_and_ipw_trims():
    e = np.array([0.01, 0.2, 0.5, 0.8, 0.99])
    t = np.array([1, 1, 0, 0, 1])
    ow = weights(e, t, "overlap")
    assert np.all((ow >= 0) & (ow <= 1))
    iw = weights(e, t, "ipw_ate")
    assert iw[0] == 0 and iw[-1] == 0 and iw[1] == pytest.approx(5.0)


def test_smd_zero_when_groups_identical():
    X = pd.DataFrame({"a": [1.0, 2.0, 1.0, 2.0]})
    assert smd(X, np.array([1, 1, 0, 0]), np.ones(4))["a"] == pytest.approx(0.0)


def test_weighted_diff_ignores_missing_outcomes():
    y = np.array([1.0, np.nan, 3.0, 5.0])
    d, m1, m0 = weighted_diff(y, np.array([1, 1, 0, 0]), np.ones(4))
    assert m1 == 1.0 and m0 == 4.0 and d == -3.0


def test_taylor_mean_reduces_to_iid_se_with_one_psu_per_unit():
    rng = np.random.default_rng(0)
    y = rng.normal(10, 2, 2000)
    w = np.ones(2000)
    m, se = taylor_mean(y, w, np.zeros(2000), np.arange(2000))
    assert m == pytest.approx(y.mean())
    assert se == pytest.approx(y.std(ddof=1) / np.sqrt(2000), rel=0.02)


def test_msk_definition_excludes_ra_gout_osteoporosis():
    codes = {c for v in MSK_GROUPS.values() for c in v}
    assert {"M54", "M17", "M25", "M79", "M51"} <= codes
    assert not {"M05", "M06", "M10", "M81"} & codes


def test_committed_results_are_coherent():
    s = json.loads((ROOT / "results" / "summary.json").read_text())
    assert s["n_pt"] < s["n_cohort"]
    assert s["max_abs_smd"]["overlap"] < s["max_abs_smd"]["unweighted"]
    eff = pd.DataFrame(s["effects"])
    assert set(eff.estimator) == {"overlap", "ipw_ate", "aipw"}
    assert (eff.ci_low <= eff["diff"]).all() and (eff["diff"] <= eff.ci_high).all()
