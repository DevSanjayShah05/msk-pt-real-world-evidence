"""Static figures for the README (matplotlib, reference palette)."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES, FIG = ROOT / "results", ROOT / "docs" / "figures"
BLUE, ORANGE, AQUA, GRAY, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#b4b2a9", "#52514e"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": "#cfcdc6", "axes.labelcolor": INK2, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})


def love_plot():
    b = pd.read_csv(RES / "balance_smd.csv", index_col=0)
    b = b.loc[b.unweighted.abs().sort_values().index]
    b = b.iloc[-25:]
    fig, ax = plt.subplots(figsize=(6.5, 7))
    y = np.arange(len(b))
    ax.scatter(b.unweighted.abs(), y, color=GRAY, s=28, label="Unweighted", zorder=3)
    ax.scatter(b.overlap.abs(), y, color=BLUE, s=28, label="Overlap-weighted", zorder=3)
    ax.axvline(0.1, color=ORANGE, ls="--", lw=1)
    ax.text(0.102, len(b) - 0.5, "0.10", color=ORANGE, fontsize=9)
    ax.set_yticks(y, b.index, fontsize=8)
    ax.set_xlabel("|Standardized mean difference|, PT vs no PT")
    ax.set_title("Covariate balance (25 most imbalanced before weighting)", fontsize=11, loc="left")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "balance_love_plot.png", dpi=160)


def ps_overlap():
    ps = pd.read_csv(RES / "propensity_scores.csv")
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    bins = np.linspace(0, 1, 26)
    ax.hist(ps.ps[ps.pt == 0], bins=bins, color=GRAY, alpha=0.9, label="No PT", density=True)
    ax.hist(ps.ps[ps.pt == 1], bins=bins, color=BLUE, alpha=0.75, label="PT", density=True)
    ax.set_xlabel("Estimated propensity to use PT")
    ax.set_ylabel("Density")
    ax.set_title("Propensity score overlap", fontsize=11, loc="left")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "propensity_overlap.png", dpi=160)


def forest():
    r = pd.read_csv(RES / "effects.csv")
    order = ["y_total_spend", "nc_dental_spend", "y_any_opioid", "y_any_er", "y_any_inpatient", "y_opioid_fills",
             "y_pcs", "y_pain_interference"]
    panels = [("Dollars", ["y_total_spend", "nc_dental_spend"], 1, "Difference in 2022 spending ($)"),
              ("Probabilities", ["y_any_opioid", "y_any_er", "y_any_inpatient"], 100, "Difference (percentage points)"),
              ("Counts and scores", ["y_opioid_fills", "y_pcs", "y_pain_interference"], 1, "Difference (units of each outcome)")]
    fig, axes = plt.subplots(3, 1, figsize=(7, 7.5), gridspec_kw={"height_ratios": [2, 3, 3]})
    colors = {"overlap": BLUE, "ipw_ate": ORANGE, "aipw": AQUA}
    names = {"overlap": "Overlap weights (primary)", "ipw_ate": "IPW, ATE", "aipw": "AIPW (doubly robust)"}
    for ax, (title, outs, scale, xlabel) in zip(axes, panels):
        for i, o in enumerate(outs):
            for j, est in enumerate(["overlap", "ipw_ate", "aipw"]):
                row = r[(r.outcome == o) & (r.estimator == est)].iloc[0]
                yy = i + (j - 1) * 0.22
                ax.errorbar(row["diff"] * scale, yy, xerr=[[(row["diff"] - row.ci_low) * scale],
                            [(row.ci_high - row["diff"]) * scale]], fmt="o", color=colors[est], ms=5, capsize=2,
                            label=names[est] if (i == 0 and ax is axes[0]) else None)
        ax.axvline(0, color=INK2, lw=0.8)
        ax.set_yticks(range(len(outs)), [r[r.outcome == o].label.iloc[0].replace("NEGATIVE CONTROL: ", "NEG. CONTROL: ")
                                          for o in outs], fontsize=8)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_title(title, fontsize=10, loc="left")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG / "effects_forest.png", dpi=160)


if __name__ == "__main__":
    FIG.mkdir(parents=True, exist_ok=True)
    love_plot()
    ps_overlap()
    forest()
    print("figures written")
