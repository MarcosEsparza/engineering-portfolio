"""Run the full analysis: read raw data, compute results, save tables and plots.

Usage:
    python analyze.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from friction import core
from friction.analysis import ROUGH_PIPE_KS_D, compute_results, load_raw, summarize

HERE = Path(__file__).parent
RESULTS = HERE / "results"
COLORS = {"13.6 mm smooth": "#2a78d6", "26.2 mm smooth": "#e3812f", "14 mm rough": "#4a9b4f"}


def plot_friction_factor(results, path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for pipe, group in results.groupby("pipe", sort=False):
        ax.errorbar(group["reynolds"], group["f_measured"],
                    yerr=group["f_measured"] * group["f_uncertainty_rel"],
                    fmt="o", capsize=3, color=COLORS[pipe], label=f"{pipe} (measured)")

    re = np.logspace(np.log10(3000), np.log10(2e4), 100)
    ax.plot(re, core.blasius(re), "k-", label="Blasius (smooth)")
    ax.plot(re, [core.colebrook(r, ROUGH_PIPE_KS_D) for r in re], "k--",
            label=f"Colebrook, ks/d = {ROUGH_PIPE_KS_D}")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Reynolds number, Re")
    ax.set_ylabel("Darcy friction factor, f")
    ax.set_title("H408 friction factor vs. Reynolds number")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8, loc="center right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_head_loss(results, summary, path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    exponents = summary.set_index("pipe")["exponent_n"]
    for pipe, group in results.groupby("pipe", sort=False):
        ax.loglog(group["velocity_m_s"], group["head_loss_m"], "o", color=COLORS[pipe],
                  label=f"{pipe}: h ∝ u^{exponents[pipe]:.2f}")
    ax.set_xlabel("Mean velocity, u (m/s)")
    ax.set_ylabel("Head loss, h (m)")
    ax.set_title("Head loss scaling (theory: n ≈ 1.75 smooth, n ≈ 2 rough)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    RESULTS.mkdir(exist_ok=True)
    results = compute_results(load_raw(HERE / "data" / "h408_raw.csv"))
    summary = summarize(results)

    results.round(5).to_csv(RESULTS / "results.csv", index=False)
    summary.round(4).to_csv(RESULTS / "summary.csv", index=False)
    plot_friction_factor(results, RESULTS / "friction_factor_vs_re.png")
    plot_head_loss(results, summary, RESULTS / "head_loss_scaling.png")

    print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
