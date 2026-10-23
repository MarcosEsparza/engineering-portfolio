"""Turn raw H408 readings into friction factors and comparisons."""

import numpy as np
import pandas as pd
from scipy.stats import linregress

from friction import core

# Relative roughness of the H408 roughened test pipe (ks/d) from the lab handout
ROUGH_PIPE_KS_D = 0.036


def load_raw(path):
    """Load the raw lab data CSV."""
    return pd.read_csv(path)


def compute_results(raw):
    """Add velocity, Re, head loss, friction factors, and uncertainty columns."""
    df = raw.copy()
    volume_m3 = df["volume_L"] / 1000
    df["head_loss_m"] = (df["upstream_mm"] - df["downstream_mm"]) / 1000
    df["velocity_m_s"] = core.mean_velocity(volume_m3, df["time_s"], df["diameter_m"])
    df["reynolds"] = core.reynolds(df["velocity_m_s"], df["diameter_m"])
    df["f_measured"] = core.darcy_friction_factor(
        df["head_loss_m"], df["velocity_m_s"], df["diameter_m"], df["length_m"])
    df["f_uncertainty_rel"] = core.friction_factor_uncertainty(
        df["head_loss_m"], df["time_s"], volume_m3)

    rel_roughness = np.where(df["pipe"].str.contains("rough"), ROUGH_PIPE_KS_D, 0.0)
    df["f_blasius"] = core.blasius(df["reynolds"])
    df["f_colebrook"] = [core.colebrook(re, r) for re, r in zip(df["reynolds"], rel_roughness)]
    df["ratio_to_colebrook"] = df["f_measured"] / df["f_colebrook"]
    return df


def head_loss_exponent(group):
    """Fit h = k * u^n on log-log axes and return (n, r^2).

    Theory: n ~ 1.75 for smooth turbulent flow, n ~ 2 for fully rough flow.
    """
    fit = linregress(np.log(group["velocity_m_s"]), np.log(group["head_loss_m"]))
    return fit.slope, fit.rvalue**2


def summarize(results):
    """One row per pipe: average agreement with theory and fitted exponent."""
    rows = []
    for pipe, group in results.groupby("pipe", sort=False):
        n, r2 = head_loss_exponent(group)
        rows.append({
            "pipe": pipe,
            "re_min": group["reynolds"].min(),
            "re_max": group["reynolds"].max(),
            "f_measured_mean": group["f_measured"].mean(),
            "f_theory_mean": group["f_colebrook"].mean(),
            "ratio_mean": group["ratio_to_colebrook"].mean(),
            "uncertainty_max_pct": 100 * group["f_uncertainty_rel"].max(),
            "exponent_n": n,
            "exponent_r2": r2,
        })
    return pd.DataFrame(rows)
