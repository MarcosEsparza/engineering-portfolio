"""Core fluid-mechanics relations used in the analysis.

All friction factors are Darcy friction factors (the convention used on the
Moody chart): h = f * (L/d) * u^2 / (2g). The Fanning factor is f / 4.
"""

import numpy as np
from scipy.optimize import brentq

G = 9.81            # gravitational acceleration [m/s^2]
NU_WATER = 1.004e-6  # kinematic viscosity of water at ~20 °C [m^2/s]


def pipe_area(d):
    """Cross-sectional area of a circular pipe [m^2]."""
    return np.pi * d**2 / 4


def mean_velocity(volume_m3, time_s, d):
    """Mean velocity from a timed volumetric collection [m/s]."""
    return volume_m3 / time_s / pipe_area(d)


def reynolds(u, d, nu=NU_WATER):
    """Reynolds number Re = u d / nu."""
    return u * d / nu


def darcy_friction_factor(head_loss_m, u, d, length_m):
    """Darcy friction factor from measured head loss (Darcy-Weisbach)."""
    return 2 * G * head_loss_m * d / (length_m * u**2)


def blasius(re):
    """Blasius correlation for smooth pipes, valid for roughly 4,000 < Re < 100,000."""
    return 0.316 * np.asarray(re, dtype=float) ** -0.25


def colebrook(re, rel_roughness):
    """Darcy friction factor from the Colebrook-White equation.

    Solved numerically with Brent's method. rel_roughness is ks/d
    (use 0 for a hydraulically smooth pipe).
    """
    def residual(f):
        return 1 / np.sqrt(f) + 2 * np.log10(rel_roughness / 3.7 + 2.51 / (re * np.sqrt(f)))

    return brentq(residual, 1e-4, 1.0)


def friction_factor_uncertainty(head_loss_m, time_s, volume_m3,
                                sigma_reading_m=0.001, sigma_time_s=0.5,
                                sigma_volume_m3=0.1e-3):
    """Relative uncertainty in f from the measured quantities.

    Since u = V / (t A), the friction factor scales as f ~ h t^2 / V^2.
    Independent errors combine as:
        (df/f)^2 = (dh/h)^2 + (2 dt/t)^2 + (2 dV/V)^2
    The head loss is the difference of two piezometer readings, so its
    uncertainty is sqrt(2) times the single-reading uncertainty.
    """
    sigma_h = np.sqrt(2) * sigma_reading_m
    return np.sqrt((sigma_h / head_loss_m) ** 2
                   + (2 * sigma_time_s / time_s) ** 2
                   + (2 * sigma_volume_m3 / volume_m3) ** 2)
