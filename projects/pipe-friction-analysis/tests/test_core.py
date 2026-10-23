import numpy as np
import pandas as pd
import pytest

from friction import core
from friction.analysis import compute_results, head_loss_exponent


def test_velocity_matches_lab_hand_calculation():
    # Report trial 1: 15 L in 83 s through the 13.6 mm pipe gives about 1.24 m/s
    u = core.mean_velocity(15e-3, 83, 0.0136)
    assert u == pytest.approx(1.24, abs=0.01)


def test_reynolds_matches_lab_hand_calculation():
    # Report: Re = 16,848 using the unrounded velocity from trial 1
    u = core.mean_velocity(15e-3, 83, 0.0136)
    assert core.reynolds(u, 0.0136) == pytest.approx(16848, rel=1e-3)


def test_darcy_factor_inverts_darcy_weisbach():
    f, d, length, u = 0.03, 0.02, 1.0, 1.5
    h = f * (length / d) * u**2 / (2 * core.G)
    assert core.darcy_friction_factor(h, u, d, length) == pytest.approx(f)


def test_blasius_known_value():
    # 0.316 * 10,000^-0.25 = 0.0316
    assert core.blasius(1e4) == pytest.approx(0.0316, rel=1e-3)


@pytest.mark.parametrize("re, rel_roughness, expected", [
    (1e5, 0.0, 0.0180),     # smooth pipe, Moody chart
    (1e5, 0.001, 0.0222),   # moderately rough
    (1e6, 0.01, 0.0380),    # nearly fully rough
])
def test_colebrook_matches_moody_chart(re, rel_roughness, expected):
    assert core.colebrook(re, rel_roughness) == pytest.approx(expected, rel=0.01)


def test_colebrook_close_to_blasius_for_smooth_pipe():
    # The two smooth-pipe correlations agree within a few percent in Blasius's range
    for re in (5e3, 2e4, 5e4):
        assert core.colebrook(re, 0.0) == pytest.approx(core.blasius(re), rel=0.05)


def test_uncertainty_grows_as_head_loss_shrinks():
    big = core.friction_factor_uncertainty(0.25, 80, 15e-3)
    small = core.friction_factor_uncertainty(0.01, 80, 15e-3)
    assert small > big


def test_head_loss_exponent_recovers_power_law():
    u = np.array([0.2, 0.4, 0.8, 1.6])
    group = pd.DataFrame({"velocity_m_s": u, "head_loss_m": 0.05 * u**1.75})
    n, r2 = head_loss_exponent(group)
    assert n == pytest.approx(1.75)
    assert r2 == pytest.approx(1.0)


def test_compute_results_adds_expected_columns():
    raw = pd.DataFrame([{
        "pipe": "13.6 mm smooth", "diameter_m": 0.0136, "length_m": 0.912,
        "volume_L": 15, "time_s": 83, "upstream_mm": 628, "downstream_mm": 379,
    }])
    results = compute_results(raw)
    row = results.iloc[0]
    assert row["head_loss_m"] == pytest.approx(0.249)
    assert row["f_colebrook"] > 0
    assert 0 < row["f_uncertainty_rel"] < 1
