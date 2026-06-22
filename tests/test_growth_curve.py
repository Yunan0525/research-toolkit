# =============================================================================
# tests/test_growth_curve.py  –  Unit Tests for analysis/growth_curve.py
# =============================================================================
# HOW TO RUN
#   pytest tests/test_growth_curve.py -v
#
# PHASE 4 STATUS
#   Tests are specifications. xfail until Phase 4.
# =============================================================================

import pytest
import numpy as np
import pandas as pd

from analysis.growth_curve import fit_growth_model, compute_group_statistics


def _logistic(t, ODmax, mu_max, lag):
    """Reference logistic model for generating synthetic data."""
    return ODmax / (1 + np.exp(-mu_max * (t - lag)))


# ── Synthetic data: known parameters ─────────────────────────────────────────
TRUE_ODMAX = 1.5
TRUE_MU_MAX = 0.4
TRUE_LAG = 2.0
TIME = np.linspace(0, 24, 50)
OD_CLEAN = _logistic(TIME, TRUE_ODMAX, TRUE_MU_MAX, TRUE_LAG)


@pytest.mark.xfail(reason="Implemented in Phase 4", raises=NotImplementedError, strict=True)
def test_fit_recovers_odmax():
    """Fitted ODmax should be within 5% of the true value on noise-free data."""
    result = fit_growth_model(TIME, OD_CLEAN, model="logistic")
    assert abs(result["ODmax"] - TRUE_ODMAX) / TRUE_ODMAX < 0.05


@pytest.mark.xfail(reason="Implemented in Phase 4", raises=NotImplementedError, strict=True)
def test_fit_recovers_mu_max():
    """Fitted μmax should be within 5% of the true value."""
    result = fit_growth_model(TIME, OD_CLEAN, model="logistic")
    assert abs(result["mu_max"] - TRUE_MU_MAX) / TRUE_MU_MAX < 0.05


@pytest.mark.xfail(reason="Implemented in Phase 4", raises=NotImplementedError, strict=True)
def test_fit_recovers_lag():
    """Fitted lag phase should be within 10% of the true value."""
    result = fit_growth_model(TIME, OD_CLEAN, model="logistic")
    assert abs(result["lag"] - TRUE_LAG) / TRUE_LAG < 0.10


@pytest.mark.xfail(reason="Implemented in Phase 4", raises=NotImplementedError, strict=True)
def test_r_squared_near_one():
    """R² should be > 0.99 on synthetic data with no noise."""
    result = fit_growth_model(TIME, OD_CLEAN, model="logistic")
    assert result["r_squared"] > 0.99
