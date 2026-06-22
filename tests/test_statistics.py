# =============================================================================
# tests/test_statistics.py  –  Unit Tests for analysis/statistics.py
# =============================================================================
# HOW TO RUN
#   pytest tests/test_statistics.py -v
#
# PHASE 5 STATUS
#   Tests are specifications. xfail until Phase 5.
# =============================================================================

import pytest
import numpy as np
import pandas as pd

from analysis.statistics import (
    run_comparison_test,
    apply_multiple_testing_correction,
    run_correlation_analysis,
)

# ── Synthetic data ────────────────────────────────────────────────────────────
np.random.seed(42)
CONTROL = np.random.normal(loc=1.0, scale=0.2, size=10)
TREATMENT = np.random.normal(loc=2.0, scale=0.2, size=10)  # clearly different

LONG_DF = pd.DataFrame({
    "Group": ["Control"] * 10 + ["Treatment"] * 10,
    "Value": np.concatenate([CONTROL, TREATMENT]),
})


@pytest.mark.xfail(reason="Implemented in Phase 5", raises=NotImplementedError, strict=True)
def test_ttest_detects_difference():
    """t-test should return p < 0.05 for clearly separated groups."""
    result = run_comparison_test(LONG_DF, "Value", "Group", test="t_test")
    p = result["p_value"].values[0]
    assert p < 0.05, f"Expected p < 0.05, got {p}"


@pytest.mark.xfail(reason="Implemented in Phase 5", raises=NotImplementedError, strict=True)
def test_bonferroni_never_increases_pvalues():
    """Bonferroni-corrected p-values must be ≥ raw p-values."""
    raw = pd.Series([0.01, 0.04, 0.20, 0.50])
    corrected = apply_multiple_testing_correction(raw, method="bonferroni")
    assert (corrected >= raw).all()


@pytest.mark.xfail(reason="Implemented in Phase 5", raises=NotImplementedError, strict=True)
def test_pearson_perfect_correlation():
    """Perfectly correlated data should return r ≈ 1.0."""
    x = np.arange(10, dtype=float)
    df = pd.DataFrame({"x": x, "y": x * 2})
    result = run_correlation_analysis(df, "x", "y", method="pearson")
    assert np.isclose(result["r"], 1.0, atol=1e-6)
