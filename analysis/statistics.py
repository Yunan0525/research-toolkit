# =============================================================================
# analysis/statistics.py  –  Statistical Analysis Logic
# =============================================================================
# PURPOSE
#   Pure Python functions for hypothesis testing, effect sizes, and
#   multiple-testing correction.
#   NO Streamlit code here.
#
# IMPLEMENTED IN PHASE 5
#
# TESTS SUPPORTED (Phase 5)
#   Parametric:    Student's t-test, paired t-test, one-way ANOVA
#   Non-parametric: Wilcoxon signed-rank, Kruskal-Wallis
#   Correlation:   Pearson, Spearman
#
# CORRECTIONS SUPPORTED
#   Bonferroni, Benjamini-Hochberg (FDR)
# =============================================================================

import pandas as pd
import numpy as np


def run_comparison_test(
    df: pd.DataFrame,
    value_col: str,
    group_col: str,
    test: str,
    subject_col: str = None,
) -> pd.DataFrame:
    """
    Run a pairwise or multi-group comparison test.

    Parameters
    ----------
    df : pd.DataFrame
        Input data in long format.
    value_col : str
        Column containing numeric measurements.
    group_col : str
        Column containing group labels.
    test : str
        Statistical test to run. One of:
        't_test', 'paired_t_test', 'wilcoxon', 'anova', 'kruskal'.
    subject_col : str, optional
        Column with subject identifiers (required for paired_t_test).

    Returns
    -------
    pd.DataFrame
        Results table with columns: Group1, Group2, Statistic, p_value,
        effect_size, test_used

    Raises
    ------
    NotImplementedError
        Until Phase 5 is implemented.
    """
    # ── PHASE 5: implement this function ─────────────────────────────────────
    raise NotImplementedError("run_comparison_test will be implemented in Phase 5.")


def apply_multiple_testing_correction(
    p_values: pd.Series,
    method: str = "fdr_bh",
) -> pd.Series:
    """
    Apply multiple-testing correction to a series of p-values.

    Parameters
    ----------
    p_values : pd.Series
        Raw p-values from statistical tests.
    method : str
        Correction method. Options: 'bonferroni', 'fdr_bh' (Benjamini-Hochberg).

    Returns
    -------
    pd.Series
        Adjusted p-values in the same order as the input.

    Raises
    ------
    NotImplementedError
        Until Phase 5 is implemented.
    """
    # ── PHASE 5: implement this function ─────────────────────────────────────
    raise NotImplementedError(
        "apply_multiple_testing_correction will be implemented in Phase 5."
    )


def run_correlation_analysis(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    method: str = "pearson",
) -> dict:
    """
    Compute correlation between two numeric variables.

    Parameters
    ----------
    df : pd.DataFrame
        Input data.
    x_col : str
        Column name for the first variable.
    y_col : str
        Column name for the second variable.
    method : str
        Correlation method. Options: 'pearson', 'spearman'.

    Returns
    -------
    dict
        Keys: 'r', 'p_value', 'ci_lower', 'ci_upper', 'n', 'method'

    Raises
    ------
    NotImplementedError
        Until Phase 5 is implemented.
    """
    # ── PHASE 5: implement this function ─────────────────────────────────────
    raise NotImplementedError(
        "run_correlation_analysis will be implemented in Phase 5."
    )
