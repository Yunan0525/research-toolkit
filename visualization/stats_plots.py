# =============================================================================
# visualization/stats_plots.py  –  Statistics Plotly Figure Builders
# =============================================================================
# PURPOSE
#   Functions that build Plotly figures for the Statistical Analysis module.
#
# IMPLEMENTED IN PHASE 5
#
# FIGURES PLANNED
#   plot_box_violin        – box + individual points, significance brackets
#   plot_correlation       – scatter plot with regression line and R²
#   plot_pvalue_heatmap    – pairwise p-value matrix heatmap
# =============================================================================

import plotly.graph_objects as go
import pandas as pd


PALETTE = [
    "#2E8B8B",
    "#E07B6A",
    "#4C72B0",
    "#8DA0CB",
    "#66C2A5",
]


def plot_box_violin(df: pd.DataFrame, value_col: str, group_col: str) -> go.Figure:
    """
    Box plot + violin + individual data points for group comparisons.

    Parameters
    ----------
    df : pd.DataFrame
        Input data in long format.
    value_col : str
        Column with numeric values.
    group_col : str
        Column with group labels.

    Returns
    -------
    go.Figure

    Raises
    ------
    NotImplementedError
        Until Phase 5 is implemented.
    """
    # ── PHASE 5: implement this function ─────────────────────────────────────
    raise NotImplementedError("plot_box_violin will be implemented in Phase 5.")


def plot_correlation(df: pd.DataFrame, x_col: str, y_col: str) -> go.Figure:
    """
    Scatter plot with Pearson/Spearman regression line and confidence band.

    Raises
    ------
    NotImplementedError
        Until Phase 5 is implemented.
    """
    # ── PHASE 5: implement this function ─────────────────────────────────────
    raise NotImplementedError("plot_correlation will be implemented in Phase 5.")
