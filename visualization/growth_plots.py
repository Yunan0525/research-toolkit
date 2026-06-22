# =============================================================================
# visualization/growth_plots.py  –  Growth Curve Plotly Figure Builders
# =============================================================================
# PURPOSE
#   Functions that accept OD600 DataFrames and fitted model parameters and
#   return Plotly Figure objects.
#
# IMPLEMENTED IN PHASE 4
#
# FIGURES PLANNED
#   plot_growth_curves     – multi-group OD vs time with mean ± SD ribbon
#   plot_growth_parameters – bar chart comparing ODmax, mu_max, lag across groups
# =============================================================================

import plotly.graph_objects as go
import pandas as pd
import numpy as np


PALETTE = [
    "#2E8B8B",
    "#E07B6A",
    "#4C72B0",
    "#8DA0CB",
    "#66C2A5",
]


def plot_growth_curves(df: pd.DataFrame, fitted_df: pd.DataFrame = None) -> go.Figure:
    """
    Plot raw OD600 data with optional fitted growth model curves.

    Parameters
    ----------
    df : pd.DataFrame
        Raw OD600 data. Required columns: Time_h, OD600, Group, Replicate.
    fitted_df : pd.DataFrame, optional
        Dense time/OD arrays for fitted curves (one row per time point per group).

    Returns
    -------
    go.Figure

    Raises
    ------
    NotImplementedError
        Until Phase 4 is implemented.
    """
    # ── PHASE 4: implement this function ─────────────────────────────────────
    raise NotImplementedError("plot_growth_curves will be implemented in Phase 4.")


def plot_growth_parameters(summary_df: pd.DataFrame) -> go.Figure:
    """
    Bar chart comparing extracted growth parameters across groups.

    Parameters
    ----------
    summary_df : pd.DataFrame
        Output of analysis.growth_curve.compute_group_statistics().
        Required columns: Group, ODmax, mu_max, lag.

    Returns
    -------
    go.Figure

    Raises
    ------
    NotImplementedError
        Until Phase 4 is implemented.
    """
    # ── PHASE 4: implement this function ─────────────────────────────────────
    raise NotImplementedError("plot_growth_parameters will be implemented in Phase 4.")
