# =============================================================================
# analysis/growth_curve.py  –  Bacterial Growth Curve Analysis Logic
# =============================================================================
# PURPOSE
#   Pure Python functions for fitting growth models and extracting parameters
#   from OD600 time-series data.
#   NO Streamlit code here.
#
# IMPLEMENTED IN PHASE 4
#
# MODELS SUPPORTED (Phase 4)
#   Logistic model:  OD(t) = ODmax / (1 + exp(-r*(t - t_mid)))
#   Gompertz model:  OD(t) = ODmax * exp(-exp(r*exp(1)/ODmax*(lag-t)+1))
#
# KEY PARAMETERS EXTRACTED
#   ODmax   – carrying capacity (maximum optical density)
#   mu_max  – maximum specific growth rate (h⁻¹)
#   lag     – duration of the lag phase (h)
# =============================================================================

import pandas as pd
import numpy as np


def fit_growth_model(
    time: np.ndarray,
    od: np.ndarray,
    model: str = "logistic",
) -> dict:
    """
    Fit a growth model to a single OD600 time series.

    Parameters
    ----------
    time : np.ndarray
        Array of time points (hours).
    od : np.ndarray
        Array of OD600 measurements corresponding to `time`.
    model : str
        Growth model to fit. Options: 'logistic', 'gompertz'.
        Default: 'logistic'.

    Returns
    -------
    dict
        Keys: 'ODmax', 'mu_max', 'lag', 'r_squared', 'model'

    Raises
    ------
    NotImplementedError
        Until Phase 4 is implemented.
    """
    # ── PHASE 4: implement this function ─────────────────────────────────────
    raise NotImplementedError("fit_growth_model will be implemented in Phase 4.")


def compute_group_statistics(
    df: pd.DataFrame,
    time_col: str,
    od_col: str,
    group_col: str,
    replicate_col: str,
    model: str = "logistic",
) -> pd.DataFrame:
    """
    Fit growth models to all groups and replicates; return a summary table.

    Parameters
    ----------
    df : pd.DataFrame
        Input data with time, OD, group, and replicate columns.
    time_col : str
        Column name for time values.
    od_col : str
        Column name for OD600 values.
    group_col : str
        Column name for group labels.
    replicate_col : str
        Column name for replicate labels.
    model : str
        Growth model. Options: 'logistic', 'gompertz'.

    Returns
    -------
    pd.DataFrame
        One row per (group, replicate) with columns:
        Group, Replicate, ODmax, mu_max, lag, r_squared

    Raises
    ------
    NotImplementedError
        Until Phase 4 is implemented.
    """
    # ── PHASE 4: implement this function ─────────────────────────────────────
    raise NotImplementedError("compute_group_statistics will be implemented in Phase 4.")
