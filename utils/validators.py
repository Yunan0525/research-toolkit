# =============================================================================
# utils/validators.py  –  Input Validation and User-Friendly Error Messages
# =============================================================================
# PURPOSE
#   Validate DataFrames and user selections BEFORE running any analysis.
#   Each function returns (is_valid: bool, message: str).
#   The calling page decides how to display the message (st.error / st.warning).
#
# DESIGN PRINCIPLE
#   Never let a raw Python traceback reach the user.
#   Every validation failure should tell the user:
#     1. What went wrong.
#     2. How to fix it.
# =============================================================================

import pandas as pd
import numpy as np
from typing import Tuple


def check_required_columns(
    df: pd.DataFrame,
    required: list[str],
) -> Tuple[bool, str]:
    """
    Verify that all required column names exist in the DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    required : list[str]
        Column names that must be present.

    Returns
    -------
    (True, "") if all columns present.
    (False, error_message) if any column is missing.
    """
    missing = [col for col in required if col not in df.columns]
    if missing:
        return False, (
            f"❌ Missing columns: **{', '.join(missing)}**. "
            f"Your file has: {', '.join(df.columns.tolist())}. "
            "Please check the expected format above and re-upload your file."
        )
    return True, ""


def check_numeric_column(
    df: pd.DataFrame,
    column: str,
) -> Tuple[bool, str]:
    """
    Verify that a column contains numeric (non-string) values.

    Parameters
    ----------
    df : pd.DataFrame
    column : str
        Column to check.

    Returns
    -------
    (True, "") if numeric.
    (False, error_message) if non-numeric values found.
    """
    non_numeric = pd.to_numeric(df[column], errors="coerce").isna().sum()
    total = len(df[column].dropna())

    if non_numeric > 0:
        return False, (
            f"❌ Column **{column}** contains {non_numeric} non-numeric value(s) "
            f"out of {total} rows. "
            "Make sure Ct values are numbers (e.g. 24.5) with no text, "
            "units, or extra spaces."
        )
    return True, ""


def check_minimum_rows(
    df: pd.DataFrame,
    minimum: int = 2,
) -> Tuple[bool, str]:
    """
    Verify that the DataFrame has enough rows to perform analysis.

    Returns
    -------
    (True, "") if row count ≥ minimum.
    (False, error_message) otherwise.
    """
    if len(df) < minimum:
        return False, (
            f"❌ Your file contains only {len(df)} data row(s), "
            f"but at least {minimum} rows are required for analysis."
        )
    return True, ""


def check_no_all_nan_columns(df: pd.DataFrame) -> Tuple[bool, str]:
    """
    Warn if any column is entirely empty (all NaN).

    Returns
    -------
    (True, "") if no fully-empty columns.
    (False, warning_message) if any column is all NaN.
    """
    empty_cols = [col for col in df.columns if df[col].isna().all()]
    if empty_cols:
        return False, (
            f"⚠️ The following column(s) are completely empty: "
            f"**{', '.join(empty_cols)}**. "
            "They will be ignored during analysis."
        )
    return True, ""


def check_group_has_replicates(
    df: pd.DataFrame,
    group_col: str,
    min_replicates: int = 2,
) -> Tuple[bool, str]:
    """
    Verify that every group has at least `min_replicates` rows.

    Parameters
    ----------
    df : pd.DataFrame
    group_col : str
        Column containing group labels.
    min_replicates : int
        Minimum number of rows required per group.

    Returns
    -------
    (True, "") if all groups have enough replicates.
    (False, warning_message) listing under-represented groups.
    """
    counts = df[group_col].value_counts()
    small_groups = counts[counts < min_replicates].index.tolist()

    if small_groups:
        return False, (
            f"⚠️ The following group(s) have fewer than {min_replicates} "
            f"replicate(s): **{', '.join(str(g) for g in small_groups)}**. "
            "Statistical results for these groups may be unreliable."
        )
    return True, ""
