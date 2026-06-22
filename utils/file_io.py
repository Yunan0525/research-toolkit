# =============================================================================
# utils/file_io.py  –  File Upload and Parsing Utilities
# =============================================================================
# PURPOSE
#   Handles reading CSV and Excel files uploaded via Streamlit's file_uploader.
#   Returns clean pandas DataFrames.
#   Called by every tool page before any analysis runs.
#
# WHY A SEPARATE FILE?
#   All three tools accept uploads in the same way. Centralising this logic
#   means a bug fix or format change only needs to be made in one place.
# =============================================================================

import pandas as pd
import streamlit as st
from io import BytesIO


def load_uploaded_file(uploaded_file) -> pd.DataFrame | None:
    """
    Read a Streamlit UploadedFile object into a pandas DataFrame.

    Supports:
        - .csv  (comma or semicolon delimited, auto-detected)
        - .xlsx (first sheet)
        - .xls  (first sheet, legacy Excel)

    Parameters
    ----------
    uploaded_file : streamlit.runtime.uploaded_file_manager.UploadedFile
        The object returned by st.file_uploader().
        Pass None if no file has been uploaded yet.

    Returns
    -------
    pd.DataFrame
        Parsed data, or None if uploaded_file is None.

    Notes
    -----
    - Trailing whitespace is stripped from all string values.
    - Column names are stripped and whitespace-normalised.
    - On parse error, a user-friendly st.error() message is shown
      and None is returned (the calling page should check for None).
    """
    if uploaded_file is None:
        return None

    file_name = uploaded_file.name.lower()
    raw_bytes = BytesIO(uploaded_file.read())

    try:
        if file_name.endswith(".csv"):
            # Try comma delimiter first; fall back to semicolon
            df = pd.read_csv(raw_bytes)
            if df.shape[1] == 1:
                raw_bytes.seek(0)
                df = pd.read_csv(raw_bytes, sep=";")

        elif file_name.endswith((".xlsx", ".xls")):
            df = pd.read_excel(raw_bytes, sheet_name=0)

        else:
            st.error(
                f"❌ Unsupported file type: **{uploaded_file.name}**. "
                "Please upload a `.csv` or `.xlsx` file."
            )
            return None

    except Exception as exc:
        st.error(
            f"❌ Could not read **{uploaded_file.name}**. "
            f"Make sure the file is not corrupted or password-protected.\n\n"
            f"Technical detail: `{exc}`"
        )
        return None

    # ── Clean up column names ─────────────────────────────────────────────────
    df.columns = df.columns.str.strip()

    # ── Strip leading/trailing whitespace from string cells ───────────────────
    str_cols = df.select_dtypes(include="object").columns
    df[str_cols] = df[str_cols].apply(lambda col: col.str.strip())

    return df


def load_sample_file(file_path: str) -> pd.DataFrame | None:
    """
    Load one of the bundled sample datasets from the data/ directory.

    Parameters
    ----------
    file_path : str
        Relative path from the project root, e.g. 'data/sample_qpcr.csv'.

    Returns
    -------
    pd.DataFrame or None
    """
    try:
        return pd.read_csv(file_path)
    except FileNotFoundError:
        st.warning(f"Sample file not found at `{file_path}`.")
        return None
    except Exception as exc:
        st.warning(f"Could not load sample file: `{exc}`")
        return None
