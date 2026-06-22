# =============================================================================
# components/upload_widget.py  –  Reusable File Upload Component
# =============================================================================
# PURPOSE
#   A standardised file uploader widget used on every tool page.
#   It combines:
#     1. A Streamlit file_uploader
#     2. A preview of the uploaded data (first 5 rows)
#     3. A "Download sample file" button so users can see the expected format
#
# HOW TO USE
#   from components.upload_widget import render_upload_widget
#
#   uploaded_file = render_upload_widget(
#       tool_name="qPCR Analyzer",
#       sample_file_path="data/sample_qpcr.csv",
#       accepted_types=["csv", "xlsx"],
#   )
#   if uploaded_file is not None:
#       df = load_uploaded_file(uploaded_file)
# =============================================================================

import sys
from pathlib import Path

# ── Path bootstrap ────────────────────────────────────────────────────────────
# components/ is one level below the project root
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import streamlit as st
import pandas as pd
from utils.file_io import load_uploaded_file


def render_upload_widget(
    tool_name: str,
    sample_file_path: str,
    accepted_types: list[str] = None,
    show_preview_rows: int = 5,
):
    """
    Render a file uploader with a data preview and sample file download.

    Parameters
    ----------
    tool_name : str
        Human-readable name of the tool (used in the uploader label).
    sample_file_path : str
        Relative path to the sample dataset file in data/.
    accepted_types : list[str], optional
        List of allowed file extensions without dots, e.g. ['csv', 'xlsx'].
        Defaults to ['csv', 'xlsx'].
    show_preview_rows : int
        Number of rows to show in the data preview table.

    Returns
    -------
    UploadedFile or None
        The raw Streamlit UploadedFile object, or None if nothing uploaded yet.
    """
    if accepted_types is None:
        accepted_types = ["csv", "xlsx"]

    col_upload, col_sample = st.columns([3, 1])

    with col_upload:
        uploaded_file = st.file_uploader(
            label=f"Upload your {tool_name} data",
            type=accepted_types,
            help=(
                f"Accepted formats: {', '.join('.' + t for t in accepted_types)}. "
                f"Maximum file size: 200 MB."
            ),
        )

    with col_sample:
        st.markdown("**Need a template?**")
        _render_sample_download(sample_file_path)

    # ── Data preview ──────────────────────────────────────────────────────────
    if uploaded_file is not None:
        df = load_uploaded_file(uploaded_file)
        if df is not None:
            st.success(
                f"✅ Loaded **{uploaded_file.name}** — "
                f"{len(df):,} rows × {len(df.columns)} columns"
            )
            with st.expander(f"Preview (first {show_preview_rows} rows)", expanded=True):
                st.dataframe(df.head(show_preview_rows), use_container_width=True)

        # Reset the stream position so downstream readers get the full file
        uploaded_file.seek(0)

    return uploaded_file


def _render_sample_download(sample_file_path: str) -> None:
    """
    Show a download button for the sample file if it exists on disk.

    Parameters
    ----------
    sample_file_path : str
        Path relative to the project root.
    """
    path = Path(sample_file_path)
    if path.exists():
        with open(path, "rb") as f:
            file_bytes = f.read()
        st.download_button(
            label=f"⬇ Sample file",
            data=file_bytes,
            file_name=path.name,
            mime="text/csv",
            help="Download this sample file, inspect the format, then upload your own data.",
        )
    else:
        st.caption(f"Sample file not found: `{sample_file_path}`")
