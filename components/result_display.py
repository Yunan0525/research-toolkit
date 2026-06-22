# =============================================================================
# components/result_display.py  –  Reusable Results Display Component
# =============================================================================
# PURPOSE
#   Renders a consistent layout for showing analysis results:
#   a data table on one side and a Plotly figure on the other.
#   Used by all three tool pages.
#
# HOW TO USE
#   from components.result_display import render_results
#
#   render_results(
#       df=summary_df,
#       fig=fold_change_fig,
#       table_title="Fold Change Summary",
#       figure_title="Fold Change Plot",
#       download_filename_prefix="qpcr_results",
#   )
# =============================================================================

import sys
from pathlib import Path

# ── Path bootstrap ────────────────────────────────────────────────────────────
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from utils.export import (
    download_dataframe_csv,
    download_dataframe_excel,
    download_figure_png,
    download_figure_svg,
)


def render_results(
    df: pd.DataFrame,
    fig: go.Figure,
    table_title: str = "Results Table",
    figure_title: str = "Figure",
    download_filename_prefix: str = "results",
    table_col_width: int = 1,
    figure_col_width: int = 2,
) -> None:
    """
    Display a results table and a Plotly figure side by side with download buttons.

    Parameters
    ----------
    df : pd.DataFrame
        The results table to display.
    fig : go.Figure
        The Plotly figure to display.
    table_title : str
        Section heading for the table column.
    figure_title : str
        Section heading for the figure column.
    download_filename_prefix : str
        Prefix for downloaded files (e.g. "qpcr_results" → "qpcr_results.csv").
    table_col_width : int
        Relative width of the table column.
    figure_col_width : int
        Relative width of the figure column (larger = wider plot).
    """
    col_table, col_fig = st.columns([table_col_width, figure_col_width])

    # ── Table ─────────────────────────────────────────────────────────────────
    with col_table:
        st.subheader(table_title)
        st.dataframe(df, use_container_width=True)

        st.markdown("**Download table:**")
        dl_col1, dl_col2 = st.columns(2)
        with dl_col1:
            download_dataframe_csv(df, filename=f"{download_filename_prefix}.csv")
        with dl_col2:
            download_dataframe_excel(df, filename=f"{download_filename_prefix}.xlsx")

    # ── Figure ────────────────────────────────────────────────────────────────
    with col_fig:
        st.subheader(figure_title)
        st.plotly_chart(fig, use_container_width=True)

        st.markdown("**Download figure:**")
        dl_col3, dl_col4 = st.columns(2)
        with dl_col3:
            download_figure_png(fig, filename=f"{download_filename_prefix}.png")
        with dl_col4:
            download_figure_svg(fig, filename=f"{download_filename_prefix}.svg")
