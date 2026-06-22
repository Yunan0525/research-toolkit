# =============================================================================
# pages/03_Statistical_Analysis.py  –  Statistical Analysis Module
# =============================================================================
# PURPOSE
#   Streamlit page for statistical hypothesis testing and correlation analysis.
#   This file handles ONLY the user interface.
#
# ANALYSIS LOGIC lives in:  analysis/statistics.py
# PLOT FUNCTIONS live in:   visualization/stats_plots.py
#
# PHASE 5 will fully implement this page.
# =============================================================================

import sys
from pathlib import Path

# ── Path bootstrap ────────────────────────────────────────────────────────────
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import streamlit as st
from components.sidebar import render_sidebar
from components.upload_widget import render_upload_widget

st.set_page_config(
    page_title="Statistical Analysis · Research Toolkit",
    page_icon="📊",
    layout="wide",
)

render_sidebar()

# ── Header ────────────────────────────────────────────────────────────────────
st.title("📊 Statistical Analysis Module")
st.markdown(
    """
    Upload grouped numeric data to run hypothesis tests, calculate effect sizes,
    and apply multiple-testing corrections.

    **Available tests:** t-test · paired t-test · Wilcoxon · one-way ANOVA ·
    Kruskal-Wallis · Pearson / Spearman correlation

    **Corrections:** Bonferroni · Benjamini-Hochberg (FDR)
    """
)

st.info(
    "🚧 **Phase 5 — Coming soon.** "
    "Full implementation arrives in Phase 5.",
    icon="ℹ️",
)

st.divider()

# ── Expected format ───────────────────────────────────────────────────────────
st.subheader("Expected input format")
st.markdown(
    """
    | Column | Description | Example values |
    |--------|-------------|----------------|
    | Group | Group label for comparison | `Control`, `Treatment_A` |
    | Value | Numeric measurement | `1.23`, `4.56` |
    | Subject | Subject ID for paired tests (optional) | `P01`, `P02` |
    """
)

# ── Placeholder uploader ──────────────────────────────────────────────────────
st.subheader("Upload your data")
uploaded_file = render_upload_widget(
    tool_name="Statistical Analysis",
    sample_file_path="data/sample_stats_data.csv",
    accepted_types=["csv", "xlsx"],
)

if uploaded_file is not None:
    st.warning(
        "File received! Full analysis will be available in Phase 5.",
        icon="⏳",
    )
