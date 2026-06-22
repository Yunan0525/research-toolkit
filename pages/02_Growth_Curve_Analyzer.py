# =============================================================================
# pages/02_Growth_Curve_Analyzer.py  –  Growth Curve Analyzer
# =============================================================================
# PURPOSE
#   Streamlit page for OD600 growth curve analysis.
#   This file handles ONLY the user interface.
#
# ANALYSIS LOGIC lives in:  analysis/growth_curve.py
# PLOT FUNCTIONS live in:   visualization/growth_plots.py
#
# PHASE 4 will fully implement this page.
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
    page_title="Growth Curve Analyzer · Research Toolkit",
    page_icon="📈",
    layout="wide",
)

render_sidebar()

# ── Header ────────────────────────────────────────────────────────────────────
st.title("📈 Growth Curve Analyzer")
st.markdown(
    """
    Upload OD₆₀₀ time-series measurements to visualise bacterial growth curves
    and extract key parameters: **maximum OD**, **growth rate (μmax)**,
    and **lag phase duration**.
    """
)

st.info(
    "🚧 **Phase 4 — Coming soon.** "
    "Full implementation arrives in Phase 4.",
    icon="ℹ️",
)

st.divider()

# ── Expected format ───────────────────────────────────────────────────────────
st.subheader("Expected input format")
st.markdown(
    """
    | Column | Description | Example values |
    |--------|-------------|----------------|
    | Time_h | Time in hours | `0`, `1`, `2`, `12` |
    | OD600 | Optical density at 600 nm | `0.05`, `0.12`, `0.45` |
    | Replicate | Biological replicate label | `R1`, `R2`, `R3` |
    | Group | Treatment / strain label | `WT`, `Mutant_A` |
    """
)

# ── Placeholder uploader ──────────────────────────────────────────────────────
st.subheader("Upload your data")
uploaded_file = render_upload_widget(
    tool_name="Growth Curve Analyzer",
    sample_file_path="data/sample_growth_curve.csv",
    accepted_types=["csv", "xlsx"],
)

if uploaded_file is not None:
    st.warning(
        "File received! Full analysis will be available in Phase 4.",
        icon="⏳",
    )
