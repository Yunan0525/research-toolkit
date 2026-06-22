# =============================================================================
# app.py  –  Research Data Analysis Toolkit
# =============================================================================
# PURPOSE
#   Entry point for the Streamlit application and the homepage.
#   Streamlit Cloud runs this file first; all other pages live in pages/.
#
# STRUCTURE
#   1. Page config (must be the very first Streamlit call)
#   2. Shared sidebar
#   3. Hero header and description
#   4. Tool cards grid
#   5. Getting-started instructions
#   6. Footer
#
# ADDING A NEW TOOL
#   1. Add a new entry to the TOOLS list below.
#   2. Create the corresponding page file in pages/.
#   3. Create the analysis module in analysis/.
#   4. Create the visualisation module in visualization/.
# =============================================================================

import streamlit as st

# ── Page configuration ────────────────────────────────────────────────────────
# This MUST be the first Streamlit call in the script.
st.set_page_config(
    page_title="Research Data Analysis Toolkit",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Shared components ─────────────────────────────────────────────────────────
from components.sidebar import render_sidebar

render_sidebar()

# ── Hero ──────────────────────────────────────────────────────────────────────
st.title("🔬 Research Data Analysis Toolkit")
st.markdown(
    """
    A modular, open-source platform for **biomedical and microbiome data analysis**.
    Select a tool from the sidebar, upload your data, and download publication-quality
    results — no coding required.
    """
)

st.divider()

# ── Tool registry ─────────────────────────────────────────────────────────────
# To add a new tool:
#   1. Append a dict to this list.
#   2. Set "status" to "coming_soon" until the tool is implemented.
TOOLS = [
    {
        "icon": "🧬",
        "title": "qPCR ΔΔCt Analyzer",
        "description": (
            "Calculate ΔCt, ΔΔCt, Fold Change, and log₂ Fold Change from "
            "your qPCR Ct values. Supports multiple reference genes and "
            "generates annotated bar plots ready for publication."
        ),
        "status": "available",
    },
    {
        "icon": "📈",
        "title": "Growth Curve Analyzer",
        "description": (
            "Fit logistic or Gompertz growth models to OD₆₀₀ time-series data. "
            "Extracts maximum OD, growth rate (μmax), and lag phase. Compare "
            "treatment groups on a single interactive plot."
        ),
        "status": "available",
    },
    {
        "icon": "📊",
        "title": "Statistical Analysis Module",
        "description": (
            "Run t-tests, paired t-tests, Wilcoxon, ANOVA, Kruskal-Wallis, and "
            "correlation analyses. Applies Bonferroni or Benjamini-Hochberg "
            "correction for multiple comparisons."
        ),
        "status": "available",
    },
]

# ── Tool cards ────────────────────────────────────────────────────────────────
cols = st.columns(len(TOOLS), gap="large")

for col, tool in zip(cols, TOOLS):
    with col:
        with st.container(border=True):
            st.markdown(f"## {tool['icon']} {tool['title']}")
            st.markdown(tool["description"])

            if tool["status"] == "available":
                st.success("✅ Available — open from the sidebar")
            else:
                st.info("🚧 Coming soon")

st.divider()

# ── Getting-started guide ─────────────────────────────────────────────────────
st.subheader("Getting started")

with st.expander("How to use this toolkit", expanded=True):
    st.markdown(
        """
        1. **Choose a tool** using the navigation links in the left sidebar.
        2. **Upload your data file** (CSV or Excel). Each tool shows a sample
           file you can download to understand the expected format.
        3. **Configure the analysis** with the dropdowns and controls that
           appear after your file loads.
        4. **Review results** in the interactive tables and plots.
        5. **Download** your results as CSV (tables) or PNG/SVG (figures).
        """
    )

with st.expander("Supported file formats"):
    st.markdown(
        """
        | Format | Extension | Notes |
        |--------|-----------|-------|
        | Comma-separated values | `.csv` | UTF-8 encoding preferred |
        | Excel workbook | `.xlsx` | First sheet used by default |
        """
    )

with st.expander("Citing this toolkit"):
    st.markdown(
        """
        Each tool page shows the exact parameters used in your analysis.
        Copy those into your Methods section for full reproducibility.

        Please cite:
        > *Research Data Analysis Toolkit* (v0.1.0).
        > Available at: https://github.com/YOUR_USERNAME/research-toolkit
        """
    )

# ── Footer ────────────────────────────────────────────────────────────────────
st.divider()
st.caption(
    "Research Data Analysis Toolkit · v0.1.0 · "
    "Built with [Streamlit](https://streamlit.io) · "
    "[View source on GitHub](https://github.com/YOUR_USERNAME/research-toolkit)"
)
