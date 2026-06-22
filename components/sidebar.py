# =============================================================================
# components/sidebar.py  –  Shared Sidebar Component
# =============================================================================
# PURPOSE
#   Renders a consistent sidebar across every page of the application.
#   Called at the top of app.py and every page in pages/.
#
# WHAT IT SHOWS
#   - App logo / name
#   - Version badge
#   - Navigation hint
#   - Links: GitHub, documentation, bug report
#
# HOW TO USE
#   from components.sidebar import render_sidebar
#   render_sidebar()   # call once near the top of each page
# =============================================================================

import streamlit as st


def render_sidebar() -> None:
    """
    Render the shared sidebar content.

    This function is idempotent — calling it multiple times on the same page
    produces the same sidebar without duplicating elements.
    """
    with st.sidebar:
        # ── App identity ──────────────────────────────────────────────────────
        st.markdown("## 🔬 Research Toolkit")
        st.caption("v0.1.0  ·  Biomedical Data Analysis")
        st.divider()

        # ── Navigation guide ──────────────────────────────────────────────────
        st.markdown("### 🧭 Tools")
        st.markdown(
            """
            Use the links above ↑ to navigate between tools.
            Each tool opens on its own page.
            """
        )
        st.divider()

        # ── Resources ─────────────────────────────────────────────────────────
        st.markdown("### 📎 Resources")
        st.markdown(
            """
            - [📖 GitHub Repository](https://github.com/YOUR_USERNAME/research-toolkit)
            - [🐛 Report a Bug](https://github.com/YOUR_USERNAME/research-toolkit/issues)
            - [📋 Input Format Guide](#)
            """
        )
        st.divider()

        # ── Footer caption ────────────────────────────────────────────────────
        st.caption(
            "Built with [Streamlit](https://streamlit.io). "
            "MIT License."
        )
