# =============================================================================
# pages/01_qPCR_Analyzer.py  –  qPCR ΔΔCt Analyzer  (corrected pipeline)
# =============================================================================
# PURPOSE
#   Streamlit UI for the qPCR ΔΔCt analysis pipeline with technical replicate
#   handling and statistical testing.
#
# MANDATORY PIPELINE ORDER
#   1. average_technical_replicates()  – ALWAYS first; collapses tech reps
#   2. calculate_delta_ct()            – ΔCt at biological-sample level
#   3. calculate_delta_delta_ct()      – ΔΔCt vs control group mean ΔCt
#   4. summarise_results()             – mean ± SD across biological replicates
#   5. run_statistical_tests()         – tests on ΔCt biological replicates
#
# All mathematics live in analysis/qpcr.py.
# All figure construction lives in visualization/qpcr_plots.py.
# =============================================================================

import sys
from pathlib import Path

# ── Path bootstrap ────────────────────────────────────────────────────────────
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import streamlit as st
import pandas as pd
import numpy as np

from components.sidebar import render_sidebar
from utils.file_io import load_uploaded_file, load_sample_file
from utils.export import (
    download_dataframe_csv,
    download_dataframe_excel,
    download_figure_png,
    download_figure_png_600dpi,
    download_figure_svg,
)
from analysis.qpcr import (
    average_technical_replicates,
    calculate_delta_ct,
    calculate_delta_delta_ct,
    summarise_results,
    run_statistical_tests,
)
from visualization.qpcr_plots import (
    plot_fold_change_bar,
    plot_log2fc_bar,
    plot_deltact_heatmap,
    plot_stats_results,
)

# =============================================================================
# Page setup
# =============================================================================
st.set_page_config(
    page_title="qPCR Analyzer · Research Toolkit",
    page_icon="🧬",
    layout="wide",
)
render_sidebar()

st.title("🧬 qPCR ΔΔCt Analyzer")
st.markdown(
    """
    Handles **technical replicates** correctly: multiple Ct rows for the same
    sample × gene are averaged first, then ΔCt, ΔΔCt, and Fold Change are
    calculated at the **biological replicate** level.
    Statistical tests compare biological replicate ΔCt values between groups.
    """
)
st.divider()


def _stop(msg: str) -> None:
    """Show an error message and halt page rendering."""
    st.error(msg)
    st.stop()


def _guess(candidates: list, columns: list) -> str | None:
    """Return the first column name that case-insensitively matches a candidate."""
    for c in candidates:
        for col in columns:
            if c.lower() == col.lower():
                return col
    return None


# =============================================================================
# STEP 1 — Load data
# =============================================================================
st.subheader("Step 1 — Load your data")

input_mode = st.radio(
    "Data source",
    ["Upload my own file", "Use the built-in sample dataset"],
    horizontal=True,
    help=(
        "The sample dataset contains technical triplicates: each biological "
        "sample × gene has 3 Ct rows, simulating real qPCR plate data."
    ),
)

raw_df: pd.DataFrame | None = None

if input_mode == "Upload my own file":
    up_col, hint_col = st.columns([3, 1])
    with up_col:
        uploaded = st.file_uploader(
            "Upload CSV or Excel",
            type=["csv", "xlsx", "xls"],
            help="Max 200 MB. First Excel sheet used.",
        )
    with hint_col:
        p = _PROJECT_ROOT / "data" / "sample_qpcr.csv"
        if p.exists():
            st.markdown("**Format template:**")
            st.download_button(
                "⬇ Sample CSV",
                data=p.read_bytes(),
                file_name="sample_qpcr.csv",
                mime="text/csv",
            )
    if uploaded:
        raw_df = load_uploaded_file(uploaded)
else:
    raw_df = load_sample_file(str(_PROJECT_ROOT / "data" / "sample_qpcr.csv"))
    if raw_df is not None:
        st.success(
            "✅ Sample dataset loaded — 2 groups (Control, Treatment), "
            "3 biological replicates each, **3 technical replicates per well**, "
            "3 target genes (IL6, TNF, CXCL10), reference gene GAPDH."
        )

if raw_df is None:
    st.info(
        "👆 Upload a file or select the sample dataset to begin.",
        icon="ℹ️",
    )
    st.stop()

with st.expander(
    f"📋 Raw data preview — {len(raw_df):,} rows × {raw_df.shape[1]} columns "
    "(includes technical replicates)",
    expanded=True,
):
    st.dataframe(raw_df.head(24), use_container_width=True)

st.divider()

# =============================================================================
# STEP 2 — Map columns
# =============================================================================
st.subheader("Step 2 — Map your columns")
st.markdown(
    "Select which column in your file corresponds to each field. "
    "Column names do **not** need to match exactly."
)

cols = raw_df.columns.tolist()
c1, c2, c3, c4 = st.columns(4)

with c1:
    sample_col = st.selectbox(
        "Sample column",
        cols,
        index=cols.index(_guess(["Sample", "sample_id", "SampleID", "ID"], cols) or cols[0]),
        help=(
            "Biological sample identifier (e.g. WT_R1). "
            "Multiple rows with the same Sample + Gene = technical replicates."
        ),
    )
with c2:
    gene_col = st.selectbox(
        "Gene column",
        cols,
        index=cols.index(_guess(["Gene", "gene_name", "Target", "Assay"], cols) or cols[0]),
        help="Gene names — include both target genes and the reference gene.",
    )
with c3:
    ct_col = st.selectbox(
        "Ct column",
        cols,
        index=cols.index(_guess(["Ct", "CT", "ct_value", "Cq", "Cp"], cols) or cols[0]),
        help="Raw Ct / Cq values. Non-numeric values (e.g. 'Undetermined') will be excluded.",
    )
with c4:
    group_col = st.selectbox(
        "Group / Treatment column",
        cols,
        index=cols.index(_guess(["Group", "Treatment", "Condition", "Genotype"], cols) or cols[0]),
        help="Experimental group (e.g. Control, Treatment).",
    )

# Duplicate column guard
if len({sample_col, gene_col, ct_col, group_col}) < 4:
    st.warning("⚠️ Two or more fields map to the same column. Please select distinct columns.")

# Numeric Ct check
n_bad = pd.to_numeric(raw_df[ct_col], errors="coerce").isna().sum() - raw_df[ct_col].isna().sum()
if n_bad > 0:
    st.warning(
        f"⚠️ **{ct_col}** contains {n_bad} non-numeric value(s) "
        "(e.g. 'Undetermined'). These rows will be excluded from averaging."
    )

# Show tech-replicate structure
with st.expander("🔍 Technical replicate structure detected in your data"):
    rep_counts = (
        raw_df.groupby([sample_col, gene_col])
        .size()
        .reset_index(name="N_tech_reps")
    )
    min_reps = int(rep_counts["N_tech_reps"].min())
    max_reps = int(rep_counts["N_tech_reps"].max())
    n_bio = raw_df[sample_col].nunique()
    n_genes = raw_df[gene_col].nunique()

    info_cols = st.columns(4)
    info_cols[0].metric("Biological samples", n_bio)
    info_cols[1].metric("Unique genes", n_genes)
    info_cols[2].metric("Min tech reps / well", min_reps)
    info_cols[3].metric("Max tech reps / well", max_reps)

    if min_reps != max_reps:
        st.warning(
            "⚠️ Unequal numbers of technical replicates detected. "
            "Means are still computed correctly; wells with fewer reps have less weight."
        )
    st.dataframe(rep_counts, use_container_width=True, hide_index=True)

st.divider()

# =============================================================================
# STEP 3 — Analysis parameters
# =============================================================================
st.subheader("Step 3 — Set analysis parameters")

available_genes  = sorted(raw_df[gene_col].dropna().unique().tolist())
available_groups = raw_df[group_col].dropna().unique().tolist()

p_col, c_col, s_col, corr_col = st.columns(4)

with p_col:
    reference_gene = st.selectbox(
        "Reference (housekeeping) gene",
        available_genes,
        index=available_genes.index(
            _guess(["GAPDH", "ACTB", "B2M", "HPRT1", "18S", "RPL13A"], available_genes)
            or available_genes[0]
        ) if _guess(["GAPDH", "ACTB", "B2M", "HPRT1", "18S", "RPL13A"], available_genes) else 0,
        help="Used to normalise all target genes. ΔCt = Ct(target) − Ct(this gene).",
    )

with c_col:
    control_group = st.selectbox(
        "Control group",
        available_groups,
        help="ΔΔCt = ΔCt(sample) − mean ΔCt of this group's biological replicates.",
    )

with s_col:
    stat_test = st.selectbox(
        "Statistical test",
        ["auto", "welch", "t-test", "wilcoxon", "anova", "kruskal"],
        index=0,
        help=(
            "**auto**: Welch's t-test for 2 groups, one-way ANOVA for ≥3 groups. "
            "Tests are run on ΔCt biological replicate values (NOT on Fold Change)."
        ),
        format_func=lambda x: {
            "auto":     "Auto (Welch / ANOVA)",
            "welch":    "Welch's t-test (unequal variance)",
            "t-test":   "Student's t-test (equal variance)",
            "wilcoxon": "Wilcoxon rank-sum (non-parametric)",
            "anova":    "One-way ANOVA",
            "kruskal":  "Kruskal-Wallis (non-parametric)",
        }[x],
    )

with corr_col:
    correction = st.selectbox(
        "Multiple testing correction",
        ["fdr_bh", "bonferroni", "none"],
        index=0,
        help=(
            "Applied across all tested genes. "
            "**fdr_bh** = Benjamini-Hochberg (recommended); "
            "**bonferroni** = conservative; "
            "**none** = raw p-values only."
        ),
        format_func=lambda x: {
            "fdr_bh":    "Benjamini-Hochberg (FDR)",
            "bonferroni":"Bonferroni",
            "none":      "None (raw p-values)",
        }[x],
    )

target_genes = [g for g in available_genes if g != reference_gene]
if not target_genes:
    _stop(
        f"❌ No target genes found after excluding '{reference_gene}'. "
        "Your data must contain at least one gene other than the reference."
    )

with st.expander("Genes that will be analysed"):
    gc1, gc2 = st.columns(2)
    gc1.markdown(f"**Reference gene:** {reference_gene}")
    gc2.markdown(f"**Target genes ({len(target_genes)}):** {', '.join(target_genes)}")

st.divider()

# =============================================================================
# STEP 4 — Run
# =============================================================================
st.subheader("Step 4 — Run the analysis")

if not st.button("▶ Run ΔΔCt Analysis", type="primary"):
    st.info("Configure the settings above, then click **▶ Run ΔΔCt Analysis**.", icon="ℹ️")
    st.stop()

with st.spinner("Running analysis…"):
    try:
        # ── 1. Average technical replicates ───────────────────────────────────
        avg_df = average_technical_replicates(
            raw_df, sample_col, gene_col, ct_col, group_col
        )

        # ── 2. ΔCt at biological-sample level ────────────────────────────────
        delta_ct_df = calculate_delta_ct(
            avg_df, sample_col, gene_col, ct_col, reference_gene
        )

        # ── 3. ΔΔCt + Fold Change ────────────────────────────────────────────
        detail_df = calculate_delta_delta_ct(
            delta_ct_df, group_col, control_group
        )

        # ── 4. Biological-replicate summary (mean ± SD) ───────────────────────
        summary_df = summarise_results(detail_df, group_col)

        # ── 5. Statistics on ΔCt biological replicates ───────────────────────
        stats_df = run_statistical_tests(
            delta_ct_df, group_col,
            test=stat_test,
            correction_method=correction,
        )

    except ValueError as e:
        _stop(f"❌ Analysis error: {e}")
    except Exception as e:
        _stop(
            f"❌ Unexpected error: {e}\n\n"
            "Check your column selections and data format."
        )

st.success("✅ Analysis complete!", icon="✅")

# Quick averaging summary
n_raw_rows   = len(raw_df)
n_avg_rows   = len(avg_df)
avg_tech_n   = avg_df["Ct_Tech_N"].mean() if "Ct_Tech_N" in avg_df.columns else "?"
st.info(
    f"ℹ️ **Technical replicate averaging:** {n_raw_rows} raw rows → "
    f"{n_avg_rows} biological sample × gene entries "
    f"(average {avg_tech_n:.1f} tech reps per well).",
    icon="ℹ️",
)

st.divider()

# =============================================================================
# GROUP SELECTION  –  filter which groups appear in plots and summary tables
# =============================================================================
st.subheader("Step 5 — Select groups to display")

all_result_groups = detail_df[group_col].unique().tolist()

selected_groups = st.multiselect(
    "Groups to include in plots and summary tables",
    options=all_result_groups,
    default=all_result_groups,
    help=(
        "Choose which experimental groups appear in the bar charts and "
        "per-replicate dot overlay. All groups are always included in the "
        "raw data tables (tabs 2 and 3). "
        "Select at least 1 group to generate plots."
    ),
)

if len(selected_groups) == 0:
    st.warning(
        "⚠️ No groups selected. Please select at least one group to display plots.",
        icon="⚠️",
    )
    st.stop()

if len(selected_groups) < len(all_result_groups):
    st.info(
        f"ℹ️ Showing {len(selected_groups)} of {len(all_result_groups)} groups: "
        f"{', '.join(selected_groups)}.",
        icon="ℹ️",
    )

# Filter all result DataFrames to selected groups
plot_summary_df = summary_df[summary_df["Group"].isin(selected_groups)].copy()
plot_detail_df  = detail_df[detail_df[group_col].isin(selected_groups)].copy()
plot_stats_df   = (
    stats_df[
        stats_df["Groups_Compared"].apply(
            lambda s: any(g in str(s) for g in selected_groups)
        )
    ].copy()
    if stats_df is not None and not stats_df.empty and "Groups_Compared" in stats_df.columns
    else stats_df
)

st.divider()

# =============================================================================
# STEP 6 — Results
# =============================================================================
st.subheader("Step 6 — Results")

(
    tab_summary, tab_avg, tab_detail,
    tab_fc, tab_log2, tab_heatmap, tab_stats,
) = st.tabs([
    "📋 Biological Summary",
    "🔬 Averaged Ct Table",
    "📄 Per-Replicate Detail",
    "📊 Fold Change Plot",
    "📊 log₂FC Plot",
    "🟥 ΔCt Heatmap",
    "📈 Statistics",
])

# ── Tab 1: Biological replicate summary ──────────────────────────────────────
with tab_summary:
    st.markdown(
        "**Mean ± SD** across biological replicates (technical replicates already averaged). "
        "SD quantifies **biological variability** between independent experiments."
    )

    display_sum = summary_df.rename(columns={
        "N_BioReps":             "N (bio reps)",
        "Mean_FC":               "Mean FC",
        "SD_FC":                 "SD FC",
        "SEM_FC":                "SEM FC",
        "Mean_log2FC":           "Mean log₂FC",
        "SD_log2FC":             "SD log₂FC",
        "SEM_log2FC":            "SEM log₂FC",
        "Mean_Delta_Ct":         "Mean ΔCt",
        "SD_Delta_Ct":           "SD ΔCt",
        "SEM_Delta_Ct":          "SEM ΔCt",
        "Mean_Delta_Delta_Ct":   "Mean ΔΔCt",
        "SD_Delta_Delta_Ct":     "SD ΔΔCt",
    })
    st.dataframe(display_sum, use_container_width=True, hide_index=True)

    d1, d2 = st.columns(2)
    with d1:
        download_dataframe_csv(summary_df, "qpcr_summary.csv", "⬇ Download CSV")
    with d2:
        download_dataframe_excel(summary_df, "qpcr_summary.xlsx", "Summary", "⬇ Download Excel")

# ── Tab 2: Averaged Ct table (after tech-rep collapse) ───────────────────────
with tab_avg:
    st.markdown(
        "Mean Ct per **biological sample × gene** after averaging technical replicates. "
        "**Ct_Tech_SD > 0.5** may indicate a problematic well — check your raw data."
    )

    display_avg = avg_df.round(4)
    # Highlight high tech-rep SD
    def _highlight_sd(val):
        if isinstance(val, float) and val > 0.5:
            return "background-color: #FFF3CD"
        return ""

    if "Ct_Tech_SD" in display_avg.columns:
        st.dataframe(
            display_avg.style.applymap(_highlight_sd, subset=["Ct_Tech_SD"]),
            use_container_width=True,
            hide_index=True,
        )
    else:
        st.dataframe(display_avg, use_container_width=True, hide_index=True)

    d3, d4 = st.columns(2)
    with d3:
        download_dataframe_csv(avg_df, "qpcr_averaged_ct.csv", "⬇ Download CSV")
    with d4:
        download_dataframe_excel(avg_df, "qpcr_averaged_ct.xlsx", "Averaged Ct", "⬇ Download Excel")

# ── Tab 3: Per-biological-replicate detail ───────────────────────────────────
with tab_detail:
    st.markdown(
        "All calculated values per biological replicate × target gene. "
        "Each row = **one biological replicate** (technical replicates already averaged)."
    )

    show_cols = [
        c for c in [
            sample_col, gene_col, group_col,
            ct_col, "Ref_Ct", "Delta_Ct",
            "Mean_Delta_Ct_Control", "N_Control_BioReps",
            "Delta_Delta_Ct", "Fold_Change", "log2FC",
        ]
        if c in detail_df.columns
    ]
    st.dataframe(detail_df[show_cols].round(4), use_container_width=True, hide_index=True)

    d5, d6 = st.columns(2)
    with d5:
        download_dataframe_csv(detail_df, "qpcr_detail.csv", "⬇ Download CSV")
    with d6:
        download_dataframe_excel(detail_df, "qpcr_detail.xlsx", "Detail", "⬇ Download Excel")

# ── Tab 4: Fold Change bar + dot plot ────────────────────────────────────────
with tab_fc:
    st.markdown(
        "**Bars** = mean Fold Change across biological replicates ± SD.  "
        "**Dots** = individual biological replicate values, aligned to their group bar.  "
        "Significance stars shown for adjusted p < 0.05 (2-group comparisons)."
    )
    try:
        fig_fc = plot_fold_change_bar(
            plot_summary_df, plot_detail_df,
            group_col, gene_col,
            stats_df=plot_stats_df,
        )
        st.plotly_chart(fig_fc, use_container_width=True)

        d7, d8, d9 = st.columns(3)
        with d7:
            download_figure_png(
                fig_fc, "qpcr_fold_change.png",
                "⬇ PNG (300 dpi)"
            )
        with d8:
            download_figure_png_600dpi(
                fig_fc, "qpcr_fold_change_600dpi.png",
                "⬇ PNG (600 dpi, publication)"
            )
        with d9:
            download_figure_svg(fig_fc, "qpcr_fold_change.svg", "⬇ SVG (vector)")
    except Exception as e:
        st.error(f"Could not render plot: {e}")

# ── Tab 5: log₂FC bar + dot plot ─────────────────────────────────────────────
with tab_log2:
    st.markdown(
        "**Bars** = mean log₂ Fold Change ± SD.  "
        "**Dots** = individual biological replicate values, aligned to their group bar.  "
        "Bars above 0 = upregulated; below 0 = downregulated."
    )
    try:
        fig_log2 = plot_log2fc_bar(
            plot_summary_df, plot_detail_df,
            group_col, gene_col,
            stats_df=plot_stats_df,
        )
        st.plotly_chart(fig_log2, use_container_width=True)

        d10, d11, d12 = st.columns(3)
        with d10:
            download_figure_png(
                fig_log2, "qpcr_log2fc.png",
                "⬇ PNG (300 dpi)"
            )
        with d11:
            download_figure_png_600dpi(
                fig_log2, "qpcr_log2fc_600dpi.png",
                "⬇ PNG (600 dpi, publication)"
            )
        with d12:
            download_figure_svg(fig_log2, "qpcr_log2fc.svg", "⬇ SVG (vector)")
    except Exception as e:
        st.error(f"Could not render plot: {e}")

# ── Tab 6: ΔCt heatmap ───────────────────────────────────────────────────────
with tab_heatmap:
    st.markdown(
        "Each cell = ΔCt for one **biological sample** (mean Ct of target − mean Ct of reference). "
        "Blue = lower ΔCt (higher relative expression); Red = higher ΔCt."
    )
    try:
        # Filter heatmap to selected groups via the detail_df
        heatmap_df = delta_ct_df[delta_ct_df[group_col].isin(selected_groups)]
        fig_heat = plot_deltact_heatmap(heatmap_df, sample_col, gene_col)
        st.plotly_chart(fig_heat, use_container_width=True)

        d13, d14, d15 = st.columns(3)
        with d13:
            download_figure_png(
                fig_heat, "qpcr_deltact_heatmap.png",
                "⬇ PNG (300 dpi)"
            )
        with d14:
            download_figure_png_600dpi(
                fig_heat, "qpcr_deltact_heatmap_600dpi.png",
                "⬇ PNG (600 dpi, publication)"
            )
        with d15:
            download_figure_svg(
                fig_heat, "qpcr_deltact_heatmap.svg",
                "⬇ SVG (vector)"
            )
    except Exception as e:
        st.error(f"Could not render heatmap: {e}")

# ── Tab 7: Statistics ─────────────────────────────────────────────────────────
with tab_stats:
    n_groups = len(available_groups)
    test_name_map = {
        "auto":     "Welch's t-test (auto)" if n_groups == 2 else "One-way ANOVA (auto)",
        "welch":    "Welch's t-test",
        "t-test":   "Student's t-test",
        "wilcoxon": "Wilcoxon rank-sum (Mann-Whitney U)",
        "anova":    "One-way ANOVA",
        "kruskal":  "Kruskal-Wallis",
    }
    correction_name_map = {
        "fdr_bh":    "Benjamini-Hochberg FDR",
        "bonferroni":"Bonferroni",
        "none":      "None (raw p-values)",
    }
    st.markdown(
        f"**Test:** {test_name_map.get(stat_test, stat_test)}  |  "
        f"**Correction:** {correction_name_map.get(correction, correction)}  \n"
        "Tests compare **ΔCt biological replicate values** between groups "
        "(not Fold Change). Significance: \\*\\*\\*\\* p<0.0001, \\*\\*\\* p<0.001, "
        "\\*\\* p<0.01, \\* p<0.05, ns = not significant."
    )

    # Colour-coded display
    def _colour_sig(val):
        colours = {
            "****": "background-color:#C3E6CB; font-weight:bold",
            "***":  "background-color:#C3E6CB; font-weight:bold",
            "**":   "background-color:#FFF3CD; font-weight:bold",
            "*":    "background-color:#FFF3CD; font-weight:bold",
            "ns":   "color:#888888",
        }
        return colours.get(val, "")

    styled = stats_df.style.applymap(
        _colour_sig, subset=["Significance"]
    ) if "Significance" in stats_df.columns else stats_df

    st.dataframe(styled, use_container_width=True, hide_index=True)

    d16, d17, d18 = st.columns(3)
    with d16:
        download_dataframe_csv(stats_df, "qpcr_statistics.csv", "⬇ Download CSV")
    with d17:
        download_dataframe_excel(stats_df, "qpcr_statistics.xlsx", "Statistics", "⬇ Download Excel")
    with d18:
        try:
            fig_stats_tbl = plot_stats_results(stats_df)
            download_figure_png(
                fig_stats_tbl, "qpcr_statistics_table.png",
                "⬇ Stats table PNG"
            )
        except Exception:
            pass

# =============================================================================
# Methods text
# =============================================================================
st.divider()
with st.expander("📝 Methods text — copy into your manuscript"):
    n_bio = int(summary_df["N_BioReps"].max()) if "N_BioReps" in summary_df.columns else "n"
    tech_n_val = f"{avg_tech_n:.0f}" if isinstance(avg_tech_n, float) else str(avg_tech_n)
    test_full = test_name_map.get(stat_test, stat_test)
    corr_full = correction_name_map.get(correction, correction)

    st.markdown(
        f"""
        Gene expression was quantified by RT-qPCR. For each biological sample,
        {tech_n_val} technical replicates were averaged to obtain a single mean Ct per gene.
        Ct values were normalised to *{reference_gene}* as the reference gene
        (ΔCt = mean Ct_target − mean Ct_{reference_gene}).
        Relative expression was calculated using the comparative ΔΔCt method
        (Livak & Schmittgen, 2001): ΔΔCt = ΔCt_sample − mean(ΔCt_{control_group}),
        and Fold Change = 2^(−ΔΔCt).
        Data are presented as mean ± SD of n = {n_bio} biological replicates.
        Statistical comparisons between groups were performed using {test_full}
        on ΔCt values, with {corr_full} correction for multiple comparisons.
        Analysis was performed with the Research Data Analysis Toolkit (v0.1.0).
        """
    )

st.divider()
st.caption(
    "Reference: Livak KJ & Schmittgen TD (2001). Analysis of relative gene "
    "expression data using real-time quantitative PCR and the 2^(−ΔΔCT) method. "
    "*Methods* 25(4): 402–408."
)
