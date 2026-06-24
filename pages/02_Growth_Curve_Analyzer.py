# =============================================================================
# pages/02_Growth_Curve_Analyzer.py  –  Growth Curve Analyzer
# =============================================================================
# PURPOSE
#   Streamlit UI for bacterial growth curve analysis.
#   Wide-format OD data → background correction → tech-rep averaging →
#   growth metrics → group statistics → plots → export.
#
# PIPELINE ORDER (enforced by session_state gating)
#   1. parse_wide_format()           – wide CSV → long DataFrame
#   2. subtract_blank()              – per-well background correction
#   3. average_technical_replicates() – collapse tech reps → bio samples
#   4. calculate_growth_metrics()    – 9 metrics per biological sample
#   5. summarise_metrics()           – mean ± SD per group
#   6. run_statistical_tests()       – stats on biological replicates
#
# SESSION STATE KEYS  (all prefixed gc_)
#   gc_long_df          – parsed long-format raw data
#   gc_sample_names     – list of all detected sample names (incl. blanks)
#   gc_results          – dict of all analysis DataFrames + metadata
#   gc_blank_assignment – dict: sample → blank name
#   gc_neg_handling     – str: negative value handling
#   gc_group_assignment – dict: sample → group name
#   gc_sel_groups       – list: groups selected for plotting
#   gc_sel_samples      – list: samples selected for plotting
#   gc_group_order      – list: display order of groups
#   gc_palette          – str: colour palette name
#   gc_show_reps        – bool: show individual replicate lines
#   gc_show_sig         – bool: show significance stars
#   gc_sel_metrics      – list: metrics selected for bar plots
#
# All mathematics live in analysis/growth_curve.py.
# All figure construction lives in visualization/growth_plots.py.
# =============================================================================

import sys
from pathlib import Path

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
from analysis.growth_curve import (
    parse_wide_format,
    get_detected_samples,
    subtract_blank,
    average_technical_replicates,
    calculate_growth_metrics,
    summarise_metrics,
    run_statistical_tests,
)
from visualization.growth_plots import (
    plot_growth_curves,
    plot_metric_bars,
    plot_metric_panel,
    plot_stats_table,
    PALETTES,
)

# =============================================================================
# Page config
# =============================================================================
st.set_page_config(
    page_title="Growth Curve Analyzer · Research Toolkit",
    page_icon="📈",
    layout="wide",
)
render_sidebar()

st.title("📈 Growth Curve Analyzer")
st.markdown(
    """
    Upload wide-format OD₆₀₀ time-series data.
    The pipeline handles **blank media subtraction**, **technical replicate averaging**,
    **growth metric calculation**, and **publication-quality plots** — all with persistent
    session state so changing plot options never resets your analysis.
    """
)
st.divider()


# =============================================================================
# Helpers
# =============================================================================

def _stop(msg: str) -> None:
    st.error(msg)
    st.stop()


# One-time defaults for scalar session_state keys
for _key, _default in [
    ("gc_palette",   "Default"),
    ("gc_show_reps", True),
    ("gc_show_sig",  True),
    ("gc_neg_handling", "zero"),
]:
    if _key not in st.session_state:
        st.session_state[_key] = _default

METRIC_LABELS = {
    "Max_OD":          "Maximum OD₆₀₀",
    "Final_OD":        "Final OD₆₀₀",
    "AUC":             "Area Under Curve (AUC)",
    "Max_Growth_Rate": "Max Growth Rate (OD/h)",
    "Lag_Phase_h":     "Lag Phase (h)",
    "Doubling_Time_h": "Doubling Time (h)",
    "Time_to_Max_Rate":"Time to Max Rate (h)",
    "Time_to_OD_0.2":  "Time to OD 0.2 (h)",
    "Time_to_OD_0.5":  "Time to OD 0.5 (h)",
}

ALL_METRICS = list(METRIC_LABELS.keys())


# =============================================================================
# STEP 1 — Load data
# =============================================================================
st.subheader("Step 1 — Load your data")

input_mode = st.radio(
    "Data source",
    ["Upload my own file", "Use the built-in sample dataset"],
    horizontal=True,
    key="gc_input_mode",
    help=(
        "Wide format: first column = Time_h, remaining columns = sample wells. "
        "Duplicate column names indicate technical replicates of the same biological sample."
    ),
)

raw_df: pd.DataFrame | None = None
sample_path = _PROJECT_ROOT / "data" / "sample_growth_curve.csv"

if input_mode == "Upload my own file":
    up_col, hint_col = st.columns([3, 1])
    with up_col:
        uploaded = st.file_uploader(
            "Upload CSV or Excel (wide format)",
            type=["csv", "xlsx", "xls"],
            key="gc_file_uploader",
            help="First column must be Time_h. Remaining columns are sample wells.",
        )
    with hint_col:
        if sample_path.exists():
            st.markdown("**Format template:**")
            st.download_button(
                "⬇ Sample CSV",
                data=sample_path.read_bytes(),
                file_name="sample_growth_curve.csv",
                mime="text/csv",
                key="gc_dl_sample_template",
            )
    if uploaded:
        raw_df = load_uploaded_file(uploaded)
else:
    raw_df = load_sample_file(str(sample_path))
    if raw_df is not None:
        st.success(
            "✅ Sample dataset loaded — 2 groups (WT × 3 bio reps, Mutant × 3 bio reps), "
            "3 technical replicates each, 2 blank types (BHI_Blank × 2, MRS_Blank × 1), "
            "13 timepoints (0–12 h)."
        )

if raw_df is None:
    st.info("👆 Upload a file or select the sample dataset to begin.", icon="ℹ️")
    st.stop()

# ── Parse wide format ─────────────────────────────────────────────────────────
# Detect which column is time (first column, or named Time_h / Time / time)
time_col_candidates = [c for c in raw_df.columns if "time" in c.lower()]
time_col = time_col_candidates[0] if time_col_candidates else raw_df.columns[0]

try:
    long_df = parse_wide_format(raw_df, time_col=time_col)
except Exception as e:
    _stop(f"❌ Could not parse data: {e}")

sample_info = get_detected_samples(long_df)
all_sample_names = list(sample_info.keys())

# Store parsed data in session_state
st.session_state["gc_long_df"]      = long_df
st.session_state["gc_sample_names"] = all_sample_names

with st.expander(
    f"📋 Data preview — {len(raw_df):,} timepoints × {len(all_sample_names)} samples "
    f"({long_df['TechRep_Index'].max() + 1} tech reps detected per sample on average)",
    expanded=True,
):
    sc1, sc2 = st.columns(2)
    with sc1:
        st.markdown("**Detected samples and tech-rep counts:**")
        sdf = pd.DataFrame(
            [(s, n) for s, n in sample_info.items()],
            columns=["Sample", "N_TechReps"]
        )
        st.dataframe(sdf, use_container_width=True, hide_index=True)
    with sc2:
        st.markdown("**Long-format preview (first 20 rows):**")
        st.dataframe(long_df.head(20), use_container_width=True, hide_index=True)

st.divider()


# =============================================================================
# STEP 2 — Identify blanks and assign to samples
# =============================================================================
st.subheader("Step 2 — Identify blank media controls")
st.markdown(
    "Select which samples are **blank media controls**. "
    "You will then assign each blank to the samples it should correct."
)

blank_samples = st.multiselect(
    "Blank / media control samples",
    options=all_sample_names,
    default=[s for s in all_sample_names if any(
        kw in s.lower() for kw in ["blank", "media", "bhi", "mrs", "lmrs", "control"]
    )],
    key="gc_blank_samples",
    help=(
        "Select all wells that are blank media controls. "
        "These will be used for background correction and excluded from analysis."
    ),
)

biological_samples = [s for s in all_sample_names if s not in blank_samples]

if not biological_samples:
    st.warning(
        "⚠️ All samples are marked as blanks. "
        "Please deselect at least one sample as a biological sample."
    )
    st.stop()

st.divider()


# =============================================================================
# STEP 3 — Assign blanks to samples + negative value handling
# =============================================================================
st.subheader("Step 3 — Assign blank to each sample & configure correction")

if blank_samples:
    st.markdown(
        "For each biological sample, select which blank media to subtract. "
        "Leave as **'None (no correction)'** to skip correction for that sample."
    )

    blank_options = ["None (no correction)"] + blank_samples
    blank_assignment: dict = {}

    # Lay out assignment selectors in a 3-column grid
    n_bio = len(biological_samples)
    cols_per_row = 3
    for row_start in range(0, n_bio, cols_per_row):
        row_samples = biological_samples[row_start: row_start + cols_per_row]
        cols = st.columns(cols_per_row)
        for ci, sample in enumerate(row_samples):
            with cols[ci]:
                # Try to guess the right blank by partial name match
                default_blank = "None (no correction)"
                for b in blank_samples:
                    # simple heuristic: if blank name is a substring of sample or vice versa
                    b_base = b.lower().replace("_blank", "").replace("blank", "").strip("_")
                    if b_base and b_base in sample.lower():
                        default_blank = b
                        break
                    # or if only one blank exists, use it
                    if len(blank_samples) == 1:
                        default_blank = blank_samples[0]

                chosen = st.selectbox(
                    f"Blank for **{sample}**",
                    options=blank_options,
                    index=blank_options.index(default_blank),
                    key=f"gc_blank_assign_{sample}",
                )
                if chosen != "None (no correction)":
                    blank_assignment[sample] = chosen

    st.markdown("---")
    neg_handling = st.selectbox(
        "How to handle negative corrected OD values",
        options=["zero", "keep", "na"],
        index=["zero", "keep", "na"].index(st.session_state.get("gc_neg_handling", "zero")),
        key="gc_neg_handling",
        format_func=lambda x: {
            "zero": "Set to zero (recommended — biologically sensible)",
            "keep": "Keep negative values",
            "na":   "Replace with NaN (excluded from averages)",
        }[x],
        help=(
            "Negative values arise when sample OD is slightly below blank OD due to noise. "
            "Setting to zero is standard practice."
        ),
    )
else:
    st.info(
        "ℹ️ No blank samples identified. Analysis will proceed without background correction.",
        icon="ℹ️",
    )
    blank_assignment = {}
    neg_handling     = "zero"

st.divider()


# =============================================================================
# STEP 4 — Group assignment
# =============================================================================
st.subheader("Step 4 — Assign biological samples to groups")
st.markdown(
    "Assign each biological sample to an experimental group. "
    "Samples in the same group are treated as **biological replicates** for statistics."
)

# Auto-detect groups from sample names (strip trailing _R1, _R2, _R3, _1, _2, _3)
import re as _re

def _auto_group(sample_name: str) -> str:
    """Strip trailing replicate suffix to guess a group name."""
    stripped = _re.sub(r"[_\-\s]*(R\d+|\d+)$", "", sample_name, flags=_re.IGNORECASE).strip()
    return stripped if stripped else sample_name

group_assignment: dict = {}
detected_auto_groups = {s: _auto_group(s) for s in biological_samples}
all_auto_group_names = sorted(set(detected_auto_groups.values()))

# Show a compact editable table: sample → group
st.markdown(
    "Groups are auto-detected from sample names. "
    "Edit the text boxes below to override."
)

n_bio = len(biological_samples)
cols_per_row = 3
for row_start in range(0, n_bio, cols_per_row):
    row_samples = biological_samples[row_start: row_start + cols_per_row]
    cols = st.columns(cols_per_row)
    for ci, sample in enumerate(row_samples):
        with cols[ci]:
            auto_grp = detected_auto_groups[sample]
            grp = st.text_input(
                f"Group for **{sample}**",
                value=auto_grp,
                key=f"gc_group_{sample}",
            )
            group_assignment[sample] = grp if grp.strip() else auto_grp

# Summary
final_groups = sorted(set(group_assignment.values()))
with st.expander("🔍 Group assignment summary"):
    gdf = pd.DataFrame(
        [(s, g) for s, g in group_assignment.items()],
        columns=["Sample", "Group"]
    )
    st.dataframe(gdf, use_container_width=True, hide_index=True)
    st.markdown(f"**{len(final_groups)} groups detected:** {', '.join(final_groups)}")

st.divider()


# =============================================================================
# STEP 5 — Analysis parameters + Run
# =============================================================================
st.subheader("Step 5 — Analysis parameters")

a1, a2 = st.columns(2)
with a1:
    stat_test = st.selectbox(
        "Statistical test",
        ["auto", "welch", "t-test", "wilcoxon", "anova", "kruskal"],
        index=0,
        key="gc_stat_test",
        format_func=lambda x: {
            "auto":     "Auto (Welch / ANOVA)",
            "welch":    "Welch's t-test (unequal variance)",
            "t-test":   "Student's t-test (equal variance)",
            "wilcoxon": "Wilcoxon rank-sum (non-parametric)",
            "anova":    "One-way ANOVA",
            "kruskal":  "Kruskal-Wallis (non-parametric)",
        }[x],
        help="Tests compare biological replicate metric values between groups.",
    )
with a2:
    correction = st.selectbox(
        "Multiple testing correction",
        ["fdr_bh", "bonferroni", "none"],
        index=0,
        key="gc_correction",
        format_func=lambda x: {
            "fdr_bh":    "Benjamini-Hochberg (FDR)",
            "bonferroni":"Bonferroni",
            "none":      "None (raw p-values)",
        }[x],
    )

run_clicked = st.button("▶ Run Growth Curve Analysis", type="primary", key="gc_run_btn")

if run_clicked:
    with st.spinner("Running analysis pipeline…"):
        try:
            # ── Step 2: Background correction ────────────────────────────────
            corrected_df = subtract_blank(
                long_df,
                blank_assignment=blank_assignment,
                negative_handling=st.session_state["gc_neg_handling"],
            )

            # ── Step 3: Average technical replicates ──────────────────────────
            avg_df = average_technical_replicates(
                corrected_df,
                group_assignment=group_assignment,
            )

            # ── Step 4: Growth metrics ────────────────────────────────────────
            metrics_df = calculate_growth_metrics(avg_df)

            # ── Step 5: Group summary ─────────────────────────────────────────
            summary_df = summarise_metrics(metrics_df)

            # ── Step 6: Statistics ────────────────────────────────────────────
            try:
                stats_df = run_statistical_tests(
                    metrics_df,
                    test=stat_test,
                    correction_method=correction,
                )
            except ValueError as ve:
                stats_df = pd.DataFrame()
                st.warning(f"⚠️ Statistical testing skipped: {ve}")

            # ── Store results in session_state ────────────────────────────────
            st.session_state["gc_results"] = {
                "corrected_df":    corrected_df,
                "avg_df":          avg_df,
                "metrics_df":      metrics_df,
                "summary_df":      summary_df,
                "stats_df":        stats_df,
                "blank_assignment": blank_assignment,
                "group_assignment": group_assignment,
                "stat_test":       stat_test,
                "correction":      correction,
                "neg_handling":    st.session_state["gc_neg_handling"],
            }

            # Reset plot-option selections to "all" for a fresh run
            all_groups  = sorted(set(group_assignment.values()))
            all_samples = list(group_assignment.keys())
            st.session_state["gc_sel_groups"]  = all_groups
            st.session_state["gc_sel_samples"] = all_samples
            st.session_state["gc_group_order"] = all_groups
            st.session_state["gc_sel_metrics"] = ALL_METRICS

        except Exception as e:
            _stop(f"❌ Analysis error: {e}\n\nCheck your data format and group assignments.")

# ── Gate: require results before showing anything below ───────────────────────
if "gc_results" not in st.session_state:
    st.info("Configure settings above, then click **▶ Run Growth Curve Analysis**.", icon="ℹ️")
    st.stop()

# ── Unpack results ────────────────────────────────────────────────────────────
_r            = st.session_state["gc_results"]
corrected_df  = _r["corrected_df"]
avg_df        = _r["avg_df"]
metrics_df    = _r["metrics_df"]
summary_df    = _r["summary_df"]
stats_df      = _r["stats_df"]
stat_test     = _r["stat_test"]
correction    = _r["correction"]

all_result_groups  = sorted(metrics_df["Group"].unique().tolist())
all_result_samples = sorted(metrics_df["Sample"].unique().tolist())

st.success("✅ Analysis complete!", icon="✅")

n_bio_samples = metrics_df["Sample"].nunique()
n_timepoints  = avg_df["Time_h"].nunique()
st.info(
    f"ℹ️ **{n_bio_samples} biological samples** across "
    f"**{len(all_result_groups)} groups** · "
    f"**{n_timepoints} timepoints** · "
    f"Metrics calculated per biological sample.",
    icon="ℹ️",
)

st.divider()


# =============================================================================
# STEP 6 — Plot options (all session_state-backed with stable keys)
# =============================================================================
st.subheader("Step 6 — Plot options")
st.divider()

# Sanitise session_state lists against current results
for _k, _all in [("gc_sel_groups", all_result_groups), ("gc_sel_samples", all_result_samples),
                 ("gc_group_order", all_result_groups), ("gc_sel_metrics", ALL_METRICS)]:
    if _k not in st.session_state:
        st.session_state[_k] = _all
    else:
        st.session_state[_k] = [v for v in st.session_state[_k] if v in _all] or _all

# ── A: Group selection & order ────────────────────────────────────────────────
st.markdown("**Group display & order**")
st.caption("Select groups and the order they appear in plots (first selected = leftmost).")

st.multiselect(
    "Groups to display (select in desired order)",
    options=all_result_groups,
    default=st.session_state["gc_sel_groups"],
    key="gc_sel_groups",
    help="Download tables always contain all groups.",
)
selected_groups = st.session_state["gc_sel_groups"]
group_order     = selected_groups

if not selected_groups:
    st.warning("⚠️ No groups selected.", icon="⚠️")
    st.stop()

# ── B: Sample selection ───────────────────────────────────────────────────────
st.markdown("**Sample selection**")
samples_in_selected_groups = [
    s for s in all_result_samples
    if metrics_df.loc[metrics_df["Sample"] == s, "Group"].iloc[0] in selected_groups
]
# Sanitise
valid_sel_samples = [s for s in st.session_state["gc_sel_samples"] if s in samples_in_selected_groups]
if not valid_sel_samples:
    valid_sel_samples = samples_in_selected_groups
st.session_state["gc_sel_samples"] = valid_sel_samples

st.multiselect(
    "Biological samples to display in growth curve plot",
    options=samples_in_selected_groups,
    default=st.session_state["gc_sel_samples"],
    key="gc_sel_samples",
    help="Affects growth curve plot only. All samples included in metric and stat tables.",
)
selected_samples = st.session_state["gc_sel_samples"]

st.divider()

# ── C / D / E: Metric selection, palette, toggles ─────────────────────────────
opt1, opt2, opt3 = st.columns([2, 1, 1])

with opt1:
    st.markdown("**Metrics to display in bar plots**")
    st.multiselect(
        "Growth metrics shown in plots",
        options=ALL_METRICS,
        default=st.session_state["gc_sel_metrics"],
        key="gc_sel_metrics",
        format_func=lambda m: METRIC_LABELS.get(m, m),
        help="Affects bar charts only. All metrics are always available for download.",
    )
    selected_metrics = st.session_state["gc_sel_metrics"]
    if not selected_metrics:
        st.warning("⚠️ Select at least one metric.")
        st.stop()

with opt2:
    st.markdown("**Color palette**")
    st.selectbox(
        "Color scheme",
        options=list(PALETTES.keys()),
        index=list(PALETTES.keys()).index(st.session_state.get("gc_palette", "Default")),
        key="gc_palette",
        help=(
            "**Default** — teal/coral.\n"
            "**Nature-style** — Nature journal colors.\n"
            "**Colorblind-friendly** — Wong (2011).\n"
            "**Pastel** / **High contrast** / **Viridis**."
        ),
    )
    palette_name = st.session_state["gc_palette"]

with opt3:
    st.markdown("**Display options**")
    st.checkbox(
        "Show individual biological replicates",
        value=st.session_state.get("gc_show_reps", True),
        key="gc_show_reps",
        help="Overlay thin dotted lines for each biological sample in the growth curve plot.",
    )
    st.checkbox(
        "Show significance stars",
        value=st.session_state.get("gc_show_sig", True),
        key="gc_show_sig",
        help="Display significance annotations on metric bar plots (2-group comparisons).",
    )
    show_reps = st.session_state["gc_show_reps"]
    show_sig  = st.session_state["gc_show_sig"]

st.divider()


# =============================================================================
# STEP 7 — Results tabs
# =============================================================================
st.subheader("Step 7 — Results")

(
    tab_curves, tab_metrics_panel, tab_metric_detail,
    tab_corrected, tab_averaged, tab_summary,
    tab_stats, tab_methods,
) = st.tabs([
    "📈 Growth Curves",
    "📊 Metrics Overview",
    "📊 Individual Metric",
    "🔬 Corrected OD Table",
    "📋 Averaged OD Table",
    "📄 Group Summary",
    "📈 Statistics",
    "📝 Methods",
])

# ── Tab 1: Growth curves ──────────────────────────────────────────────────────
with tab_curves:
    st.markdown(
        "**Thick lines** = group mean OD ± SD (shaded ribbon).  "
        "**Dotted lines** = individual biological replicates (toggle above)."
    )
    try:
        # Filter avg_df to selected samples
        plot_avg_df = avg_df[
            avg_df["Group"].isin(selected_groups) &
            avg_df["Sample"].isin(selected_samples)
        ]
        fig_curves = plot_growth_curves(
            plot_avg_df,
            selected_groups=selected_groups,
            selected_samples=selected_samples,
            palette=palette_name,
            show_replicates=show_reps,
        )
        st.plotly_chart(fig_curves, use_container_width=True)

        c1, c2, c3 = st.columns(3)
        with c1:
            download_figure_png(fig_curves, "growth_curves.png", "⬇ PNG (300 dpi)")
        with c2:
            download_figure_png_600dpi(fig_curves, "growth_curves_600dpi.png", "⬇ PNG (600 dpi)")
        with c3:
            download_figure_svg(fig_curves, "growth_curves.svg", "⬇ SVG (vector)")
    except Exception as e:
        st.error(f"Could not render growth curves: {e}")

# ── Tab 2: Metrics overview panel ────────────────────────────────────────────
with tab_metrics_panel:
    st.markdown(
        "Bar + dot plots for all selected metrics. "
        "**Bars** = group mean ± SD. **Dots** = individual biological replicates."
    )
    try:
        plot_metrics_df = metrics_df[metrics_df["Group"].isin(selected_groups)]
        fig_panel = plot_metric_panel(
            plot_metrics_df,
            selected_metrics=selected_metrics,
            group_order=group_order,
            palette=palette_name,
            show_annotations=show_sig,
            stats_df=stats_df if not stats_df.empty else None,
        )
        st.plotly_chart(fig_panel, use_container_width=True)

        p1, p2, p3 = st.columns(3)
        with p1:
            download_figure_png(fig_panel, "metrics_panel.png", "⬇ PNG (300 dpi)")
        with p2:
            download_figure_png_600dpi(fig_panel, "metrics_panel_600dpi.png", "⬇ PNG (600 dpi)")
        with p3:
            download_figure_svg(fig_panel, "metrics_panel.svg", "⬇ SVG (vector)")
    except Exception as e:
        st.error(f"Could not render metrics panel: {e}")

# ── Tab 3: Single metric deep-dive ───────────────────────────────────────────
with tab_metric_detail:
    st.markdown("Select a single metric for a detailed bar + dot plot.")

    chosen_metric = st.selectbox(
        "Select metric",
        options=selected_metrics,
        format_func=lambda m: METRIC_LABELS.get(m, m),
        key="gc_chosen_metric",
    )
    try:
        plot_metrics_df2 = metrics_df[metrics_df["Group"].isin(selected_groups)]
        fig_single = plot_metric_bars(
            plot_metrics_df2,
            metric=chosen_metric,
            group_order=group_order,
            palette=palette_name,
            show_annotations=show_sig,
            stats_df=stats_df if not stats_df.empty else None,
            y_label=METRIC_LABELS.get(chosen_metric, chosen_metric),
            title=METRIC_LABELS.get(chosen_metric, chosen_metric),
        )
        st.plotly_chart(fig_single, use_container_width=True)

        m1, m2, m3 = st.columns(3)
        with m1:
            download_figure_png(fig_single, f"{chosen_metric}.png", "⬇ PNG (300 dpi)")
        with m2:
            download_figure_png_600dpi(fig_single, f"{chosen_metric}_600dpi.png", "⬇ PNG (600 dpi)")
        with m3:
            download_figure_svg(fig_single, f"{chosen_metric}.svg", "⬇ SVG (vector)")
    except Exception as e:
        st.error(f"Could not render {chosen_metric}: {e}")

# ── Tab 4: Background-corrected OD table ─────────────────────────────────────
with tab_corrected:
    st.markdown(
        "Per-well OD after background subtraction (technical replicates still separate). "
        "Use this to verify blank assignment."
    )
    display_corr = corrected_df[corrected_df["Blank_Used"] != "self"].copy()
    st.dataframe(display_corr.round(4), use_container_width=True, hide_index=True)

    tc1, tc2 = st.columns(2)
    with tc1:
        download_dataframe_csv(corrected_df, "corrected_od.csv", "⬇ Download CSV")
    with tc2:
        download_dataframe_excel(corrected_df, "corrected_od.xlsx", "Corrected OD", "⬇ Download Excel")

# ── Tab 5: Averaged OD table ─────────────────────────────────────────────────
with tab_averaged:
    st.markdown(
        "Mean OD per **biological sample × timepoint** after averaging technical replicates."
    )
    st.dataframe(avg_df.round(4), use_container_width=True, hide_index=True)

    ta1, ta2 = st.columns(2)
    with ta1:
        download_dataframe_csv(avg_df, "averaged_od.csv", "⬇ Download CSV")
    with ta2:
        download_dataframe_excel(avg_df, "averaged_od.xlsx", "Averaged OD", "⬇ Download Excel")

# ── Tab 6: Group summary ─────────────────────────────────────────────────────
with tab_summary:
    st.markdown(
        "Mean ± SD per **group × metric** across biological replicates. "
        "Use this table for publication-ready summary statistics."
    )

    # Pivot to wide format for readability: rows = groups, cols = metrics
    try:
        pivot = summary_df.pivot_table(
            index="Group", columns="Metric",
            values=["Mean", "SD", "N"],
            aggfunc="first",
        )
        pivot.columns = [f"{m}_{stat}" for stat, m in pivot.columns]
        st.dataframe(pivot.round(4), use_container_width=True)
    except Exception:
        st.dataframe(summary_df.round(4), use_container_width=True, hide_index=True)

    ts1, ts2 = st.columns(2)
    with ts1:
        download_dataframe_csv(summary_df, "group_summary.csv", "⬇ Download CSV")
    with ts2:
        download_dataframe_excel(summary_df, "group_summary.xlsx", "Group Summary", "⬇ Download Excel")

    # Per-biological-sample metrics table (for download)
    st.markdown("---")
    st.markdown("**Per-biological-sample metrics (all groups):**")
    st.dataframe(metrics_df.round(4), use_container_width=True, hide_index=True)
    tm1, tm2 = st.columns(2)
    with tm1:
        download_dataframe_csv(metrics_df, "per_sample_metrics.csv", "⬇ Download CSV")
    with tm2:
        download_dataframe_excel(metrics_df, "per_sample_metrics.xlsx",
                                 "Per-Sample Metrics", "⬇ Download Excel")

# ── Tab 7: Statistics ─────────────────────────────────────────────────────────
with tab_stats:
    test_name_map = {
        "auto":     "Welch's t-test (auto)" if len(all_result_groups) == 2 else "One-way ANOVA (auto)",
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

    if stats_df is None or stats_df.empty:
        st.warning(
            "⚠️ Statistical results not available. "
            "Ensure you have ≥2 groups with ≥2 biological replicates each."
        )
    else:
        st.markdown(
            f"**Test:** {test_name_map.get(stat_test, stat_test)}  |  "
            f"**Correction:** {correction_name_map.get(correction, correction)}  \n"
            "Tests compare **metric values across biological replicates** between groups. "
            "Significance: \\*\\*\\*\\* p<0.0001, \\*\\*\\* p<0.001, "
            "\\*\\* p<0.01, \\* p<0.05, ns = not significant."
        )

        def _colour_sig(val):
            colours = {
                "****": "background-color:#C3E6CB; font-weight:bold",
                "***":  "background-color:#C3E6CB; font-weight:bold",
                "**":   "background-color:#FFF3CD; font-weight:bold",
                "*":    "background-color:#FFF3CD; font-weight:bold",
                "ns":   "color:#888888",
            }
            return colours.get(val, "")

        styled = (
            stats_df.style.applymap(_colour_sig, subset=["Significance"])
            if "Significance" in stats_df.columns else stats_df
        )
        st.dataframe(styled, use_container_width=True, hide_index=True)

        sts1, sts2, sts3 = st.columns(3)
        with sts1:
            download_dataframe_csv(stats_df, "growth_statistics.csv", "⬇ Download CSV")
        with sts2:
            download_dataframe_excel(stats_df, "growth_statistics.xlsx",
                                     "Statistics", "⬇ Download Excel")
        with sts3:
            try:
                fig_stats_tbl = plot_stats_table(stats_df)
                download_figure_png(fig_stats_tbl, "statistics_table.png", "⬇ Stats table PNG")
            except Exception:
                pass

# ── Tab 8: Methods text ───────────────────────────────────────────────────────
with tab_methods:
    n_bio_total   = metrics_df["Sample"].nunique()
    n_grps        = len(all_result_groups)
    blank_desc    = (
        f"Background subtraction was performed using {len(set(blank_assignment.values()))} "
        f"blank media control(s) ({', '.join(set(blank_assignment.values()))}), "
        f"corrected values below zero were set to {st.session_state['gc_results']['neg_handling']}. "
        if blank_assignment else
        "No blank media correction was applied. "
    )
    test_full = test_name_map.get(stat_test, stat_test)
    corr_full = correction_name_map.get(correction, correction)

    st.markdown(
        f"""
        **Methods text** (copy into manuscript):

        OD₆₀₀ measurements were recorded at {n_timepoints} timepoints over the growth experiment.
        {blank_desc}
        Where multiple technical replicates were present, they were averaged per biological
        sample and timepoint prior to analysis. Growth metrics were calculated for each of
        {n_bio_total} biological samples across {n_grps} group(s):
        maximum OD₆₀₀ (Max_OD), area under the growth curve (AUC, trapezoidal rule),
        maximum growth rate (Max_Growth_Rate, OD₆₀₀/h from the first derivative of
        Savitzky-Golay smoothed data), lag phase duration (Lag_Phase_h, tangent-line method),
        and doubling time (Doubling_Time_h = ln(2)/μmax).
        Statistical comparisons between groups were performed using {test_full}
        on per-biological-sample metric values, with {corr_full} correction for multiple comparisons.
        Analysis was performed with the Research Data Analysis Toolkit (v0.1.0).
        """
    )

st.divider()

# ── Optional debug expander ───────────────────────────────────────────────────
with st.expander("🛠 Debug — session state (optional)", expanded=False):
    gc_keys = {k: v for k, v in st.session_state.items() if k.startswith("gc_")}
    scalar_keys = {
        k: (str(v)[:120] if not isinstance(v, (dict, pd.DataFrame))
            else ("{dict}" if isinstance(v, dict) else f"DataFrame {getattr(v, 'shape', '')}"))
        for k, v in gc_keys.items()
    }
    st.json(scalar_keys)

st.caption(
    "Reference: Zwietering MH et al. (1990). Modeling of the bacterial growth curve. "
    "*Applied and Environmental Microbiology* 56(6): 1875–1881."
)
