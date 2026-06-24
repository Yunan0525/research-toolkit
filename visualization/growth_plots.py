# =============================================================================
# visualization/growth_plots.py  –  Growth Curve Plotly Figure Builders
# =============================================================================
# PURPOSE
#   Build publication-quality interactive Plotly figures from DataFrames
#   produced by analysis/growth_curve.py.  No Streamlit code; no analysis math.
#
# INPUT CONTRACT
#   All plot functions receive data AFTER technical replicates have been
#   averaged.  One trace per biological sample or per group mean.
#
# FIGURES
#   plot_growth_curves       – OD vs time: group mean ± SD ribbon + optional
#                              individual biological replicate lines
#   plot_metric_bars         – bar + dot chart for a single growth metric
#                              (mirrors qPCR fold-change style)
#   plot_metric_panel        – grid of bar+dot subplots for all selected metrics
#   plot_stats_table         – Plotly Table of statistical results
#
# STYLE (Nature-journal conventions, matching qPCR module)
#   Font   : Arial, sans-serif  (14 pt title, 13 pt axis, 11 pt ticks)
#   Errors : SD of biological replicates, cap width 8 px
#   Grid   : horizontal only, #EEEEEE, 1 px
#   BG     : white plot area, transparent paper
# =============================================================================

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np


# =============================================================================
# Shared style constants  (kept in sync with qpcr_plots.py)
# =============================================================================

PALETTES = {
    "Default": [
        "#2E8B8B", "#E07B6A", "#4C72B0", "#DD8452",
        "#55A868", "#C44E52", "#8172B3", "#937860",
    ],
    "Nature-style": [
        "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
        "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
    ],
    "Colorblind-friendly": [
        "#0072B2", "#E69F00", "#009E73", "#CC79A7",
        "#56B4E9", "#D55E00", "#F0E442", "#000000",
    ],
    "Pastel": [
        "#AEC6CF", "#FFD1DC", "#B5EAD7", "#FFDAC1",
        "#C7CEEA", "#E2F0CB", "#F8C8D4", "#BEE3F8",
    ],
    "High contrast": [
        "#000000", "#E6194B", "#3CB44B", "#4363D8",
        "#F58231", "#911EB4", "#42D4F4", "#F032E6",
    ],
    "Viridis": [
        "#440154", "#31688E", "#35B779", "#FDE725",
        "#21908C", "#5DC863", "#443983", "#90D743",
    ],
}

PALETTE = PALETTES["Default"]   # convenience alias

_FONT = "Arial, sans-serif"

_AXIS = dict(
    showgrid=False,
    linecolor="#CCCCCC",
    linewidth=1,
    ticks="outside",
    tickfont=dict(family=_FONT, size=11, color="#2C2C2C"),
    title_font=dict(family=_FONT, size=13, color="#2C2C2C"),
)

_LAYOUT = dict(
    font=dict(family=_FONT, size=12, color="#2C2C2C"),
    plot_bgcolor="#FFFFFF",
    paper_bgcolor="rgba(0,0,0,0)",
    margin=dict(l=70, r=40, t=80, b=60),
    legend=dict(
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="#CCCCCC",
        borderwidth=1,
        font=dict(family=_FONT, size=11),
    ),
)


def _title_style(text: str) -> dict:
    return dict(text=text, font=dict(family=_FONT, size=14, color="#1F4E79"), x=0.05)


def _y_axis_style(title: str, **extra) -> dict:
    base = dict(
        title=title,
        showgrid=True,
        gridcolor="#EEEEEE",
        gridwidth=1,
        **{k: v for k, v in _AXIS.items() if k != "showgrid"},
    )
    base.update(extra)
    return base


# =============================================================================
# Figure 1 — Growth curves: mean ± SD ribbon + optional individual lines
# =============================================================================

def plot_growth_curves(
    avg_df: pd.DataFrame,
    group_col: str = "Group",
    show_replicates: bool = True,
    selected_groups: list | None = None,
    selected_samples: list | None = None,
    palette: str = "Default",
    title: str = "OD₆₀₀ Growth Curves  —  mean ± SD",
) -> go.Figure:
    """
    Plot OD vs time for each group.

    Mean line + SD ribbon (shaded) per group.
    Optional: overlay individual biological replicate traces.

    Parameters
    ----------
    avg_df : pd.DataFrame
        Output of average_technical_replicates(). Required columns:
        Time_h, Sample, Group, Mean_OD, SD_OD.
    group_col : str
        Column name for group labels.
    show_replicates : bool
        If True, draw thin semi-transparent lines for each biological sample.
    selected_groups : list or None
        Groups to include. None = all groups.
    selected_samples : list or None
        Specific samples to include. None = all samples in selected groups.
    palette : str
        Palette name from PALETTES dict.
    title : str
        Figure title.

    Returns
    -------
    go.Figure
    """
    df = avg_df.copy()

    # Apply group filter
    if selected_groups is not None:
        df = df[df[group_col].isin(selected_groups)]
    if selected_samples is not None:
        df = df[df["Sample"].isin(selected_samples)]

    groups = df[group_col].unique().tolist()
    active_palette = PALETTES.get(palette, PALETTES["Default"])
    fig = go.Figure()

    for i, group in enumerate(groups):
        color     = active_palette[i % len(active_palette)]
        color_rgb = _hex_to_rgba(color, 0.15)   # ribbon fill
        grp       = df[df[group_col] == group].sort_values("Time_h")

        # ── Group mean and SD per timepoint ───────────────────────────────────
        grp_mean = (
            grp.groupby("Time_h", sort=True)
            .agg(mean_od=("Mean_OD", "mean"), sd_od=("Mean_OD", "std"))
            .reset_index()
        )

        time_fwd  = grp_mean["Time_h"].tolist()
        mean_fwd  = grp_mean["mean_od"].tolist()
        sd_vals   = grp_mean["sd_od"].fillna(0).tolist()
        upper     = [m + s for m, s in zip(mean_fwd, sd_vals)]
        lower     = [max(0, m - s) for m, s in zip(mean_fwd, sd_vals)]

        # SD ribbon (filled area between upper and lower)
        fig.add_trace(go.Scatter(
            x=time_fwd + time_fwd[::-1],
            y=upper + lower[::-1],
            fill="toself",
            fillcolor=color_rgb,
            line=dict(color="rgba(0,0,0,0)"),
            showlegend=False,
            hoverinfo="skip",
            legendgroup=group,
        ))

        # Mean line
        fig.add_trace(go.Scatter(
            x=grp_mean["Time_h"],
            y=grp_mean["mean_od"],
            mode="lines",
            name=group,
            line=dict(color=color, width=2.5),
            legendgroup=group,
            hovertemplate=(
                f"<b>{group}</b><br>"
                "Time: %{x:.1f} h<br>"
                "Mean OD: %{y:.4f}<br>"
                "<extra></extra>"
            ),
        ))

        # Individual biological replicate lines (optional)
        if show_replicates:
            samples_in_group = grp["Sample"].unique().tolist()
            for sample in samples_in_group:
                sdf = grp[grp["Sample"] == sample].sort_values("Time_h")
                fig.add_trace(go.Scatter(
                    x=sdf["Time_h"],
                    y=sdf["Mean_OD"],
                    mode="lines+markers",
                    name=sample,
                    line=dict(color=color, width=1, dash="dot"),
                    marker=dict(size=5, color=color, opacity=0.7),
                    opacity=0.55,
                    legendgroup=group,
                    showlegend=True,
                    hovertemplate=(
                        f"<b>{sample}</b> ({group})<br>"
                        "Time: %{x:.1f} h<br>"
                        "OD: %{y:.4f}<br>"
                        "<extra></extra>"
                    ),
                ))

    fig.update_layout(
        xaxis=dict(title="Time (h)", **_AXIS),
        yaxis=_y_axis_style("OD₆₀₀", rangemode="tozero"),
        title=_title_style(title),
        **_LAYOUT,
    )
    return fig


# =============================================================================
# Figure 2 — Bar + dot chart for a single growth metric
# =============================================================================

def plot_metric_bars(
    metrics_df: pd.DataFrame,
    metric: str,
    group_col: str = "Group",
    group_order: list | None = None,
    palette: str = "Default",
    show_annotations: bool = True,
    stats_df: pd.DataFrame | None = None,
    y_label: str | None = None,
    title: str | None = None,
) -> go.Figure:
    """
    Bar + dot plot for a single growth metric, mirroring the qPCR style.

    Bars  — group mean.
    Error — SD across biological replicates.
    Dots  — individual biological sample values.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        Output of calculate_growth_metrics(). One row per biological sample.
    metric : str
        Column name of the metric to plot (e.g. 'Max_OD', 'Lag_Phase_h').
    group_col : str
        Column name for group labels.
    group_order : list or None
        Explicit display order for groups.
    palette : str
        Palette name from PALETTES dict.
    show_annotations : bool
        Draw significance stars if stats_df is provided.
    stats_df : pd.DataFrame or None
        Output of run_statistical_tests(). Used for significance stars.
    y_label : str or None
        Y-axis label. Defaults to metric name.
    title : str or None
        Figure title. Defaults to metric name.

    Returns
    -------
    go.Figure
    """
    df = metrics_df.dropna(subset=[metric]).copy()
    if df.empty:
        return _empty_figure(f"No data for {metric}")

    active_palette = PALETTES.get(palette, PALETTES["Default"])

    # Resolve group order
    available_groups = df[group_col].unique().tolist()
    if group_order:
        groups = [g for g in group_order if g in available_groups]
    else:
        groups = available_groups
    if not groups:
        groups = available_groups

    # Compute group stats
    summary = (
        df.groupby(group_col)[metric]
        .agg(mean_val="mean", sd_val="std", n_val="count")
        .reindex(groups)
        .reset_index()
    )

    fig = go.Figure()

    for i, group in enumerate(groups):
        color = active_palette[i % len(active_palette)]
        row   = summary[summary[group_col] == group]
        if row.empty:
            continue
        mean_val = float(row["mean_val"].iloc[0])
        sd_val   = float(row["sd_val"].iloc[0]) if not pd.isna(row["sd_val"].iloc[0]) else 0

        # Bar
        fig.add_trace(go.Bar(
            name=group,
            x=[group],
            y=[mean_val],
            error_y=dict(type="data", array=[sd_val], visible=True,
                         color="#444444", thickness=1.8, width=8),
            marker=dict(color=color, opacity=0.82,
                        line=dict(color="#333333", width=0.8)),
            legendgroup=group,
            hovertemplate=(
                f"<b>{group}</b><br>"
                f"{metric}: %{{y:.4f}} ± SD<br>"
                "<extra></extra>"
            ),
        ))

        # Dots (individual biological replicates)
        bio_vals = df[df[group_col] == group][metric].dropna().tolist()
        n_pts    = len(bio_vals)
        if n_pts > 0:
            jitter = np.linspace(-0.08, 0.08, n_pts) if n_pts > 1 else [0.0]
            fig.add_trace(go.Scatter(
                x=[group] * n_pts,
                y=bio_vals,
                mode="markers",
                name=group,
                legendgroup=group,
                showlegend=False,
                marker=dict(color=color, size=8, opacity=0.95,
                            line=dict(color="white", width=1.2)),
                hovertemplate=(
                    f"<b>{group}</b><br>"
                    f"{metric}: %{{y:.4f}}<br>"
                    "Biological replicate<extra></extra>"
                ),
            ))

    # Significance annotations (2-group only)
    if show_annotations and stats_df is not None and not stats_df.empty and len(groups) == 2:
        metric_row = stats_df[stats_df["Metric"] == metric]
        if not metric_row.empty:
            stars = metric_row["Significance"].iloc[0]
            if stars not in ("ns", "n/a"):
                top = summary["mean_val"].fillna(0) + summary["sd_val"].fillna(0)
                y_annot = top.max() * 1.12 if top.max() > 0 else 0.3
                mid_x   = groups[0] if len(groups) == 1 else 0.5
                fig.add_annotation(
                    x=mid_x if len(groups) != 2 else groups[1],
                    y=y_annot,
                    text=stars,
                    showarrow=False,
                    font=dict(family=_FONT, size=14, color="#222222"),
                    xanchor="center",
                )

    y_title = y_label or metric.replace("_", " ")
    fig.update_layout(
        xaxis=dict(title="Group", **_AXIS),
        yaxis=_y_axis_style(y_title, rangemode="tozero"),
        title=_title_style(title or metric.replace("_", " ")),
        barmode="group",
        **_LAYOUT,
    )
    return fig


# =============================================================================
# Figure 3 — Multi-metric panel (subplots grid)
# =============================================================================

def plot_metric_panel(
    metrics_df: pd.DataFrame,
    selected_metrics: list,
    group_col: str = "Group",
    group_order: list | None = None,
    palette: str = "Default",
    show_annotations: bool = True,
    stats_df: pd.DataFrame | None = None,
    ncols: int = 3,
) -> go.Figure:
    """
    Subplot grid showing bar+dot plots for each selected metric.

    Each subplot is one metric (matching plot_metric_bars style).
    Useful for a single publication figure overview.

    Parameters
    ----------
    metrics_df   : output of calculate_growth_metrics()
    selected_metrics : list of metric column names to include
    group_col    : group column
    group_order  : explicit group display order
    palette      : colour palette name
    show_annotations : draw significance stars
    stats_df     : output of run_statistical_tests()
    ncols        : number of subplot columns (default 3)

    Returns
    -------
    go.Figure
    """
    n        = len(selected_metrics)
    nrows    = int(np.ceil(n / ncols))
    active_palette = PALETTES.get(palette, PALETTES["Default"])

    available_groups = metrics_df[group_col].unique().tolist()
    if group_order:
        groups = [g for g in group_order if g in available_groups]
    else:
        groups = available_groups

    fig = make_subplots(
        rows=nrows, cols=ncols,
        subplot_titles=[m.replace("_", " ") for m in selected_metrics],
        horizontal_spacing=0.10,
        vertical_spacing=0.16,
    )

    for idx, metric in enumerate(selected_metrics):
        row = idx // ncols + 1
        col = idx % ncols + 1

        df_m = metrics_df.dropna(subset=[metric])
        summary = (
            df_m.groupby(group_col)[metric]
            .agg(mean_val="mean", sd_val="std")
            .reindex(groups)
            .reset_index()
        )

        for i, group in enumerate(groups):
            color = active_palette[i % len(active_palette)]
            r     = summary[summary[group_col] == group]
            if r.empty:
                continue
            mean_val = float(r["mean_val"].iloc[0])
            sd_val   = float(r["sd_val"].iloc[0]) if not pd.isna(r["sd_val"].iloc[0]) else 0

            show_leg = (idx == 0)   # show legend only for first subplot

            fig.add_trace(go.Bar(
                name=group,
                x=[group],
                y=[mean_val],
                error_y=dict(type="data", array=[sd_val], visible=True,
                             color="#444444", thickness=1.5, width=6),
                marker=dict(color=color, opacity=0.82,
                            line=dict(color="#333333", width=0.8)),
                legendgroup=group,
                showlegend=show_leg,
                hovertemplate=f"<b>{group}</b><br>{metric}: %{{y:.4f}}<extra></extra>",
            ), row=row, col=col)

            # Dots
            bio_vals = df_m[df_m[group_col] == group][metric].dropna().tolist()
            if bio_vals:
                fig.add_trace(go.Scatter(
                    x=[group] * len(bio_vals),
                    y=bio_vals,
                    mode="markers",
                    name=group,
                    legendgroup=group,
                    showlegend=False,
                    marker=dict(color=color, size=7, opacity=0.9,
                                line=dict(color="white", width=1)),
                    hovertemplate=f"<b>{group}</b><br>{metric}: %{{y:.4f}}<extra></extra>",
                ), row=row, col=col)

    fig.update_layout(
        height=320 * nrows,
        barmode="group",
        font=dict(family=_FONT, size=11, color="#2C2C2C"),
        plot_bgcolor="#FFFFFF",
        paper_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=60, r=30, t=60, b=40),
        legend=dict(bgcolor="rgba(255,255,255,0.9)", bordercolor="#CCCCCC",
                    borderwidth=1, font=dict(family=_FONT, size=11)),
        title=_title_style("Growth Metrics Overview"),
    )
    fig.update_xaxes(showgrid=False, linecolor="#CCCCCC",
                     tickfont=dict(family=_FONT, size=10))
    fig.update_yaxes(showgrid=True, gridcolor="#EEEEEE",
                     linecolor="#CCCCCC", tickfont=dict(family=_FONT, size=10))
    return fig


# =============================================================================
# Figure 4 — Statistics results table
# =============================================================================

def plot_stats_table(
    stats_df: pd.DataFrame,
    title: str = "Statistical Test Results",
) -> go.Figure:
    """
    Render the statistics results table as a Plotly Table figure.
    Colour coding: green = p<0.01, yellow = p<0.05, white = ns.
    """
    if stats_df is None or stats_df.empty:
        return _empty_figure("No statistical results available.")

    fill_colors = []
    for _, row in stats_df.iterrows():
        sig = row.get("Significance", "ns")
        if sig in ("****", "***"):
            fill_colors.append("#D4EDDA")
        elif sig in ("**", "*"):
            fill_colors.append("#FFF3CD")
        else:
            fill_colors.append("#FFFFFF")

    def fmt_p(v):
        try:
            f = float(v)
            return "< 0.0001" if f < 0.0001 else f"{f:.4f}"
        except Exception:
            return str(v)

    display = stats_df.copy()
    for col in ["p_value", "p_adj"]:
        if col in display.columns:
            display[col] = display[col].apply(fmt_p)

    fig = go.Figure(go.Table(
        header=dict(
            values=[f"<b>{c}</b>" for c in display.columns],
            fill_color="#2E8B8B",
            font=dict(color="white", family=_FONT, size=12),
            align="left",
            height=32,
        ),
        cells=dict(
            values=[display[c].tolist() for c in display.columns],
            fill_color=[fill_colors] * len(display.columns),
            font=dict(family=_FONT, size=11, color="#2C2C2C"),
            align="left",
            height=28,
        ),
    ))
    fig.update_layout(
        title=_title_style(title),
        margin=dict(l=20, r=20, t=60, b=20),
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


# =============================================================================
# Kept from original stub — now implemented
# =============================================================================

def plot_growth_curves_legacy(df: pd.DataFrame, fitted_df: pd.DataFrame = None) -> go.Figure:
    """Legacy alias — delegates to plot_growth_curves."""
    return plot_growth_curves(df)


def plot_growth_parameters(summary_df: pd.DataFrame) -> go.Figure:
    """Legacy alias — delegates to plot_stats_table."""
    return plot_stats_table(summary_df)


# =============================================================================
# Internal helpers
# =============================================================================

def _hex_to_rgba(hex_color: str, alpha: float) -> str:
    """Convert '#RRGGBB' to 'rgba(r,g,b,alpha)' string."""
    hex_color = hex_color.lstrip("#")
    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _empty_figure(message: str) -> go.Figure:
    """Return a blank figure with a centred message."""
    fig = go.Figure()
    fig.add_annotation(
        text=message, showarrow=False, x=0.5, y=0.5,
        font=dict(size=14, color="#888888"),
    )
    fig.update_layout(
        plot_bgcolor="#FFFFFF", paper_bgcolor="rgba(0,0,0,0)",
        xaxis=dict(visible=False), yaxis=dict(visible=False),
    )
    return fig
