# =============================================================================
# visualization/qpcr_plots.py  –  qPCR Plotly Figure Builders
# =============================================================================
# PURPOSE
#   Build publication-quality interactive Plotly figures from DataFrames
#   produced by analysis/qpcr.py.  No Streamlit code; no analysis math.
#
# INPUT CONTRACT
#   Every function receives data that has ALREADY had technical replicates
#   averaged.  Each dot = ONE biological replicate.
#
# DOT POSITIONING (grouped bars)
# ─────────────────────────────────────────────────────────────────────────────
#   Plotly's grouped bar chart uses a categorical x-axis.  The physical pixel
#   position of each group's bar within a gene category is determined by
#   bargap and bargroupgap.  We cannot attach scatter dots to a categorical
#   axis and expect them to land on the correct bar.
#
#   Solution: overlay a NUMERIC x-axis (xaxis2) on top of the categorical
#   axis (xaxis).  Gene categories map to integer positions 0, 1, 2, …
#   Each group's bar center is computed with the formula:
#
#       bar_width = (1 − BARGAP) / N_groups × (1 − BARGROUPGAP)
#       offset_i  = (i − (N_groups−1)/2) × (1−BARGAP)/N_groups
#       bar_center = gene_index + offset_i
#
#   Dots are placed at bar_center ± small_jitter (≤18% of bar_width).
#   This formula is kept in sync with the layout constants BARGAP / BARGROUPGAP.
#
# FIGURES
#   plot_fold_change_bar  – bars mean FC ± SD, dots = bio replicates
#   plot_log2fc_bar       – same on log₂ scale
#   plot_deltact_heatmap  – ΔCt matrix: biological samples × genes
#   plot_stats_results    – statistics table as Plotly Table
#
# STYLE (Nature-journal conventions)
#   Font   : Arial, sans-serif  (14 pt title, 13 pt axis, 11 pt ticks)
#   Errors : SD of biological replicates, cap width 8 px
#   Grid   : horizontal only, #EEEEEE, 1 px
#   BG     : white plot area, transparent paper
# =============================================================================

import plotly.graph_objects as go
import pandas as pd
import numpy as np


# =============================================================================
# Shared style constants
# =============================================================================

# =============================================================================
# Color palettes — selectable via the UI (palette= parameter on plot functions)
# =============================================================================

PALETTES = {
    "Default": [
        "#2E8B8B",   # Teal
        "#E07B6A",   # Coral
        "#4C72B0",   # Steel blue
        "#DD8452",   # Amber
        "#55A868",   # Sage green
        "#C44E52",   # Crimson
        "#8172B3",   # Purple
        "#937860",   # Taupe
    ],
    "Nature-style": [
        "#E64B35",   # Red         — Nature red
        "#4DBBD5",   # Cyan        — Nature blue
        "#00A087",   # Green       — Nature teal
        "#3C5488",   # Navy        — Nature navy
        "#F39B7F",   # Salmon
        "#8491B4",   # Lavender blue
        "#91D1C2",   # Mint
        "#DC0000",   # Bright red
    ],
    "Colorblind-friendly": [
        "#0072B2",   # Blue        — Wong palette
        "#E69F00",   # Orange
        "#009E73",   # Green
        "#CC79A7",   # Pink/Purple
        "#56B4E9",   # Sky blue
        "#D55E00",   # Vermillion
        "#F0E442",   # Yellow
        "#000000",   # Black
    ],
    "Pastel": [
        "#AEC6CF",   # Pastel blue
        "#FFD1DC",   # Pastel pink
        "#B5EAD7",   # Pastel mint
        "#FFDAC1",   # Pastel peach
        "#C7CEEA",   # Pastel lavender
        "#E2F0CB",   # Pastel lime
        "#F8C8D4",   # Pastel rose
        "#BEE3F8",   # Pastel sky
    ],
    "High contrast": [
        "#000000",   # Black
        "#E6194B",   # Red
        "#3CB44B",   # Green
        "#4363D8",   # Blue
        "#F58231",   # Orange
        "#911EB4",   # Purple
        "#42D4F4",   # Cyan
        "#F032E6",   # Magenta
    ],
    "Viridis": [
        "#440154",   # Deep purple
        "#31688E",   # Blue
        "#35B779",   # Green
        "#FDE725",   # Yellow
        "#21908C",   # Teal
        "#5DC863",   # Light green
        "#443983",   # Indigo
        "#90D743",   # Lime
    ],
}

# Keep PALETTE as the default for any internal fallback
PALETTE = PALETTES["Default"]

_FONT = "Arial, sans-serif"

# Bar geometry constants — MUST stay in sync with _compute_bar_geometry()
# and with the matplotlib renderer in utils/export.py (_plotly_to_matplotlib_600dpi).
BARGAP      = 0.25   # gap between gene clusters (fraction of slot width)
BARGROUPGAP = 0.08   # gap between bars within a cluster

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
    margin=dict(l=70, r=40, t=80, b=80),
    legend=dict(
        bgcolor="rgba(255,255,255,0.9)",
        bordercolor="#CCCCCC",
        borderwidth=1,
        font=dict(family=_FONT, size=11),
    ),
)


def _title_style(text: str) -> dict:
    return dict(
        text=text,
        font=dict(family=_FONT, size=14, color="#1F4E79"),
        x=0.05,
    )


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
# Bar geometry helper
# =============================================================================

def _compute_bar_geometry(n_groups: int) -> tuple[list[float], float, float]:
    """
    Compute the x-axis offset and width of each group's bar.

    Returns
    -------
    offsets : list[float]
        x-offset from the gene's integer center for each group (0-indexed).
        e.g. for 2 groups: [-0.1875, +0.1875]
    bar_width : float
        Width of each individual bar in numeric axis units.
    jitter_max : float
        Maximum horizontal jitter applied to dots (18 % of bar_width).
        Dots span ≤36 % of the bar width, staying visually inside the bar.
    """
    bar_width  = (1 - BARGAP) / n_groups * (1 - BARGROUPGAP)
    offsets    = [
        (i - (n_groups - 1) / 2) * (1 - BARGAP) / n_groups
        for i in range(n_groups)
    ]
    jitter_max = bar_width * 0.18
    return offsets, bar_width, jitter_max


# =============================================================================
# Shared bar + dot plot builder
# =============================================================================

def _bar_dot_plot(
    summary_df: pd.DataFrame,
    detail_df: pd.DataFrame,
    group_col: str,
    gene_col: str,
    y_mean: str,
    y_sd: str,
    y_detail: str,
    y_axis_title: str,
    ref_line: float,
    ref_label: str,
    title: str,
    stats_df: pd.DataFrame | None = None,
    group_order: list | None = None,
    palette: str = "Default",
    show_annotations: bool = True,
) -> go.Figure:
    """
    Shared builder: grouped bar chart + individual dot overlay.

    Bars  — one per group per gene; height = group mean; error = SD.
    Dots  — one per biological replicate; positioned ON the corresponding bar
            using a numeric overlay axis (xaxis2), not on the categorical center.
    Stats — significance stars drawn above bar pairs (2-group only).

    Parameters
    ----------
    summary_df  : output of summarise_results()
    detail_df   : output of calculate_delta_delta_ct()
    group_col   : name of group column in detail_df
    gene_col    : name of gene column in detail_df
    y_mean      : column name for group mean in summary_df
    y_sd        : column name for SD in summary_df
    y_detail    : column name for individual values in detail_df
    y_axis_title: label for the y-axis
    ref_line    : y-value for the horizontal reference line (1.0 or 0.0)
    ref_label   : annotation text for the reference line
    title       : figure title (HTML allowed)
    stats_df         : optional statistics table for significance annotations
    group_order      : list of group names in desired display order; if None,
                       uses the order found in summary_df. Groups not present
                       in summary_df are silently ignored.
    palette          : name of color palette from PALETTES dict (default "Default")
    show_annotations : if False, significance stars are not drawn even when
                       stats_df is provided
    """
    _validate(summary_df, ["Gene", "Group", y_mean, y_sd])
    _validate(detail_df,  [gene_col, group_col, y_detail])

    genes  = summary_df["Gene"].unique().tolist()

    # ── Resolve group order ───────────────────────────────────────────────────
    # If caller supplies an explicit order, use it (filtering to groups that
    # actually exist in the data).  Otherwise fall back to natural order.
    available_groups = summary_df["Group"].unique().tolist()
    if group_order is not None:
        groups = [g for g in group_order if g in available_groups]
        if not groups:                          # safety: nothing matched
            groups = available_groups
    else:
        groups = available_groups

    # ── Resolve color palette ─────────────────────────────────────────────────
    active_palette = PALETTES.get(palette, PALETTES["Default"])

    n_genes  = len(genes)
    n_groups = len(groups)

    offsets, bar_width, jitter_max = _compute_bar_geometry(n_groups)

    # Map gene names to integer positions: IL6→0, TNF→1, …
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    fig = go.Figure()

    # ── Categorical x-axis bars ───────────────────────────────────────────────
    for i, group in enumerate(groups):
        color   = active_palette[i % len(active_palette)]
        grp_sum = summary_df[summary_df["Group"] == group]

        # Reorder to match `genes` list so bars align with dot overlay axis
        grp_sum = grp_sum.set_index("Gene").reindex(genes).reset_index()

        fig.add_trace(
            go.Bar(
                name=group,
                x=grp_sum["Gene"],
                y=grp_sum[y_mean],
                xaxis="x",       # categorical axis (gene names)
                error_y=dict(
                    type="data",
                    array=grp_sum[y_sd].tolist(),
                    visible=True,
                    color="#444444",
                    thickness=1.8,
                    width=8,
                ),
                marker=dict(
                    color=color,
                    opacity=0.82,
                    line=dict(color="#333333", width=0.8),
                ),
                legendgroup=group,
                showlegend=True,
                hovertemplate=(
                    "<b>%{x}</b><br>"
                    f"Group: {group}<br>"
                    f"{y_axis_title}: %{{y:.3f}} ± SD<br>"
                    "<extra></extra>"
                ),
            )
        )

    # ── Numeric overlay axis dots (xaxis2) ────────────────────────────────────
    # xaxis2 shares the same plot area but has a numeric domain [−0.5, n_genes−0.5]
    # so that integer 0 aligns with the centre of the first gene cluster,
    # integer 1 with the second, etc.  This is guaranteed by setting
    # xaxis2.range to [−0.5, n_genes − 0.5].

    for i, group in enumerate(groups):
        color   = active_palette[i % len(active_palette)]
        grp_det = detail_df[detail_df[group_col] == group]

        # Collect all (numeric_x, y_value) pairs for this group's dots
        all_x = []
        all_y = []
        all_hover_gene = []

        for gene in genes:
            gene_det = grp_det[grp_det[gene_col] == gene]
            if gene_det.empty:
                continue

            y_vals = gene_det[y_detail].tolist()
            n_pts  = len(y_vals)
            gene_idx = gene_to_idx[gene]

            # Deterministic symmetric jitter centred on the bar x-position
            if n_pts == 1:
                jitter = np.array([0.0])
            else:
                jitter = np.linspace(-jitter_max, jitter_max, n_pts)

            x_positions = gene_idx + offsets[i] + jitter

            all_x.extend(x_positions.tolist())
            all_y.extend(y_vals)
            all_hover_gene.extend([gene] * n_pts)

        if not all_x:
            continue

        # Build per-point hover text
        hover_texts = [
            f"<b>{g}</b><br>Group: {group}<br>"
            f"{y_axis_title}: {y:.3f}<br>Biological replicate<extra></extra>"
            for g, y in zip(all_hover_gene, all_y)
        ]

        fig.add_trace(
            go.Scatter(
                x=all_x,
                y=all_y,
                mode="markers",
                name=group,
                xaxis="x2",      # numeric overlay axis
                legendgroup=group,
                showlegend=False,  # dots share bar legend entry
                marker=dict(
                    color=color,
                    size=8,
                    opacity=0.95,
                    line=dict(color="white", width=1.2),
                    symbol="circle",
                ),
                hovertemplate=hover_texts,
            )
        )

    # ── Reference line ────────────────────────────────────────────────────────
    fig.add_hline(
        y=ref_line,
        line_dash="dash",
        line_color="#AAAAAA",
        line_width=1,
        annotation_text=ref_label,
        annotation_position="right",
        annotation_font=dict(size=10, color="#AAAAAA"),
    )

    # ── Significance annotations ──────────────────────────────────────────────
    # Only drawn when the user has enabled them AND stats results are available.
    if show_annotations and stats_df is not None and not stats_df.empty:
        _add_significance_annotations(fig, stats_df, summary_df, groups, y_mean, y_sd)

    # ── Layout with dual x-axes ───────────────────────────────────────────────
    fig.update_layout(
        barmode="group",
        bargap=BARGAP,
        bargroupgap=BARGROUPGAP,
        # Primary categorical axis for bars
        xaxis=dict(
            title="Gene",
            **_AXIS,
        ),
        # Numeric overlay axis for dots — same physical space, hidden labels
        xaxis2=dict(
            overlaying="x",
            side="bottom",
            range=[-0.5, n_genes - 0.5],   # aligns integer 0 with first gene cluster
            visible=False,                  # hides ticks/labels; bars already label x
            fixedrange=True,
        ),
        yaxis=_y_axis_style(y_axis_title),
        title=_title_style(title),
        **_LAYOUT,
    )
    return fig


# =============================================================================
# Significance annotation helper
# =============================================================================

def _add_significance_annotations(
    fig: go.Figure,
    stats_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    groups: list,
    y_mean: str,
    y_sd: str,
) -> None:
    """
    Add significance star annotations above bars (2-group comparisons only).

    For ≥3 groups, stars would require pairwise brackets which clutter the
    figure; the Statistics tab is the primary output in that case.
    """
    if len(groups) != 2:
        return

    for _, row in stats_df.iterrows():
        gene  = row.get("Gene", "")
        stars = row.get("Significance", "ns")
        if stars in ("ns", "n/a") or pd.isna(row.get("p_adj", np.nan)):
            continue

        gene_sum = summary_df[summary_df["Gene"] == gene]
        if gene_sum.empty:
            continue

        # Place star above the tallest bar top (mean + SD)
        tops = gene_sum[y_mean] + gene_sum[y_sd].fillna(0)
        max_top = tops.max()
        y_annot = max_top * 1.12 if max_top > 0 else abs(max_top) * 0.12 + 0.3

        fig.add_annotation(
            x=gene,
            y=y_annot,
            text=stars,
            showarrow=False,
            font=dict(family=_FONT, size=14, color="#222222"),
            xanchor="center",
        )


# =============================================================================
# Public plot functions
# =============================================================================

def plot_fold_change_bar(
    summary_df: pd.DataFrame,
    detail_df: pd.DataFrame,
    group_col: str,
    gene_col: str,
    stats_df: pd.DataFrame | None = None,
    title: str = "Fold Change (2<sup>−ΔΔCt</sup>)  —  mean ± SD, biological replicates",
    group_order: list | None = None,
    palette: str = "Default",
    show_annotations: bool = True,
) -> go.Figure:
    """
    Bar + dot plot of Fold Change.

    Bars : mean FC ± SD (biological replicates).
    Dots : individual biological replicate FC values, aligned to their bar.

    Parameters
    ----------
    summary_df       : output of summarise_results()
                       Required columns: Gene, Group, Mean_FC, SD_FC
    detail_df        : output of calculate_delta_delta_ct()
                       Required columns: gene_col, group_col, Fold_Change
    group_col        : group column name in detail_df
    gene_col         : gene column name in detail_df
    stats_df         : optional output of run_statistical_tests()
                       If provided, significance stars added above bars (2 groups only)
    title            : figure title (HTML ok)
    group_order      : explicit group display order (list of group name strings)
    palette          : color palette name — one of PALETTES keys
    show_annotations : whether to draw significance stars
    """
    return _bar_dot_plot(
        summary_df       = summary_df,
        detail_df        = detail_df,
        group_col        = group_col,
        gene_col         = gene_col,
        y_mean           = "Mean_FC",
        y_sd             = "SD_FC",
        y_detail         = "Fold_Change",
        y_axis_title     = "Fold Change (mean ± SD)",
        ref_line         = 1.0,
        ref_label        = "FC = 1",
        title            = title,
        stats_df         = stats_df,
        group_order      = group_order,
        palette          = palette,
        show_annotations = show_annotations,
    )


def plot_log2fc_bar(
    summary_df: pd.DataFrame,
    detail_df: pd.DataFrame,
    group_col: str,
    gene_col: str,
    stats_df: pd.DataFrame | None = None,
    title: str = "log₂ Fold Change (−ΔΔCt)  —  mean ± SD, biological replicates",
    group_order: list | None = None,
    palette: str = "Default",
    show_annotations: bool = True,
) -> go.Figure:
    """
    Bar + dot plot of log₂ Fold Change.

    Bars above 0 = upregulated; bars below 0 = downregulated.

    Parameters
    ----------
    summary_df       : Required columns: Gene, Group, Mean_log2FC, SD_log2FC
    detail_df        : Required column: log2FC (plus group_col, gene_col)
    group_col        : group column name in detail_df
    gene_col         : gene column name in detail_df
    stats_df         : optional statistics table
    title            : figure title
    group_order      : explicit group display order (list of group name strings)
    palette          : color palette name — one of PALETTES keys
    show_annotations : whether to draw significance stars
    """
    return _bar_dot_plot(
        summary_df       = summary_df,
        detail_df        = detail_df,
        group_col        = group_col,
        gene_col         = gene_col,
        y_mean           = "Mean_log2FC",
        y_sd             = "SD_log2FC",
        y_detail         = "log2FC",
        y_axis_title     = "log₂ Fold Change (mean ± SD)",
        ref_line         = 0.0,
        ref_label        = "log₂FC = 0",
        title            = title,
        stats_df         = stats_df,
        group_order      = group_order,
        palette          = palette,
        show_annotations = show_annotations,
    )


def plot_deltact_heatmap(
    detail_df: pd.DataFrame,
    sample_col: str,
    gene_col: str,
    title: str = "ΔCt heatmap — biological samples × genes",
) -> go.Figure:
    """
    Heatmap of ΔCt values: rows = biological samples, columns = target genes.

    Each cell = ΔCt for one biological sample (tech reps already averaged).
    Diverging RdBu scale: Blue = lower ΔCt (higher expression); Red = higher ΔCt.
    """
    if "Delta_Ct" not in detail_df.columns:
        raise ValueError("'Delta_Ct' column not found in detail_df.")

    pivot = detail_df.pivot_table(
        index=sample_col,
        columns=gene_col,
        values="Delta_Ct",
        aggfunc="mean",
    )

    fig = go.Figure(
        go.Heatmap(
            z=pivot.values,
            x=pivot.columns.tolist(),
            y=pivot.index.tolist(),
            colorscale="RdBu_r",
            zmid=float(np.nanmean(pivot.values)),
            colorbar=dict(
                title=dict(text="ΔCt", font=dict(family=_FONT, size=12)),
                tickfont=dict(family=_FONT, size=10),
                thickness=15,
                len=0.8,
            ),
            hovertemplate=(
                "Sample: %{y}<br>Gene: %{x}<br>ΔCt: %{z:.3f}<br><extra></extra>"
            ),
        )
    )
    fig.update_layout(
        xaxis=dict(title="Gene", **_AXIS),
        yaxis=dict(title="Biological Sample", autorange="reversed", **_AXIS),
        title=_title_style(title),
        **_LAYOUT,
    )
    return fig


def plot_stats_results(
    stats_df: pd.DataFrame,
    title: str = "Statistical Test Results",
) -> go.Figure:
    """
    Render the statistics results table as a Plotly Table figure.

    Colour coding: green = p<0.01, yellow = p<0.05, white = ns.
    """
    if stats_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No statistical results available.",
            showarrow=False, x=0.5, y=0.5,
            font=dict(size=14, color="#888888"),
        )
        return fig

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

    fig = go.Figure(
        go.Table(
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
        )
    )
    fig.update_layout(
        title=_title_style(title),
        margin=dict(l=20, r=20, t=60, b=20),
        paper_bgcolor="rgba(0,0,0,0)",
    )
    return fig


# =============================================================================
# Internal validation
# =============================================================================

def _validate(df: pd.DataFrame, required: list) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"DataFrame missing required columns: {missing}.")
