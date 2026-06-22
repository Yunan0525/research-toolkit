# =============================================================================
# utils/export.py  –  Download Button Helpers
# =============================================================================
# PURPOSE
#   Generates Streamlit download buttons for DataFrames and Plotly figures.
#   Each tool page calls these functions to let users export results.
#
# EXPORTS SUPPORTED
#   download_dataframe_csv()       – DataFrame → UTF-8 CSV
#   download_dataframe_excel()     – DataFrame → formatted .xlsx
#   download_figure_png()          – Plotly Figure → 300 dpi PNG (interactive display)
#   download_figure_png_600dpi()   – Plotly Figure → 600 dpi PNG via matplotlib
#                                    (publication-quality; uses correct physical size)
#   download_figure_svg()          – Plotly Figure → SVG vector
#
# 600 DPI APPROACH
#   Plotly's kaleido exporter works in CSS pixels (72 px = 1 in by convention).
#   To reach 600 dpi we would need scale=8.33 — kaleido is unreliable at such
#   high scale factors and text rendering degrades.
#   Instead, download_figure_png_600dpi() re-draws the figure with matplotlib
#   at the target physical size and dpi=600. This produces a clean, correctly
#   sized output file that satisfies journal submission requirements.
#   The interactive Plotly chart in Streamlit is unaffected.
# =============================================================================

import sys
from pathlib import Path

# ── Path bootstrap (export.py lives in utils/, one level below project root) ──
_PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import io
from datetime import datetime

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go

import matplotlib
matplotlib.use("Agg")          # non-interactive backend; safe in Streamlit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import AutoMinorLocator


# =============================================================================
# Timestamp helper
# =============================================================================

def _ts() -> str:
    """Return a compact ISO-style timestamp for filenames: YYYYMMDD_HHMMSS."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


# =============================================================================
# DataFrame exports
# =============================================================================

def download_dataframe_csv(
    df: pd.DataFrame,
    filename: str = "results.csv",
    label: str = "⬇ Download CSV",
) -> None:
    """Render a download button that exports df as UTF-8 CSV."""
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=label,
        data=csv_bytes,
        file_name=filename,
        mime="text/csv",
    )


def download_dataframe_excel(
    df: pd.DataFrame,
    filename: str = "results.xlsx",
    sheet_name: str = "Results",
    label: str = "⬇ Download Excel",
) -> None:
    """Render a download button that exports df as a formatted .xlsx file."""
    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="xlsxwriter") as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
        workbook  = writer.book
        worksheet = writer.sheets[sheet_name]
        header_fmt = workbook.add_format(
            {"bold": True, "bg_color": "#2E8B8B", "font_color": "white"}
        )
        for col_num, col_name in enumerate(df.columns):
            worksheet.write(0, col_num, col_name, header_fmt)
            max_width = max(
                len(str(col_name)),
                df[col_name].astype(str).str.len().max(),
            )
            worksheet.set_column(col_num, col_num, min(max_width + 2, 40))

    st.download_button(
        label=label,
        data=buffer.getvalue(),
        file_name=filename,
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )


# =============================================================================
# Plotly figure exports  (300 dpi via kaleido)
# =============================================================================

def download_figure_png(
    fig: go.Figure,
    filename: str = "figure.png",
    label: str = "⬇ Download PNG (300 dpi)",
    width: int = 1200,
    height: int = 700,
    scale: int = 3,
) -> None:
    """
    Render a download button exporting a Plotly figure as a 300 dpi PNG.

    Uses kaleido (Plotly's static image renderer).
    scale=3 at width=1200 px gives ~300 dpi for a typical journal figure width.

    Parameters
    ----------
    fig    : Plotly Figure
    filename : suggested download filename (without timestamp)
    label  : button text
    width  : base width in CSS pixels
    height : base height in CSS pixels
    scale  : resolution multiplier (3 ≈ 300 dpi)
    """
    # Insert timestamp before the extension
    stem, _, ext = filename.rpartition(".")
    stamped = f"{stem}_{_ts()}.{ext}"

    try:
        png_bytes = fig.to_image(format="png", width=width, height=height, scale=scale)
        st.download_button(
            label=label,
            data=png_bytes,
            file_name=stamped,
            mime="image/png",
        )
    except Exception as exc:
        st.warning(f"⚠️ PNG export unavailable: {exc}. Try the SVG download instead.")


def download_figure_svg(
    fig: go.Figure,
    filename: str = "figure.svg",
    label: str = "⬇ Download SVG (vector)",
    width: int = 1200,
    height: int = 700,
) -> None:
    """
    Render a download button exporting a Plotly figure as an SVG vector file.

    SVG is resolution-independent and can be edited in Illustrator / Inkscape.
    """
    stem, _, ext = filename.rpartition(".")
    stamped = f"{stem}_{_ts()}.{ext}"

    try:
        svg_bytes = fig.to_image(format="svg", width=width, height=height)
        st.download_button(
            label=label,
            data=svg_bytes,
            file_name=stamped,
            mime="image/svg+xml",
        )
    except Exception as exc:
        st.warning(f"⚠️ SVG export unavailable: {exc}.")


# =============================================================================
# 600 dpi matplotlib export
# =============================================================================

def download_figure_png_600dpi(
    fig: go.Figure,
    filename: str = "figure_600dpi.png",
    label: str = "⬇ Download PNG (600 dpi, publication)",
    width_inches: float = 7.0,
    height_inches: float = 4.5,
) -> None:
    """
    Export a Plotly figure as a 600 dpi PNG using matplotlib for rendering.

    Why matplotlib instead of scaling Plotly/kaleido?
    ──────────────────────────────────────────────────
    Kaleido exports at ~72 dpi (CSS pixels). Reaching 600 dpi requires
    scale ≈ 8.3, at which point text anti-aliasing and thin lines degrade.
    Matplotlib renders vector elements at any dpi with clean results.

    This function re-draws the figure data (bars, error bars, dots, reference
    line, annotations) using matplotlib, matching the Plotly visual style as
    closely as possible: same colours, same layout, same axis labels.

    The Plotly interactive chart in Streamlit is not affected.

    Parameters
    ----------
    fig          : Plotly Figure produced by plot_fold_change_bar() or
                   plot_log2fc_bar(). Heatmaps and table figures are passed
                   through to a fallback kaleido export.
    filename     : suggested download filename (timestamp auto-appended)
    label        : button text
    width_inches : output figure width in inches (7 in = typical journal column)
    height_inches: output figure height in inches
    """
    stem, _, ext = filename.rpartition(".")
    stamped = f"{stem}_{_ts()}.{ext}"

    png_bytes = _plotly_to_matplotlib_600dpi(fig, width_inches, height_inches)

    st.download_button(
        label=label,
        data=png_bytes,
        file_name=stamped,
        mime="image/png",
    )


def _plotly_to_matplotlib_600dpi(
    fig: go.Figure,
    width_inches: float = 7.0,
    height_inches: float = 4.5,
) -> bytes:
    """
    Convert a Plotly bar+dot figure to a 600 dpi matplotlib PNG.

    Reads bar, scatter, and hline traces from the Plotly Figure object,
    reconstructs them in matplotlib, and saves at 600 dpi.

    Returns
    -------
    bytes : PNG image bytes suitable for st.download_button(data=...).
    """
    # ── Extract traces from the Plotly figure ─────────────────────────────────
    bar_traces    = [t for t in fig.data if isinstance(t, go.Bar)]
    scatter_traces= [t for t in fig.data if isinstance(t, go.Scatter)]

    # Fall back to kaleido for figures without bar traces (heatmaps, tables)
    if not bar_traces:
        return _kaleido_fallback(fig, width_inches, height_inches)

    # ── Reconstruct gene list and group list from bar traces ──────────────────
    genes  = list(bar_traces[0].x)   # categorical x values from first bar
    groups = [t.name for t in bar_traces]
    n_genes  = len(genes)
    n_groups = len(groups)

    # Extract colour from each bar trace
    colors = [t.marker.color for t in bar_traces]

    # ── Compute bar geometry (must match _bar_dot_plot constants exactly) ──────
    BARGAP      = 0.25
    BARGROUPGAP = 0.08
    bar_width = (1 - BARGAP) / n_groups * (1 - BARGROUPGAP)
    offsets   = [
        (i - (n_groups - 1) / 2) * (1 - BARGAP) / n_groups
        for i in range(n_groups)
    ]
    jitter_max = bar_width * 0.18

    gene_positions = np.arange(n_genes, dtype=float)

    # ── Extract reference line y-value from hline shapes ─────────────────────
    ref_y = None
    if hasattr(fig, "layout") and fig.layout.shapes:
        for shape in fig.layout.shapes:
            if getattr(shape, "type", None) == "line":
                ref_y = getattr(shape, "y0", None)
                break
    # Fallback: infer from y-axis title
    if ref_y is None:
        ytitle = ""
        try:
            ytitle = fig.layout.yaxis.title.text or ""
        except Exception:
            pass
        ref_y = 0.0 if "log" in ytitle.lower() else 1.0

    # ── Extract axis labels ───────────────────────────────────────────────────
    try:
        y_label = fig.layout.yaxis.title.text or ""
    except Exception:
        y_label = ""
    try:
        plot_title = fig.layout.title.text or ""
        # Strip HTML superscripts for matplotlib
        import re
        plot_title = re.sub(r"<[^>]+>", "", plot_title)
    except Exception:
        plot_title = ""

    # ── Build matplotlib figure ───────────────────────────────────────────────
    mpl_fig, ax = plt.subplots(figsize=(width_inches, height_inches), dpi=600)

    for i, bar_trace in enumerate(bar_traces):
        color  = colors[i]
        group  = groups[i]
        y_mean = list(bar_trace.y)
        y_err  = list(bar_trace.error_y.array) if bar_trace.error_y and bar_trace.error_y.array else [0]*n_genes

        bar_centers = gene_positions + offsets[i]

        # Bars with error bars
        ax.bar(
            bar_centers,
            y_mean,
            width=bar_width,
            color=color,
            alpha=0.82,
            edgecolor="#333333",
            linewidth=0.7,
            label=group,
            yerr=y_err,
            capsize=4,
            error_kw=dict(ecolor="#444444", elinewidth=1.6, capthick=1.6),
        )

        # ── Individual biological replicate dots from scatter traces ──────────
        # Match scatter traces to this group by legendgroup name
        grp_scatters = [
            t for t in scatter_traces
            if getattr(t, "legendgroup", None) == group
        ]

        for scat in grp_scatters:
            # scat.x is the gene name (categorical), scat.y is the y-value
            x_cats = list(scat.x)
            y_vals = list(scat.y)

            for x_cat, y_val in zip(x_cats, y_vals):
                if x_cat not in genes:
                    continue
                gene_idx = genes.index(x_cat)
                # Spread dots evenly across the bar
                # Count how many dots belong to this gene+group
                same_gene_y = [
                    y_v for x_c, y_v in zip(x_cats, y_vals) if x_c == x_cat
                ]
                n_pts = len(same_gene_y)
                jitter_vals = (
                    np.linspace(-jitter_max, jitter_max, n_pts)
                    if n_pts > 1 else np.array([0.0])
                )
                dot_idx = same_gene_y.index(y_val)
                x_dot = gene_positions[gene_idx] + offsets[i] + jitter_vals[dot_idx]

                ax.scatter(
                    x_dot, y_val,
                    color=color,
                    s=28,
                    zorder=5,
                    edgecolors="white",
                    linewidths=0.8,
                )

    # ── Reference line ────────────────────────────────────────────────────────
    ax.axhline(ref_y, color="#AAAAAA", linestyle="--", linewidth=0.9, zorder=1)

    # ── Significance annotations (from Plotly annotations) ───────────────────
    try:
        for annot in fig.layout.annotations:
            text = getattr(annot, "text", "")
            x    = getattr(annot, "x", None)
            y    = getattr(annot, "y", None)
            if text and x in genes and y is not None:
                ax.text(
                    genes.index(x), float(y), text,
                    ha="center", va="bottom",
                    fontsize=12, fontfamily="sans-serif",
                    color="#222222",
                )
    except Exception:
        pass

    # ── Axes and style ────────────────────────────────────────────────────────
    ax.set_xticks(gene_positions)
    ax.set_xticklabels(genes, fontsize=11, fontfamily="sans-serif")
    ax.tick_params(axis="y", labelsize=10)
    ax.set_xlabel("Gene", fontsize=13, fontfamily="sans-serif", labelpad=6)
    ax.set_ylabel(y_label, fontsize=13, fontfamily="sans-serif", labelpad=6)
    if plot_title:
        ax.set_title(plot_title, fontsize=13, fontfamily="sans-serif",
                     loc="left", pad=10, color="#1F4E79")

    ax.set_facecolor("white")
    ax.grid(axis="y", color="#EEEEEE", linewidth=0.8, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")
    ax.tick_params(direction="out", colors="#2C2C2C")

    legend = ax.legend(
        fontsize=10,
        framealpha=0.9,
        edgecolor="#CCCCCC",
        fancybox=False,
    )
    for text in legend.get_texts():
        text.set_fontfamily("sans-serif")

    mpl_fig.patch.set_facecolor("white")
    plt.tight_layout()

    # ── Render to bytes ───────────────────────────────────────────────────────
    buf = io.BytesIO()
    mpl_fig.savefig(buf, format="png", dpi=600, bbox_inches="tight",
                    facecolor="white")
    buf.seek(0)
    png_bytes = buf.read()
    plt.close(mpl_fig)
    return png_bytes


def _kaleido_fallback(
    fig: go.Figure,
    width_inches: float,
    height_inches: float,
) -> bytes:
    """
    Fallback for figure types that can't be re-drawn by matplotlib
    (heatmaps, stats tables).  Uses kaleido at the highest safe scale.
    """
    try:
        width_px  = int(width_inches * 72)
        height_px = int(height_inches * 72)
        return fig.to_image(
            format="png",
            width=width_px,
            height=height_px,
            scale=6,   # ~432 dpi — best kaleido can reliably do
        )
    except Exception as exc:
        # Last resort: return a tiny placeholder
        mpl_fig, ax = plt.subplots(figsize=(4, 2))
        ax.text(0.5, 0.5, f"Export failed:\n{exc}", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="red")
        buf = io.BytesIO()
        mpl_fig.savefig(buf, format="png", dpi=150)
        buf.seek(0)
        plt.close(mpl_fig)
        return buf.read()
