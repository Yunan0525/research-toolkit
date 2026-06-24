# =============================================================================
# analysis/growth_curve.py  –  Bacterial Growth Curve Analysis Logic
# =============================================================================
# PURPOSE
#   Pure Python / pandas / numpy / scipy functions for bacterial growth curve
#   analysis from OD600 time-series data.
#   NO Streamlit imports — works identically in notebooks, scripts, or APIs.
#
# MANDATORY PIPELINE (always call in this order)
# ─────────────────────────────────────────────
#   1. parse_wide_format()
#        Input : wide-format DataFrame (Time col + duplicate sample/blank cols)
#        Output: long-format DataFrame with columns:
#                Time_h, Sample, TechRep_Index, OD600, IsBiological
#
#   2. subtract_blank()
#        Input : long-format data + blank_assignment dict
#        Output: same shape DataFrame with Corrected_OD column added
#                Must be called BEFORE average_technical_replicates()
#
#   3. average_technical_replicates()
#        Input : background-corrected long-format data
#        Output: one row per (Sample × Timepoint); mean Corrected_OD,
#                SD, N across technical replicates
#
#   4. calculate_growth_metrics()
#        Input : averaged long-format data (one biological sample per group)
#        Output: one row per Sample with 9 metrics
#
#   5. summarise_metrics()
#        Input : per-biological-sample metrics from step 4
#        Output: mean ± SD per Group across biological replicates
#
#   6. run_statistical_tests()   [optional]
#        Input : per-biological-sample metrics from step 4
#        Output: p-values and test statistics per metric per group pair
#
# DESIGN DECISIONS
# ─────────────────
#   • Background subtraction happens BEFORE tech-rep averaging (per-timepoint,
#     per-well correction preserves within-well noise structure).
#   • Technical replicates are columns with the same name in wide format.
#     They are identified by pandas duplicate-column detection after read_csv.
#   • Statistics operate on biological replicates (one value per sample),
#     never on technical replicates.
#   • All metrics return NaN gracefully when data is insufficient.
# =============================================================================

import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
from scipy.signal import savgol_filter
from statsmodels.stats.multitest import multipletests


# =============================================================================
# STEP 1  —  Parse wide-format data
# =============================================================================

def parse_wide_format(
    df: pd.DataFrame,
    time_col: str = "Time_h",
) -> pd.DataFrame:
    """
    Convert a wide-format OD DataFrame into long format.

    Wide format has one row per timepoint and one column per well.
    Duplicate column names indicate technical replicates of the same sample.

    Parameters
    ----------
    df : pd.DataFrame
        Wide-format data. First column must be time (name given by time_col).
        All other columns are sample wells. Duplicate column names = tech reps.
        pandas read_csv renames duplicates as Name, Name.1, Name.2, … —
        this function handles both the original duplicate names AND the
        pandas-renamed variants.
    time_col : str
        Name of the time column (default "Time_h").

    Returns
    -------
    pd.DataFrame
        Long format with columns:
        • Time_h         – timepoint (numeric)
        • Sample         – sample name (duplicates stripped of .1 .2 suffix)
        • TechRep_Index  – integer 0, 1, 2, … within each sample
        • OD600          – raw OD measurement
    """
    df = df.copy()

    # Ensure time column is numeric
    df[time_col] = pd.to_numeric(df[time_col], errors="coerce")

    # Collect all non-time columns
    sample_cols = [c for c in df.columns if c != time_col]

    rows = []
    # Track which base sample name we've seen, to assign TechRep_Index
    seen_counts: dict = {}

    for col in sample_cols:
        # Strip pandas duplicate suffix (.1 .2 .3 …)
        base_name = col.rsplit(".", 1)[0] if _is_pandas_dup_suffix(col) else col

        tech_idx = seen_counts.get(base_name, 0)
        seen_counts[base_name] = tech_idx + 1

        for _, row in df.iterrows():
            rows.append({
                "Time_h":        row[time_col],
                "Sample":        base_name,
                "TechRep_Index": tech_idx,
                "OD600":         pd.to_numeric(row[col], errors="coerce"),
            })

    long_df = pd.DataFrame(rows)
    long_df = long_df.dropna(subset=["Time_h"])
    long_df = long_df.sort_values(["Sample", "TechRep_Index", "Time_h"]).reset_index(drop=True)
    return long_df


def _is_pandas_dup_suffix(col_name: str) -> bool:
    """
    Return True if col_name ends with .N where N is a positive integer,
    indicating pandas auto-renamed a duplicate column.
    e.g. 'WT.1' → True,  'WT' → False,  'BHI_Blank.2' → True
    """
    parts = col_name.rsplit(".", 1)
    if len(parts) == 2:
        return parts[1].isdigit()
    return False


def get_detected_samples(long_df: pd.DataFrame) -> dict:
    """
    Return a summary dict of detected samples and their tech-rep counts.

    Returns
    -------
    dict
        { sample_name: n_tech_reps, … }
    """
    return (
        long_df.groupby("Sample")["TechRep_Index"]
        .nunique()
        .to_dict()
    )


# =============================================================================
# STEP 2  —  Background subtraction
# =============================================================================

def subtract_blank(
    long_df: pd.DataFrame,
    blank_assignment: dict,
    negative_handling: str = "zero",
) -> pd.DataFrame:
    """
    Subtract blank media OD from each sample well, per timepoint.

    Background correction is performed BEFORE technical-replicate averaging
    so that each well is corrected individually.

    Parameters
    ----------
    long_df : pd.DataFrame
        Output of parse_wide_format(). Must contain Time_h, Sample,
        TechRep_Index, OD600.
    blank_assignment : dict
        Maps each non-blank sample name → blank sample name(s).
        The blank OD used is the mean across all tech reps of the assigned
        blank at each timepoint.
        Example:
            { 'WT': 'BHI_Blank', 'Mutant': 'MRS_Blank' }
        Blank samples themselves are NOT corrected (they remain in the
        output for reference but get Corrected_OD = NaN).
        A sample not present in blank_assignment is left uncorrected
        (Corrected_OD = OD600).
    negative_handling : str
        How to handle Corrected_OD < 0:
        'zero'  – set to 0.0 (default; biologically sensible)
        'keep'  – keep raw negative values
        'na'    – replace with NaN

    Returns
    -------
    pd.DataFrame
        Same rows as long_df with two new columns:
        • 'Corrected_OD'   – background-corrected OD
        • 'Blank_Used'     – name of the blank applied (or 'none')
    """
    df = long_df.copy()
    df["Corrected_OD"] = np.nan
    df["Blank_Used"]   = "none"

    # Pre-compute mean blank OD per (blank_sample, timepoint)
    # across all tech reps of each blank
    blank_means: dict = {}   # blank_name → {timepoint: mean_OD}
    all_blank_names = set(blank_assignment.values())

    for blank_name in all_blank_names:
        blank_rows = df[df["Sample"] == blank_name]
        if blank_rows.empty:
            continue
        bm = (
            blank_rows
            .groupby("Time_h")["OD600"]
            .mean()
            .to_dict()
        )
        blank_means[blank_name] = bm

    # Identify all blank sample names so we can skip correcting them
    all_blank_sample_names = df["Sample"].unique()
    blank_sample_set = {
        s for s in all_blank_sample_names
        if s in all_blank_names
    }

    for idx, row in df.iterrows():
        sample  = row["Sample"]
        time    = row["Time_h"]
        raw_od  = row["OD600"]

        if sample in blank_sample_set:
            # Blank wells: leave Corrected_OD as NaN (they are reference only)
            df.at[idx, "Blank_Used"] = "self"
            continue

        assigned_blank = blank_assignment.get(sample)
        if assigned_blank is None or assigned_blank not in blank_means:
            # No blank assigned: pass through uncorrected
            df.at[idx, "Corrected_OD"] = raw_od
            df.at[idx, "Blank_Used"]   = "none"
            continue

        blank_od = blank_means[assigned_blank].get(time, np.nan)
        corrected = raw_od - blank_od if not np.isnan(blank_od) else raw_od

        # Handle negative values
        if not np.isnan(corrected) and corrected < 0:
            if negative_handling == "zero":
                corrected = 0.0
            elif negative_handling == "na":
                corrected = np.nan
            # "keep": do nothing

        df.at[idx, "Corrected_OD"] = corrected
        df.at[idx, "Blank_Used"]   = assigned_blank

    return df.reset_index(drop=True)


# =============================================================================
# STEP 3  —  Average technical replicates
# =============================================================================

def average_technical_replicates(
    corrected_df: pd.DataFrame,
    sample_col: str = "Sample",
    group_assignment: dict | None = None,
) -> pd.DataFrame:
    """
    Collapse technical replicates into one mean Corrected_OD per
    biological sample × timepoint.

    Must be called AFTER subtract_blank() so that each well has already
    been individually background-corrected.

    Parameters
    ----------
    corrected_df : pd.DataFrame
        Output of subtract_blank(). Must contain Time_h, Sample,
        TechRep_Index, Corrected_OD.
    sample_col : str
        Column identifying sample names (default 'Sample').
    group_assignment : dict or None
        Maps sample name → group name (e.g. {'WT_R1': 'WT', 'WT_R2': 'WT'}).
        If provided, a 'Group' column is added to the output.
        If None, 'Group' is set equal to 'Sample'.

    Returns
    -------
    pd.DataFrame
        One row per (Sample × Timepoint). Columns:
        • Time_h          – timepoint
        • Sample          – biological sample identifier
        • Group           – experimental group
        • Mean_OD         – mean Corrected_OD across technical replicates
        • SD_OD           – SD across technical replicates
        • N_TechReps      – number of valid technical replicates averaged
        • OD_Has_NaN      – True if any tech rep was NaN at this timepoint
    """
    # Exclude blank wells (Corrected_OD is NaN for blank rows)
    df = corrected_df[corrected_df["Blank_Used"] != "self"].copy()
    df = df.dropna(subset=["Corrected_OD"])

    agg = (
        df.groupby([sample_col, "Time_h"], sort=False)["Corrected_OD"]
        .agg(
            Mean_OD="mean",
            SD_OD="std",
            N_TechReps="count",
        )
        .reset_index()
        .rename(columns={sample_col: "Sample"})
    )

    # Flag timepoints where any tech rep was NaN
    n_raw = (
        df.groupby([sample_col, "Time_h"])["Corrected_OD"]
        .size()
        .reset_index(name="_raw_n")
        .rename(columns={sample_col: "Sample"})
    )
    agg = agg.merge(n_raw, on=["Sample", "Time_h"], how="left")
    agg["OD_Has_NaN"] = agg["N_TechReps"] < agg["_raw_n"]
    agg = agg.drop(columns=["_raw_n"])

    # Attach group labels
    if group_assignment:
        agg["Group"] = agg["Sample"].map(group_assignment).fillna(agg["Sample"])
    else:
        agg["Group"] = agg["Sample"]

    agg = agg.sort_values(["Group", "Sample", "Time_h"]).reset_index(drop=True)
    return agg


# =============================================================================
# STEP 4  —  Calculate growth metrics per biological sample
# =============================================================================

def calculate_growth_metrics(
    avg_df: pd.DataFrame,
    smoothing: bool = True,
    smoothing_window: int = 5,
    od_thresholds: list | None = None,
) -> pd.DataFrame:
    """
    Calculate 9 growth metrics for each biological sample.

    Operates on the output of average_technical_replicates() — one averaged
    OD trace per biological sample.

    Metrics calculated
    ------------------
    Max_OD            – maximum OD value in the trace
    Final_OD          – OD at the last timepoint
    AUC               – area under the OD-time curve (trapezoidal rule)
    Max_Growth_Rate   – maximum instantaneous growth rate (OD/h), from
                        smoothed first derivative
    Lag_Phase_h       – lag phase duration (h): x-intercept of the tangent
                        line at the inflection point, relative to initial OD
    Doubling_Time_h   – ln(2) / max specific growth rate (h)
                        max specific growth rate = Max_Growth_Rate / OD_at_inflection
    Time_to_Max_Rate  – timepoint at which Max_Growth_Rate occurs
    Time_to_OD_X      – time to reach each OD threshold (default [0.2, 0.5])

    Parameters
    ----------
    avg_df : pd.DataFrame
        Output of average_technical_replicates(). Must contain
        Time_h, Sample, Group, Mean_OD.
    smoothing : bool
        Apply Savitzky-Golay smoothing before derivative calculation.
        Reduces noise-driven spikes in growth rate. Default True.
    smoothing_window : int
        Window length for Savitzky-Golay filter (must be odd, ≥ 5).
        Default 5.
    od_thresholds : list or None
        OD values for Time_to_OD_X metrics. Default [0.2, 0.5].

    Returns
    -------
    pd.DataFrame
        One row per biological sample. Columns:
        Sample, Group, Max_OD, Final_OD, AUC,
        Max_Growth_Rate, Lag_Phase_h, Doubling_Time_h,
        Time_to_Max_Rate, Time_to_OD_0.2, Time_to_OD_0.5
        (or other thresholds if specified)
    """
    if od_thresholds is None:
        od_thresholds = [0.2, 0.5]

    results = []

    for (sample, group), grp in avg_df.groupby(["Sample", "Group"], sort=False):
        grp = grp.sort_values("Time_h")
        time = grp["Time_h"].values.astype(float)
        od   = grp["Mean_OD"].values.astype(float)

        # Remove leading NaNs
        valid = ~np.isnan(od)
        time  = time[valid]
        od    = od[valid]

        if len(time) < 4:
            # Not enough points for any meaningful metric
            row = _empty_metrics_row(sample, group, od_thresholds)
            results.append(row)
            continue

        # ── Smooth OD for derivative calculation ─────────────────────────────
        win = min(smoothing_window, len(od))
        if win % 2 == 0:
            win -= 1       # must be odd
        if smoothing and win >= 5:
            try:
                od_smooth = savgol_filter(od, window_length=win, polyorder=3)
                od_smooth = np.clip(od_smooth, 0, None)
            except Exception:
                od_smooth = od.copy()
        else:
            od_smooth = od.copy()

        # ── Basic metrics ─────────────────────────────────────────────────────
        max_od   = float(np.nanmax(od))
        final_od = float(od[-1])
        auc      = float(np.trapezoid(od, time) if hasattr(np, 'trapezoid') else np.trapz(od, time))

        # ── Growth rate (first derivative of smoothed OD) ────────────────────
        dt        = np.diff(time)
        dod       = np.diff(od_smooth)
        rates     = dod / dt                       # OD/h, absolute
        rate_time = (time[:-1] + time[1:]) / 2    # midpoints

        max_rate_idx = int(np.nanargmax(rates))
        max_rate     = float(rates[max_rate_idx])
        time_max_rate = float(rate_time[max_rate_idx])

        # ── Specific growth rate at inflection → doubling time ────────────────
        # Inflection point OD from smoothed curve at time_max_rate
        # Approximate by interpolating smoothed OD at time_max_rate
        od_at_inflection = float(np.interp(time_max_rate, time, od_smooth))
        if od_at_inflection > 0 and max_rate > 0:
            mu_max        = max_rate / od_at_inflection   # h⁻¹
            doubling_time = float(np.log(2) / mu_max)
        else:
            doubling_time = np.nan

        # ── Lag phase ─────────────────────────────────────────────────────────
        # Tangent line at inflection point: y = max_rate*(t - time_max_rate) + od_at_inflection
        # Lag = time where tangent crosses OD = od[0]
        if max_rate > 0:
            lag = time_max_rate - (od_at_inflection - od[0]) / max_rate
            lag = max(0.0, float(lag))
        else:
            lag = np.nan

        # ── Time to OD thresholds ─────────────────────────────────────────────
        threshold_times = {}
        for thresh in od_thresholds:
            tt = _time_to_threshold(time, od, thresh)
            threshold_times[f"Time_to_OD_{thresh}"] = tt

        row = {
            "Sample":          sample,
            "Group":           group,
            "Max_OD":          round(max_od, 4),
            "Final_OD":        round(final_od, 4),
            "AUC":             round(auc, 4),
            "Max_Growth_Rate": round(max_rate, 4),
            "Lag_Phase_h":     round(lag, 4) if not np.isnan(lag) else np.nan,
            "Doubling_Time_h": round(doubling_time, 4) if not np.isnan(doubling_time) else np.nan,
            "Time_to_Max_Rate": round(time_max_rate, 4),
        }
        row.update({k: (round(v, 4) if not np.isnan(v) else np.nan)
                    for k, v in threshold_times.items()})
        results.append(row)

    return pd.DataFrame(results).reset_index(drop=True)


def _time_to_threshold(
    time: np.ndarray,
    od: np.ndarray,
    threshold: float,
) -> float:
    """
    Return the time at which OD first reaches `threshold` by linear
    interpolation between flanking timepoints.
    Returns NaN if the threshold is never reached.
    """
    for i in range(len(od) - 1):
        if od[i] <= threshold <= od[i + 1]:
            # Linear interpolation
            frac = (threshold - od[i]) / (od[i + 1] - od[i])
            return float(time[i] + frac * (time[i + 1] - time[i]))
    return np.nan


def _empty_metrics_row(sample: str, group: str, od_thresholds: list) -> dict:
    """Return a row of NaN metrics for samples with insufficient data."""
    row = {
        "Sample": sample, "Group": group,
        "Max_OD": np.nan, "Final_OD": np.nan, "AUC": np.nan,
        "Max_Growth_Rate": np.nan, "Lag_Phase_h": np.nan,
        "Doubling_Time_h": np.nan, "Time_to_Max_Rate": np.nan,
    }
    for t in od_thresholds:
        row[f"Time_to_OD_{t}"] = np.nan
    return row


# =============================================================================
# STEP 5  —  Group summary statistics
# =============================================================================

def summarise_metrics(
    metrics_df: pd.DataFrame,
    group_col: str = "Group",
) -> pd.DataFrame:
    """
    Compute mean ± SD ± SEM per group across biological replicates.

    Each row in metrics_df is ONE biological sample.
    This function must never be called on technical-replicate data.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        Output of calculate_growth_metrics(). One row per biological sample.
    group_col : str
        Column identifying experimental groups.

    Returns
    -------
    pd.DataFrame
        One row per (group × metric). Columns:
        Group, Metric, N, Mean, SD, SEM
    """
    numeric_metrics = [
        c for c in metrics_df.columns
        if c not in ("Sample", "Group") and pd.api.types.is_numeric_dtype(metrics_df[c])
    ]

    rows = []
    for group, grp in metrics_df.groupby(group_col, sort=False):
        for metric in numeric_metrics:
            vals = grp[metric].dropna().values
            n    = len(vals)
            mean = float(np.mean(vals)) if n > 0 else np.nan
            sd   = float(np.std(vals, ddof=1)) if n > 1 else np.nan
            sem  = sd / np.sqrt(n) if n > 1 else np.nan
            rows.append({
                "Group":  group,
                "Metric": metric,
                "N":      n,
                "Mean":   round(mean, 4) if not np.isnan(mean) else np.nan,
                "SD":     round(sd,   4) if not np.isnan(sd)   else np.nan,
                "SEM":    round(sem,  4) if not np.isnan(sem)  else np.nan,
            })

    return pd.DataFrame(rows).reset_index(drop=True)


# =============================================================================
# STEP 6  —  Statistical testing
# =============================================================================

def run_statistical_tests(
    metrics_df: pd.DataFrame,
    group_col: str = "Group",
    test: str = "auto",
    correction_method: str = "fdr_bh",
) -> pd.DataFrame:
    """
    Compare growth metrics between groups using biological replicate values.

    Tests operate on individual metric values (one per biological sample),
    NOT on technical replicates.

    Parameters
    ----------
    metrics_df : pd.DataFrame
        Output of calculate_growth_metrics(). One row per biological sample.
    group_col : str
        Column identifying experimental groups.
    test : str
        'auto'     – Welch's t-test for 2 groups, ANOVA for ≥3 (default)
        't-test'   – Student's t-test (equal variance)
        'welch'    – Welch's t-test (unequal variance; recommended)
        'wilcoxon' – Mann-Whitney U (non-parametric, 2 groups)
        'anova'    – one-way ANOVA
        'kruskal'  – Kruskal-Wallis (non-parametric)
    correction_method : str
        'fdr_bh'    – Benjamini-Hochberg FDR (default)
        'bonferroni' – Bonferroni
        'none'      – raw p-values

    Returns
    -------
    pd.DataFrame
        One row per metric. Columns:
        Metric, Test, Groups_Compared, N_per_group,
        Statistic, p_value, p_adj, Significance
    """
    numeric_metrics = [
        c for c in metrics_df.columns
        if c not in ("Sample", "Group") and pd.api.types.is_numeric_dtype(metrics_df[c])
    ]
    groups   = metrics_df[group_col].unique().tolist()
    n_groups = len(groups)

    if n_groups < 2:
        raise ValueError(f"Need ≥2 groups for statistical testing. Found: {groups}")

    valid_tests = {"auto", "t-test", "welch", "wilcoxon", "anova", "kruskal"}
    if test not in valid_tests:
        raise ValueError(f"Unknown test '{test}'. Choose from: {sorted(valid_tests)}")

    resolved_test = test
    if test == "auto":
        resolved_test = "welch" if n_groups == 2 else "anova"

    if resolved_test in ("t-test", "welch", "wilcoxon") and n_groups > 2:
        raise ValueError(
            f"Test '{resolved_test}' only supports 2 groups. "
            f"Your data has {n_groups} groups. Use 'anova' or 'kruskal'."
        )

    results = []
    for metric in numeric_metrics:
        group_arrays = [
            metrics_df.loc[metrics_df[group_col] == g, metric].dropna().values
            for g in groups
        ]
        n_per_group = " / ".join(str(len(a)) for a in group_arrays)
        valid = [(g, a) for g, a in zip(groups, group_arrays) if len(a) > 0]

        if len(valid) < 2:
            results.append({
                "Metric": metric, "Test": resolved_test,
                "Groups_Compared": " vs ".join(groups),
                "N_per_group": n_per_group,
                "Statistic": np.nan, "p_value": np.nan,
            })
            continue

        valid_groups, valid_arrays = zip(*valid)
        groups_label = " vs ".join(valid_groups)

        try:
            if resolved_test == "t-test":
                stat, pval = scipy_stats.ttest_ind(valid_arrays[0], valid_arrays[1], equal_var=True)
            elif resolved_test == "welch":
                stat, pval = scipy_stats.ttest_ind(valid_arrays[0], valid_arrays[1], equal_var=False)
            elif resolved_test == "wilcoxon":
                stat, pval = scipy_stats.mannwhitneyu(valid_arrays[0], valid_arrays[1], alternative="two-sided")
            elif resolved_test == "anova":
                stat, pval = scipy_stats.f_oneway(*valid_arrays)
            elif resolved_test == "kruskal":
                stat, pval = scipy_stats.kruskal(*valid_arrays)
            else:
                stat, pval = np.nan, np.nan
        except Exception:
            stat, pval = np.nan, np.nan

        results.append({
            "Metric":          metric,
            "Test":            resolved_test,
            "Groups_Compared": groups_label,
            "N_per_group":     n_per_group,
            "Statistic":       round(float(stat), 4) if not np.isnan(stat) else np.nan,
            "p_value":         float(pval) if not np.isnan(pval) else np.nan,
        })

    stats_df = pd.DataFrame(results)

    # Multiple-testing correction
    raw_pvals = stats_df["p_value"].values
    if correction_method == "none" or stats_df["p_value"].isna().all():
        stats_df["p_adj"] = raw_pvals
    else:
        mask = ~np.isnan(raw_pvals)
        p_adj = raw_pvals.copy()
        if mask.sum() > 0:
            _, p_corrected, _, _ = multipletests(raw_pvals[mask], method=correction_method, alpha=0.05)
            p_adj[mask] = p_corrected
        stats_df["p_adj"] = p_adj

    stats_df["Significance"] = stats_df["p_adj"].apply(_pval_to_stars)

    for col in ["Statistic", "p_value", "p_adj"]:
        stats_df[col] = stats_df[col].round(4)

    col_order = ["Metric", "Test", "Groups_Compared", "N_per_group",
                 "Statistic", "p_value", "p_adj", "Significance"]
    stats_df = stats_df[[c for c in col_order if c in stats_df.columns]]
    return stats_df.reset_index(drop=True)


# =============================================================================
# Internal helpers
# =============================================================================

def _pval_to_stars(p: float) -> str:
    """Convert a p-value to a significance star string."""
    if pd.isna(p):      return "n/a"
    if p < 0.0001:      return "****"
    if p < 0.001:       return "***"
    if p < 0.01:        return "**"
    if p < 0.05:        return "*"
    return "ns"
