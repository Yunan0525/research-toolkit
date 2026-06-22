# =============================================================================
# analysis/qpcr.py  –  qPCR ΔΔCt Analysis Logic
# =============================================================================
# PURPOSE
#   Pure Python / pandas / scipy functions for the comparative ΔΔCt method
#   with technical replicate handling and statistical testing.
#   NO Streamlit imports — works identically in notebooks, scripts, or APIs.
#
# MANDATORY PIPELINE (always call in this order)
# ─────────────────────────────────────────────
#   1. average_technical_replicates()
#        Input : raw data (may have multiple rows per Sample × Gene)
#        Output: one row per biological Sample × Gene (mean Ct, SD, N tech reps)
#
#   2. calculate_delta_ct()
#        Input : averaged data from step 1
#        Output: ΔCt per biological sample × target gene
#                ΔCt = mean_Ct(target, sample) − mean_Ct(reference, sample)
#
#   3. calculate_delta_delta_ct()
#        Input : ΔCt table from step 2
#        Output: ΔΔCt, Fold Change, log2FC per biological sample × gene
#                ΔΔCt = ΔCt(sample) − mean( ΔCt(all control biological replicates) )
#
#   4. summarise_results()
#        Input : per-biological-replicate results from step 3
#        Output: mean ± SD per Gene × Group across biological replicates
#                (SD chosen over SEM because N of biological replicates is small)
#
#   5. run_statistical_tests()   [optional but recommended]
#        Input : per-biological-replicate ΔCt values from step 2
#        Output: p-values and test statistics comparing groups per gene
#                Tests operate on ΔCt (the right unit for stats, not FC)
#
# DESIGN DECISIONS
# ─────────────────
#   • Averaging tech reps is MANDATORY, not optional.
#     Fold Changes and statistics must be computed on biological replicates.
#   • Statistics are performed on ΔCt values, not on Fold Change.
#     Fold Change is log-normally distributed; ΔCt is approximately normal,
#     making parametric tests (t-test, ANOVA) appropriate.
#   • SD is reported in summaries instead of SEM because biological replicate
#     numbers are typically small (n=3–6) and SD better characterises
#     biological variability for publication.
#   • SEM is still computed and included for completeness.
#
# MATHEMATICS
# ────────────
#   mean_Ct(gene, sample)  = arithmetic mean of technical replicates
#   ΔCt(target, sample)    = mean_Ct(target, sample) − mean_Ct(reference, sample)
#   mean_ΔCt(gene, control)= mean of ΔCt across all control biological replicates
#   ΔΔCt(sample)           = ΔCt(sample) − mean_ΔCt(gene, control)
#   Fold Change            = 2^(−ΔΔCt)
#   log2FC                 = −ΔΔCt
#
# ASSUMPTION: equal primer efficiency for all genes (standard ΔΔCt assumption).
# =============================================================================

import pandas as pd
import numpy as np
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests


# =============================================================================
# STEP 1  —  Average technical replicates
# =============================================================================

def average_technical_replicates(
    df: pd.DataFrame,
    sample_col: str,
    gene_col: str,
    ct_col: str,
    group_col: str,
) -> pd.DataFrame:
    """
    Collapse technical replicates into one mean Ct per biological sample × gene.

    This MUST be called before calculate_delta_ct().  It is not optional:
    statistics and fold changes must be computed on biological replicates,
    not on technical replicates.

    Parameters
    ----------
    df : pd.DataFrame
        Raw data. May contain multiple rows per Sample × Gene (tech reps).
        Required columns: sample_col, gene_col, ct_col, group_col.
    sample_col : str
        Column identifying biological samples (e.g. "Sample").
    gene_col : str
        Column identifying genes (e.g. "Gene").
    ct_col : str
        Column with raw numeric Ct values (e.g. "Ct").
    group_col : str
        Column identifying experimental groups (e.g. "Group").
        Preserved in output; used downstream by calculate_delta_delta_ct.

    Returns
    -------
    pd.DataFrame
        One row per (sample × gene). Columns:
        • sample_col          – biological sample identifier
        • gene_col            – gene name
        • group_col           – experimental group
        • ct_col              – mean Ct across technical replicates
        • "Ct_Tech_SD"        – SD of Ct across technical replicates
        • "Ct_Tech_N"         – number of technical replicates averaged
        High Ct_Tech_SD (> 0.5) may indicate a problematic well; check raw data.

    Notes
    -----
    Non-numeric Ct values (e.g. "Undetermined") are coerced to NaN and
    excluded from the average.  A warning column "Ct_Has_NaN" flags samples
    where at least one technical replicate was excluded.
    """
    df = df.copy()

    # Coerce Ct to numeric; non-numeric becomes NaN
    df[ct_col] = pd.to_numeric(df[ct_col], errors="coerce")

    # Group by biological identity: sample + gene (+ group to carry it through)
    group_keys = [sample_col, gene_col, group_col]

    agg = (
        df.groupby(group_keys, sort=False)[ct_col]
        .agg(
            **{
                ct_col:         "mean",   # mean Ct used for all downstream math
                "Ct_Tech_SD":   "std",    # SD within technical replicates
                "Ct_Tech_N":    "count",  # number of valid (non-NaN) tech reps
            }
        )
        .reset_index()
    )

    # Flag samples where any tech rep was NaN (excluded from mean)
    n_raw = df.groupby(group_keys)[ct_col].size().reset_index(name="_raw_n")
    agg = agg.merge(n_raw, on=group_keys, how="left")
    agg["Ct_Has_NaN"] = agg["Ct_Tech_N"] < agg["_raw_n"]
    agg = agg.drop(columns=["_raw_n"])

    return agg.reset_index(drop=True)


# =============================================================================
# STEP 2  —  Calculate ΔCt  (biological-sample level)
# =============================================================================

def calculate_delta_ct(
    df: pd.DataFrame,
    sample_col: str,
    gene_col: str,
    ct_col: str,
    reference_gene: str,
) -> pd.DataFrame:
    """
    Compute ΔCt = mean_Ct(target, sample) − mean_Ct(reference, sample).

    Must be called on the output of average_technical_replicates() so that
    one averaged Ct per biological sample × gene is used.

    Parameters
    ----------
    df : pd.DataFrame
        Averaged data (output of average_technical_replicates()).
        One row per biological sample × gene.
    sample_col : str
        Column with biological sample identifiers.
    gene_col : str
        Column with gene names.
    ct_col : str
        Column with mean Ct values (after tech-rep averaging).
    reference_gene : str
        Gene name to use as the housekeeping / normaliser gene.

    Returns
    -------
    pd.DataFrame
        One row per (biological sample × target gene). Reference gene rows
        are removed (ΔCt of the reference with itself = 0 by definition and
        carries no biological information). Columns added:
        • "Ref_Ct"        – mean Ct of reference gene for this sample
        • "Ref_Ct_SD"     – tech-rep SD of reference gene (if present)
        • "Delta_Ct"      – Ct(target) − Ct(reference)

    Raises
    ------
    ValueError
        If reference_gene absent, or any sample lacks a reference measurement.
    """
    df = df.copy()

    # ── Validate ──────────────────────────────────────────────────────────────
    if reference_gene not in df[gene_col].values:
        raise ValueError(
            f"Reference gene '{reference_gene}' not found in column '{gene_col}'. "
            f"Available genes: {sorted(df[gene_col].unique().tolist())}"
        )

    # ── Split reference from target rows ──────────────────────────────────────
    ref_cols = [sample_col, ct_col]
    if "Ct_Tech_SD" in df.columns:
        ref_cols.append("Ct_Tech_SD")

    ref_df = df[df[gene_col] == reference_gene][ref_cols].copy()
    rename_map = {ct_col: "Ref_Ct"}
    if "Ct_Tech_SD" in ref_df.columns:
        rename_map["Ct_Tech_SD"] = "Ref_Ct_SD"
    ref_df = ref_df.rename(columns=rename_map)

    target_df = df[df[gene_col] != reference_gene].copy()

    # ── Check every target sample has a reference measurement ─────────────────
    missing = set(target_df[sample_col].unique()) - set(ref_df[sample_col].unique())
    if missing:
        raise ValueError(
            f"Samples missing a '{reference_gene}' measurement: "
            f"{sorted(missing)}. Every biological sample must have a "
            "reference gene row."
        )

    # ── Merge and compute ΔCt ─────────────────────────────────────────────────
    merged = target_df.merge(ref_df, on=sample_col, how="left")
    merged["Delta_Ct"] = merged[ct_col] - merged["Ref_Ct"]

    return merged.reset_index(drop=True)


# =============================================================================
# STEP 3  —  Calculate ΔΔCt, Fold Change, log2FC
# =============================================================================

def calculate_delta_delta_ct(
    df: pd.DataFrame,
    group_col: str,
    control_group: str,
) -> pd.DataFrame:
    """
    Compute ΔΔCt, Fold Change, and log2FC for each biological replicate.

    ΔΔCt = ΔCt(biological replicate) − mean(ΔCt of control biological replicates)

    The control mean is computed per gene, using only biological replicates
    (technical replicates have already been collapsed in step 1).

    Parameters
    ----------
    df : pd.DataFrame
        Output of calculate_delta_ct(). One row per biological sample × gene.
    group_col : str
        Column with experimental group labels.
    control_group : str
        Value in group_col to use as the reference group.

    Returns
    -------
    pd.DataFrame
        Same rows as input, with added columns:
        • "Mean_Delta_Ct_Control"   – per-gene mean ΔCt of control bio-reps
        • "N_Control_BioReps"       – how many control bio-reps contributed
        • "Delta_Delta_Ct"          – ΔΔCt
        • "Fold_Change"             – 2^(−ΔΔCt)
        • "log2FC"                  – −ΔΔCt

    Raises
    ------
    ValueError
        If 'Delta_Ct' column is missing or control_group not found.
    """
    if "Delta_Ct" not in df.columns:
        raise ValueError("'Delta_Ct' column not found. Run calculate_delta_ct() first.")
    if control_group not in df[group_col].values:
        raise ValueError(
            f"Control group '{control_group}' not found in '{group_col}'. "
            f"Available: {sorted(df[group_col].unique().tolist())}"
        )

    df = df.copy()
    gene_col = _infer_gene_col(df, group_col)

    # Compute mean ΔCt of control group per gene, over biological replicates
    ctrl_df = df[df[group_col] == control_group]
    ctrl_stats = (
        ctrl_df.groupby(gene_col)["Delta_Ct"]
        .agg(
            Mean_Delta_Ct_Control="mean",
            N_Control_BioReps="count",
        )
        .reset_index()
    )

    df = df.merge(ctrl_stats, on=gene_col, how="left")

    df["Delta_Delta_Ct"] = df["Delta_Ct"] - df["Mean_Delta_Ct_Control"]
    df["Fold_Change"]    = 2.0 ** (-df["Delta_Delta_Ct"])
    df["log2FC"]         = -df["Delta_Delta_Ct"]

    return df.reset_index(drop=True)


# =============================================================================
# STEP 4  —  Summarise biological replicates  (mean ± SD per Gene × Group)
# =============================================================================

def summarise_results(
    df: pd.DataFrame,
    group_col: str,
) -> pd.DataFrame:
    """
    Aggregate per-biological-replicate values into a summary table.

    Reports mean ± SD (and SEM) per Gene × Group.

    SD is the primary dispersion metric because it quantifies biological
    variability between replicates, which is the scientifically meaningful
    quantity.  SEM is included for completeness.

    Parameters
    ----------
    df : pd.DataFrame
        Output of calculate_delta_delta_ct().
    group_col : str
        Column with experimental group labels.

    Returns
    -------
    pd.DataFrame
        One row per (Gene × Group). Columns:
        Gene, Group, N_BioReps,
        Mean_FC, SD_FC, SEM_FC,
        Mean_log2FC, SD_log2FC, SEM_log2FC,
        Mean_Delta_Ct, SD_Delta_Ct, SEM_Delta_Ct,
        Mean_Delta_Delta_Ct, SD_Delta_Delta_Ct

    Raises
    ------
    ValueError
        If required analysis columns are missing.
    """
    required = ["Fold_Change", "log2FC", "Delta_Ct", "Delta_Delta_Ct"]
    missing_cols = [c for c in required if c not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Missing columns: {missing_cols}. "
            "Ensure calculate_delta_delta_ct() has been called first."
        )

    gene_col = _infer_gene_col(df, group_col)

    def sd(x):
        return x.std(ddof=1) if x.count() > 1 else np.nan

    def sem(x):
        n = x.count()
        return x.std(ddof=1) / np.sqrt(n) if n > 1 else np.nan

    summary = (
        df.groupby([gene_col, group_col], sort=False)
        .agg(
            N_BioReps              = ("Fold_Change",       "count"),
            Mean_FC                = ("Fold_Change",       "mean"),
            SD_FC                  = ("Fold_Change",       sd),
            SEM_FC                 = ("Fold_Change",       sem),
            Mean_log2FC            = ("log2FC",            "mean"),
            SD_log2FC              = ("log2FC",            sd),
            SEM_log2FC             = ("log2FC",            sem),
            Mean_Delta_Ct          = ("Delta_Ct",          "mean"),
            SD_Delta_Ct            = ("Delta_Ct",          sd),
            SEM_Delta_Ct           = ("Delta_Ct",          sem),
            Mean_Delta_Delta_Ct    = ("Delta_Delta_Ct",    "mean"),
            SD_Delta_Delta_Ct      = ("Delta_Delta_Ct",    sd),
        )
        .reset_index()
        .rename(columns={gene_col: "Gene", group_col: "Group"})
    )

    float_cols = [c for c in summary.columns if summary[c].dtype == float]
    summary[float_cols] = summary[float_cols].round(4)

    return summary


# =============================================================================
# STEP 5  —  Statistical testing on biological replicates
# =============================================================================

def run_statistical_tests(
    df: pd.DataFrame,
    group_col: str,
    test: str = "auto",
    correction_method: str = "fdr_bh",
) -> pd.DataFrame:
    """
    Compare ΔCt values between experimental groups for each target gene.

    Tests are performed on ΔCt (the normalised Ct value), NOT on Fold Change.
    Rationale: ΔCt is approximately normally distributed; Fold Change is
    log-normally distributed and should not be used directly in t-tests.

    Group structure:
        2 groups  → t-test (parametric) or Wilcoxon (non-parametric)
        ≥3 groups → one-way ANOVA (parametric) or Kruskal-Wallis (non-parametric)

    Parameters
    ----------
    df : pd.DataFrame
        Output of calculate_delta_ct() (before or after calculate_delta_delta_ct).
        Must contain "Delta_Ct" and group_col.
    group_col : str
        Column with experimental group labels.
    test : str
        Statistical test to run. Options:
        "auto"       – t-test for 2 groups, ANOVA for ≥3 (default)
        "t-test"     – Student's independent t-test (2 groups only)
        "welch"      – Welch's t-test (unequal variance; recommended over Student's)
        "wilcoxon"   – Wilcoxon rank-sum (Mann-Whitney U; non-parametric, 2 groups)
        "anova"      – one-way ANOVA (≥2 groups)
        "kruskal"    – Kruskal-Wallis (non-parametric equivalent of ANOVA)
    correction_method : str
        Multiple-testing correction applied across all gene p-values.
        Options: "fdr_bh" (Benjamini-Hochberg, default), "bonferroni", "none".

    Returns
    -------
    pd.DataFrame
        One row per gene. Columns:
        • "Gene"            – gene name
        • "Test"            – name of the statistical test used
        • "Groups_Compared" – e.g. "Control vs Treatment"
        • "Statistic"       – test statistic (t, F, or H)
        • "p_value"         – raw p-value
        • "p_adj"           – multiple-testing corrected p-value
        • "Significance"    – "****" / "***" / "**" / "*" / "ns"
        • "N_per_group"     – n for each group as a string (e.g. "3 / 3")
        One additional column per group pair for post-hoc (ANOVA only):
        future extension — not implemented yet.

    Raises
    ------
    ValueError
        If fewer than 2 groups, or test name is invalid.
    """
    if "Delta_Ct" not in df.columns:
        raise ValueError(
            "'Delta_Ct' column not found. Run calculate_delta_ct() first."
        )

    gene_col = _infer_gene_col(df, group_col)
    groups   = df[group_col].unique().tolist()
    n_groups = len(groups)

    if n_groups < 2:
        raise ValueError(
            f"Statistical testing requires at least 2 groups. "
            f"Found only: {groups}"
        )

    # ── Resolve which test to run ─────────────────────────────────────────────
    valid_tests = {"auto", "t-test", "welch", "wilcoxon", "anova", "kruskal"}
    if test not in valid_tests:
        raise ValueError(
            f"Unknown test '{test}'. Choose from: {sorted(valid_tests)}"
        )

    if test == "auto":
        test = "welch" if n_groups == 2 else "anova"

    if test in ("t-test", "welch", "wilcoxon") and n_groups > 2:
        raise ValueError(
            f"Test '{test}' only supports 2 groups. "
            f"Your data has {n_groups} groups. Use 'anova' or 'kruskal'."
        )

    # ── Run test per gene ─────────────────────────────────────────────────────
    results = []

    for gene, gene_df in df.groupby(gene_col):
        group_arrays = [
            gene_df.loc[gene_df[group_col] == g, "Delta_Ct"].dropna().values
            for g in groups
        ]
        n_per_group = " / ".join(str(len(a)) for a in group_arrays)

        # Filter out groups with no data (safety)
        valid = [(g, a) for g, a in zip(groups, group_arrays) if len(a) > 0]
        if len(valid) < 2:
            results.append({
                "Gene": gene, "Test": test,
                "Groups_Compared": " vs ".join(groups),
                "Statistic": np.nan, "p_value": np.nan,
                "N_per_group": n_per_group,
            })
            continue

        valid_groups, valid_arrays = zip(*valid)
        groups_label = " vs ".join(valid_groups)

        try:
            if test == "t-test":
                stat, pval = scipy_stats.ttest_ind(
                    valid_arrays[0], valid_arrays[1], equal_var=True
                )
            elif test == "welch":
                stat, pval = scipy_stats.ttest_ind(
                    valid_arrays[0], valid_arrays[1], equal_var=False
                )
            elif test == "wilcoxon":
                stat, pval = scipy_stats.mannwhitneyu(
                    valid_arrays[0], valid_arrays[1], alternative="two-sided"
                )
            elif test == "anova":
                stat, pval = scipy_stats.f_oneway(*valid_arrays)
            elif test == "kruskal":
                stat, pval = scipy_stats.kruskal(*valid_arrays)
            else:
                stat, pval = np.nan, np.nan

        except Exception:
            stat, pval = np.nan, np.nan

        results.append({
            "Gene":             gene,
            "Test":             test,
            "Groups_Compared":  groups_label,
            "Statistic":        round(float(stat), 4) if not np.isnan(stat) else np.nan,
            "p_value":          float(pval) if not np.isnan(pval) else np.nan,
            "N_per_group":      n_per_group,
        })

    stats_df = pd.DataFrame(results)

    # ── Multiple-testing correction ───────────────────────────────────────────
    raw_pvals = stats_df["p_value"].values

    if correction_method == "none" or stats_df["p_value"].isna().all():
        stats_df["p_adj"] = raw_pvals
    else:
        # multipletests requires no NaN; mask them
        mask = ~np.isnan(raw_pvals)
        p_adj = raw_pvals.copy()
        if mask.sum() > 0:
            _, p_corrected, _, _ = multipletests(
                raw_pvals[mask],
                method=correction_method,
                alpha=0.05,
            )
            p_adj[mask] = p_corrected
        stats_df["p_adj"] = p_adj

    # ── Significance stars ────────────────────────────────────────────────────
    stats_df["Significance"] = stats_df["p_adj"].apply(_pval_to_stars)

    # ── Round floats ──────────────────────────────────────────────────────────
    for col in ["Statistic", "p_value", "p_adj"]:
        stats_df[col] = stats_df[col].round(4)

    # ── Column order ──────────────────────────────────────────────────────────
    col_order = [
        "Gene", "Test", "Groups_Compared",
        "N_per_group", "Statistic", "p_value", "p_adj", "Significance",
    ]
    stats_df = stats_df[[c for c in col_order if c in stats_df.columns]]

    return stats_df.reset_index(drop=True)


# =============================================================================
# Internal helpers
# =============================================================================

def _pval_to_stars(p: float) -> str:
    """Convert a p-value to a significance star string."""
    if pd.isna(p):
        return "n/a"
    if p < 0.0001:
        return "****"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def _infer_gene_col(df: pd.DataFrame, group_col: str) -> str:
    """
    Identify which column holds gene names.

    Priority 1: a column literally named Gene / gene / GENE.
    Priority 2: first object (string) column that is not the group column.
    """
    for candidate in df.columns:
        if candidate.lower() == "gene":
            return candidate

    candidates = [
        c for c in df.columns
        if df[c].dtype == object and c != group_col
    ]
    if candidates:
        return candidates[0]

    raise ValueError(
        "Could not detect the gene column automatically. "
        "Ensure your data has a column named 'Gene'."
    )
