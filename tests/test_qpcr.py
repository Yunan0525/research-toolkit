# =============================================================================
# tests/test_qpcr.py  –  Unit Tests for analysis/qpcr.py
# =============================================================================
# HOW TO RUN
#   pytest tests/test_qpcr.py -v
#
# PIPELINE UNDER TEST
#   average_technical_replicates()  →  calculate_delta_ct()
#   → calculate_delta_delta_ct()    →  summarise_results()
#   → run_statistical_tests()
#
# TEST DATA DESIGN
#   TECH_REP_DATA  – has technical triplicates; used to test averaging
#   SAMPLE_DATA    – pre-averaged (1 row per bio sample × gene); used for math
#
# HAND-CALCULATED REFERENCE VALUES (SAMPLE_DATA)
#   S1 (Control) : ΔCt(IL6) = 28.0 − 20.0 = 8.0
#   S2 (Control) : ΔCt(IL6) = 27.0 − 20.0 = 7.0
#   mean control ΔCt(IL6)   = 7.5
#   S3 (Treatment): ΔΔCt(IL6) = 4.0 − 7.5 = −3.5  →  FC = 2^3.5 ≈ 11.31
#   S4 (Treatment): ΔΔCt(IL6) = 3.0 − 7.5 = −4.5
# =============================================================================

import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

import pytest
import pandas as pd
import numpy as np

from analysis.qpcr import (
    average_technical_replicates,
    calculate_delta_ct,
    calculate_delta_delta_ct,
    summarise_results,
    run_statistical_tests,
)

# =============================================================================
# Shared fixtures
# =============================================================================

# Pre-averaged data — one row per biological sample × gene
SAMPLE_DATA = pd.DataFrame({
    "Sample": ["S1","S1","S2","S2","S3","S3","S4","S4"],
    "Gene":   ["IL6","GAPDH","IL6","GAPDH","IL6","GAPDH","IL6","GAPDH"],
    "Ct":     [28.0, 20.0,   27.0, 20.0,   24.0, 20.0,   23.0, 20.0],
    "Group":  ["Control","Control","Control","Control",
               "Treatment","Treatment","Treatment","Treatment"],
})

# Raw data with technical triplicates (3 rows per bio sample × gene)
TECH_REP_DATA = pd.DataFrame({
    "Sample": (["S1"]*6 + ["S2"]*6),
    "Gene":   (["IL6","IL6","IL6","GAPDH","GAPDH","GAPDH"] * 2),
    "Ct":     [28.0,28.2,27.8, 20.0,20.1,19.9,   # S1 IL6 mean=28.0, GAPDH mean=20.0
               27.0,27.3,26.7, 20.0,20.2,19.8],   # S2 IL6 mean=27.0, GAPDH mean=20.0
    "Group":  (["Control"]*6 + ["Control"]*6),
})


# =============================================================================
# Convenience pipeline runner
# =============================================================================

def _full_pipeline(df=None):
    """Run the mandatory 4-step pipeline on SAMPLE_DATA (already averaged)."""
    if df is None:
        df = SAMPLE_DATA
    dct    = calculate_delta_ct(df, "Sample", "Gene", "Ct", "GAPDH")
    detail = calculate_delta_delta_ct(dct, "Group", "Control")
    summ   = summarise_results(detail, "Group")
    return dct, detail, summ


# =============================================================================
# 1 — average_technical_replicates
# =============================================================================

class TestAverageTechReps:

    def test_row_count_after_averaging(self):
        """6 tech-rep rows per sample → 2 rows after averaging (IL6 + GAPDH)."""
        avg = average_technical_replicates(
            TECH_REP_DATA, "Sample", "Gene", "Ct", "Group"
        )
        assert len(avg) == 4, f"Expected 4 rows (2 samples × 2 genes), got {len(avg)}"

    def test_mean_ct_correct(self):
        """S1 IL6: raw = [28.0, 28.2, 27.8] → mean = 28.0."""
        avg = average_technical_replicates(
            TECH_REP_DATA, "Sample", "Gene", "Ct", "Group"
        )
        val = avg.loc[(avg["Sample"]=="S1") & (avg["Gene"]=="IL6"), "Ct"].values[0]
        assert np.isclose(val, 28.0, atol=1e-6), f"Expected 28.0, got {val}"

    def test_tech_sd_present(self):
        """Output must contain Ct_Tech_SD and Ct_Tech_N columns."""
        avg = average_technical_replicates(
            TECH_REP_DATA, "Sample", "Gene", "Ct", "Group"
        )
        assert "Ct_Tech_SD" in avg.columns, "Ct_Tech_SD missing"
        assert "Ct_Tech_N"  in avg.columns, "Ct_Tech_N missing"

    def test_tech_n_equals_three(self):
        """Each group had 3 tech reps; N should equal 3."""
        avg = average_technical_replicates(
            TECH_REP_DATA, "Sample", "Gene", "Ct", "Group"
        )
        assert (avg["Ct_Tech_N"] == 3).all(), \
            f"Expected Ct_Tech_N=3, got: {avg['Ct_Tech_N'].tolist()}"

    def test_group_col_preserved(self):
        """Group column must survive averaging."""
        avg = average_technical_replicates(
            TECH_REP_DATA, "Sample", "Gene", "Ct", "Group"
        )
        assert "Group" in avg.columns
        assert (avg["Group"] == "Control").all()

    def test_single_row_is_no_op(self):
        """If there is already 1 row per sample × gene, output equals input."""
        avg = average_technical_replicates(
            SAMPLE_DATA, "Sample", "Gene", "Ct", "Group"
        )
        assert len(avg) == len(SAMPLE_DATA)


# =============================================================================
# 2 — calculate_delta_ct
# =============================================================================

class TestCalculateDeltaCt:

    def test_arithmetic(self):
        """S1 IL6 ΔCt = 28.0 − 20.0 = 8.0."""
        dct, _, _ = _full_pipeline()
        val = dct.loc[dct["Sample"]=="S1", "Delta_Ct"].values[0]
        assert np.isclose(val, 8.0), f"Expected 8.0, got {val}"

    def test_reference_gene_excluded(self):
        """GAPDH rows must not appear in the ΔCt output."""
        dct, _, _ = _full_pipeline()
        assert "GAPDH" not in dct["Gene"].values

    def test_ref_ct_column_present(self):
        """Ref_Ct column must exist and equal 20.0 for all SAMPLE_DATA rows."""
        dct, _, _ = _full_pipeline()
        assert "Ref_Ct" in dct.columns
        assert np.allclose(dct["Ref_Ct"], 20.0)

    def test_missing_reference_gene_raises(self):
        with pytest.raises(ValueError, match="not found"):
            calculate_delta_ct(SAMPLE_DATA, "Sample", "Gene", "Ct", "ACTB")

    def test_missing_sample_reference_raises(self):
        """Sample S1 missing its GAPDH row → ValueError."""
        bad = SAMPLE_DATA[
            ~((SAMPLE_DATA["Sample"]=="S1") & (SAMPLE_DATA["Gene"]=="GAPDH"))
        ].copy()
        with pytest.raises(ValueError, match="S1"):
            calculate_delta_ct(bad, "Sample", "Gene", "Ct", "GAPDH")


# =============================================================================
# 3 — calculate_delta_delta_ct
# =============================================================================

class TestCalculateDeltaDeltaCt:

    def test_ddct_arithmetic(self):
        """S3 IL6: ΔCt=4.0, mean ctrl ΔCt=7.5 → ΔΔCt=−3.5."""
        _, detail, _ = _full_pipeline()
        val = detail.loc[detail["Sample"]=="S3", "Delta_Delta_Ct"].values[0]
        assert np.isclose(val, -3.5, atol=1e-6), f"Expected −3.5, got {val}"

    def test_fold_change_arithmetic(self):
        """FC = 2^(−ΔΔCt) = 2^3.5 ≈ 11.314."""
        _, detail, _ = _full_pipeline()
        val = detail.loc[detail["Sample"]=="S3", "Fold_Change"].values[0]
        assert np.isclose(val, 2**3.5, atol=1e-4), f"Expected {2**3.5:.4f}, got {val}"

    def test_log2fc_equals_negative_ddct(self):
        """log2FC = −ΔΔCt for every row."""
        _, detail, _ = _full_pipeline()
        assert np.allclose(detail["log2FC"], -detail["Delta_Delta_Ct"])

    def test_mean_log2fc_control_zero(self):
        """Mean log2FC of the control group = 0 (by definition of ΔΔCt)."""
        _, detail, _ = _full_pipeline()
        ctrl = detail.loc[detail["Group"]=="Control", "log2FC"].mean()
        assert np.isclose(ctrl, 0.0, atol=1e-9), f"Control mean log2FC = {ctrl}"

    def test_treatment_upregulated(self):
        """Treatment IL6 has lower Ct → FC > 1 for all treatment samples."""
        _, detail, _ = _full_pipeline()
        tx = detail.loc[detail["Group"]=="Treatment", "Fold_Change"]
        assert (tx > 1).all(), f"Expected all FC > 1, got: {tx.values}"

    def test_n_control_biorePs_correct(self):
        """N_Control_BioReps should equal 2 (S1, S2)."""
        _, detail, _ = _full_pipeline()
        assert (detail["N_Control_BioReps"] == 2).all()

    def test_missing_control_group_raises(self):
        dct = calculate_delta_ct(SAMPLE_DATA, "Sample", "Gene", "Ct", "GAPDH")
        with pytest.raises(ValueError, match="not found"):
            calculate_delta_delta_ct(dct, "Group", "Untreated")


# =============================================================================
# 4 — summarise_results  (SD, not SEM, is the primary dispersion metric)
# =============================================================================

class TestSummariseResults:

    def test_required_columns_present(self):
        _, _, summ = _full_pipeline()
        required = {
            "Gene", "Group", "N_BioReps",
            "Mean_FC", "SD_FC", "SEM_FC",
            "Mean_log2FC", "SD_log2FC", "SEM_log2FC",
            "Mean_Delta_Ct", "SD_Delta_Ct",
        }
        missing = required - set(summ.columns)
        assert not missing, f"Summary missing: {missing}"

    def test_n_biorePs_correct(self):
        """SAMPLE_DATA has 2 biological replicates per group."""
        _, _, summ = _full_pipeline()
        assert (summ["N_BioReps"] == 2).all(), \
            f"Expected N_BioReps=2, got: {summ['N_BioReps'].values}"

    def test_sd_not_sem_is_primary(self):
        """SD_FC must be present and ≥ SEM_FC (SD ≥ SEM always)."""
        _, _, summ = _full_pipeline()
        valid = summ.dropna(subset=["SD_FC", "SEM_FC"])
        assert (valid["SD_FC"] >= valid["SEM_FC"]).all()

    def test_sd_nan_for_single_replicate(self):
        """With N=1 biological replicate, SD must be NaN."""
        single = SAMPLE_DATA[SAMPLE_DATA["Sample"].isin(["S1","S3"])].copy()
        dct    = calculate_delta_ct(single, "Sample","Gene","Ct","GAPDH")
        detail = calculate_delta_delta_ct(dct, "Group", "Control")
        summ   = summarise_results(detail, "Group")
        assert summ["SD_FC"].isna().all(), "SD_FC must be NaN when N_BioReps=1"

    def test_control_mean_fc_near_one(self):
        """Mean FC of control group should be close to 1.0."""
        _, _, summ = _full_pipeline()
        ctrl_fc = summ.loc[summ["Group"]=="Control", "Mean_FC"].values[0]
        # With identical Ct references, control mean FC = exactly 1.0
        assert np.isclose(ctrl_fc, 1.0, atol=0.01), \
            f"Control mean FC = {ctrl_fc}"


# =============================================================================
# 5 — run_statistical_tests
# =============================================================================

class TestStatisticalTests:

    def _dct(self):
        return calculate_delta_ct(SAMPLE_DATA, "Sample","Gene","Ct","GAPDH")

    def test_returns_dataframe(self):
        result = run_statistical_tests(self._dct(), "Group")
        assert isinstance(result, pd.DataFrame)

    def test_one_row_per_gene(self):
        """One result row per target gene."""
        result = run_statistical_tests(self._dct(), "Group")
        n_genes = SAMPLE_DATA[SAMPLE_DATA["Gene"] != "GAPDH"]["Gene"].nunique()
        assert len(result) == n_genes, \
            f"Expected {n_genes} rows, got {len(result)}"

    def test_required_columns(self):
        result = run_statistical_tests(self._dct(), "Group")
        required = {"Gene", "Test", "p_value", "p_adj", "Significance", "Statistic"}
        missing = required - set(result.columns)
        assert not missing, f"Stats table missing: {missing}"

    def test_pvalue_between_0_and_1(self):
        result = run_statistical_tests(self._dct(), "Group")
        valid = result["p_value"].dropna()
        assert ((valid >= 0) & (valid <= 1)).all(), \
            f"p-values out of range: {valid.values}"

    def test_significance_stars_assigned(self):
        result = run_statistical_tests(self._dct(), "Group")
        valid_stars = {"****","***","**","*","ns","n/a"}
        for s in result["Significance"]:
            assert s in valid_stars, f"Unexpected significance value: {s!r}"

    def test_welch_default_for_two_groups(self):
        result = run_statistical_tests(self._dct(), "Group", test="auto")
        assert (result["Test"] == "welch").all()

    def test_fdr_correction_applied(self):
        """BH-adjusted p-values must be ≥ raw p-values."""
        result = run_statistical_tests(
            self._dct(), "Group", correction_method="fdr_bh"
        )
        valid = result.dropna(subset=["p_value","p_adj"])
        assert (valid["p_adj"] >= valid["p_value"] - 1e-9).all()

    def test_bonferroni_correction(self):
        result = run_statistical_tests(
            self._dct(), "Group", correction_method="bonferroni"
        )
        assert "p_adj" in result.columns
        assert not result["p_adj"].isna().all()

    def test_invalid_test_raises(self):
        with pytest.raises(ValueError, match="Unknown test"):
            run_statistical_tests(self._dct(), "Group", test="magic_test")

    def test_missing_delta_ct_raises(self):
        with pytest.raises(ValueError, match="Delta_Ct"):
            run_statistical_tests(SAMPLE_DATA, "Group")  # no Delta_Ct column
