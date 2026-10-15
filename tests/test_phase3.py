"""
Unit Tests for Phase 3 - Subtask 3.1: Population Diversity, Neutrality & Divergence
===================================================================================
Verifies mathematical boundary invariants, Hudson's Fst symmetry, exact callable
site accounting, and the integrity of generated tabular and visual deliverables.
"""

import os
import pytest
import numpy as np
import pandas as pd

from scripts.compute_diversity_and_divergence import (
    compute_sfs_statistics,
    compute_hudson_components,
    weighted_block_jackknife,
)


# ==============================================================================
# 1. Deliverables Existence & Completeness Tests
# ==============================================================================

def test_phase3_tables_exist():
    """Verify all 5 generated tabular datasets exist and are non-empty."""
    tables = [
        "results/tables/plastome_diversity_neutrality_summary.tsv",
        "results/tables/plastome_compartment_diversity.tsv",
        "results/tables/plastome_pairwise_divergence_fst.tsv",
        "results/tables/plastome_sliding_window_divergence.tsv",
        "results/tables/plastome_gene_selection_diversity.tsv",
    ]
    for tbl in tables:
        assert os.path.exists(tbl), f"Deliverable table {tbl} is missing."
        df = pd.read_csv(tbl, sep="\t")
        assert len(df) > 0, f"Table {tbl} is empty."


def test_phase3_figures_exist():
    """Verify all 4 publication vector figures exist in PDF and PNG format."""
    figures = [
        "results/figures/fig1_plastome_diversity_neutrality.pdf",
        "results/figures/fig1_plastome_diversity_neutrality.png",
        "results/figures/fig1_supp_plastome_doughnuts.pdf",
        "results/figures/fig1_supp_plastome_doughnuts.png",
        "results/figures/fig2_sliding_window_divergence_scan.pdf",
        "results/figures/fig2_sliding_window_divergence_scan.png",
        "results/figures/fig3_pairwise_differentiation_heatmaps.pdf",
        "results/figures/fig3_pairwise_differentiation_heatmaps.png",
        "results/figures/fig4_gene_selection_fst_tajima.pdf",
        "results/figures/fig4_gene_selection_fst_tajima.png",
    ]
    for fig in figures:
        assert os.path.exists(fig), f"Publication figure {fig} missing."
        assert os.path.getsize(fig) > 1000, f"Figure {fig} has abnormally small size."


# ==============================================================================
# 2. Mathematical Boundary Invariants
# ==============================================================================

def test_diversity_invariants():
    """Verify mathematical boundary invariants: Pi >= 0, Theta_W >= 0, Tajima's D consistency."""
    df_div = pd.read_csv("results/tables/plastome_diversity_neutrality_summary.tsv", sep="\t")

    assert (df_div["Pi_Per_Site"] >= 0.0).all(), "Found negative nucleotide diversity (Pi)."
    assert (df_div["Theta_W_Per_Site"] >= 0.0).all(), "Found negative Watterson's Theta."
    assert (df_div["Segregating_Sites_S"] >= 0).all(), "Found negative segregating site counts."

    # Tajima's D sign test: D > 0 iff Pi > Theta_W (when |Pi - Theta_W| is noticeable)
    for _, r in df_div.iterrows():
        d_val = r["Tajima_D"]
        pi_val = r["Pi_Per_Site"]
        th_val = r["Theta_W_Per_Site"]
        if np.isfinite(d_val) and abs(pi_val - th_val) > 1e-8:
            if pi_val > th_val:
                assert d_val > 0, f"Expected positive Tajima's D when Pi > Theta_W in {r['Population']}"
            else:
                assert d_val < 0, f"Expected negative Tajima's D when Pi < Theta_W in {r['Population']}"


def test_hudson_fst_permutation_symmetry():
    """Verify Hudson's Fst estimator satisfies exact symmetry under population permutation."""
    # Toy mathematical test
    np.random.seed(42)
    n1, n2 = 30, 25
    p1 = np.random.uniform(0.0, 1.0, 100)
    p2 = np.random.uniform(0.0, 1.0, 100)

    num12, den12 = compute_hudson_components(p1, p2, n1, n2)
    num21, den21 = compute_hudson_components(p2, p1, n2, n1)

    fst12 = np.sum(num12) / np.sum(den12)
    fst21 = np.sum(num21) / np.sum(den21)

    assert np.isclose(fst12, fst21, atol=1e-12), "Hudson's Fst failed population permutation symmetry."


def test_hudson_fst_bounds_in_dataset():
    """Verify that Hudson's Fst estimates in real output tables satisfy 0 <= Fst <= 1."""
    df_diff = pd.read_csv("results/tables/plastome_pairwise_divergence_fst.tsv", sep="\t")
    assert (df_diff["Hudson_Fst"] >= 0.0).all(), "Found negative Hudson's Fst."
    assert (df_diff["Hudson_Fst"] <= 1.0).all(), "Found Hudson's Fst exceeding 1.0."
    assert (df_diff["d_XY_Per_kb"] >= 0.0).all(), "Found negative d_XY."


# ==============================================================================
# 3. Denominator Integrity & Callable Accounting
# ==============================================================================

def test_callable_denominator_integrity():
    """Verify exact callable denominators match callable mask values from Phase 2."""
    df_div = pd.read_csv("results/tables/plastome_diversity_neutrality_summary.tsv", sep="\t")

    wp_rows = df_div[df_div["Genomic_Partition"] == "Whole_Plastome"]
    sc_rows = df_div[df_div["Genomic_Partition"] == "Single_Copy_Partition"]

    assert (wp_rows["Callable_Sites_bp"] == 133526).all(), "Whole Plastome callable denominator mismatch."
    assert (sc_rows["Callable_Sites_bp"] == 112814).all(), "Single Copy Partition callable denominator mismatch."


def test_inverted_repeats_symmetry():
    """Verify that IRa and IRb exhibit identical or highly conserved segregating site counts."""
    df_comp = pd.read_csv("results/tables/plastome_compartment_diversity.tsv", sep="\t")

    irb = df_comp[df_comp["Compartment"] == "IRb"].set_index("Population")
    ira = df_comp[df_comp["Compartment"] == "IRa"].set_index("Population")

    for pop in ["Japonica", "Indica_Canonical", "Northern_Margin_Russia", "All_Sativa", "All_Sativa_Canonical"]:
        if pop in irb.index and pop in ira.index:
            s_irb = irb.loc[pop, "Segregating_Sites_S"]
            s_ira = ira.loc[pop, "Segregating_Sites_S"]
            assert s_irb == s_ira, f"IR copy divergence detected for {pop}: IRb={s_irb}, IRa={s_ira}"


# ==============================================================================
# 4. Sliding Window & Candidate Gene Scan Tests
# ==============================================================================

def test_sliding_window_continuity():
    """Verify sliding window profile has no missing gaps and spans the full chromosome."""
    df_win = pd.read_csv("results/tables/plastome_sliding_window_divergence.tsv", sep="\t")
    assert len(df_win) == 130, f"Expected 130 sliding windows, got {len(df_win)}"

    # Check that windows advance monotonically by step size 1000 bp
    diffs = np.diff(df_win["Start_bp"].values)
    assert (diffs == 1000).all(), "Sliding windows do not advance with constant step size 1000 bp."


def test_hypothesis_h2_candidate_genes():
    """Verify that all core photosynthetic complexes and H2 candidate genes are evaluated."""
    df_gene = pd.read_csv("results/tables/plastome_gene_selection_diversity.tsv", sep="\t")
    assert len(df_gene) == 106, f"Expected 106 genes, got {len(df_gene)}"

    h2_genes = ["psbA", "matK", "rbcL", "ndhF", "ndhA", "ndhB", "ndhD"]
    for gene in h2_genes:
        match = df_gene[df_gene["Gene_Symbol"] == gene]
        assert not match.empty, f"Hypothesis H2 target gene {gene} missing from gene scan."
        row = match.iloc[0]
        assert row["Is_H2_Target"] or gene.startswith("ndh"), f"{gene} not flagged as H2 target."
        assert row["Callable_bp"] > 0, f"{gene} has invalid callable length."
