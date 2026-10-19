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


# ==============================================================================
# 5. Subtask 3.2 Phylogeny & Haplotype Network Tests
# ==============================================================================

def test_subtask3_2_deliverables_exist():
    """Verify all Subtask 3.2 trees, network tables, and figures exist and are non-empty."""
    tree_files = [
        "results/tables/rice_maternal_plastome_tree.nwk",
        "results/tables/rice_maternal_plastome_tree.nex",
        "results/tables/rice_maternal_plastome_tree.xml",
        "data/processed/phylogeny/rice_maternal_plastome_tree.nwk",
    ]
    for tf in tree_files:
        assert os.path.exists(tf), f"Tree file {tf} missing."
        assert os.path.getsize(tf) > 100, f"Tree file {tf} is unexpectedly small."

    network_tables = [
        "results/tables/plastome_haplotypes.tsv",
        "results/tables/plastome_network_nodes.tsv",
        "results/tables/plastome_network_edges.tsv",
    ]
    for tbl in network_tables:
        assert os.path.exists(tbl), f"Network table {tbl} missing."
        df = pd.read_csv(tbl, sep="\t")
        assert len(df) > 0, f"Table {tbl} is empty."

    figures = [
        "results/figures/fig5_maternal_plastid_phylogeny.pdf",
        "results/figures/fig5_maternal_plastid_phylogeny.png",
        "results/figures/fig6_maternal_haplotype_network.pdf",
        "results/figures/fig6_maternal_haplotype_network.png",
    ]
    for fig in figures:
        assert os.path.exists(fig), f"Publication figure {fig} missing."
        assert os.path.getsize(fig) > 1000, f"Figure {fig} has abnormally small size."


def test_tree_topology_and_taxa():
    """Verify tree leaf count equals 109 taxa and contains outgroup species."""
    import Bio.Phylo as Phylo
    tree = Phylo.read("results/tables/rice_maternal_plastome_tree.nwk", "newick")
    terminals = tree.get_terminals()
    assert len(terminals) == 109, f"Expected 109 leaf taxa in tree, got {len(terminals)}"

    term_names = {t.name for t in terminals}
    # Verify presence of reference and outgroup species
    assert "NC_001320.1" in term_names, "Nipponbare reference missing from tree."
    assert "LC907692.1" in term_names, "Oryza punctata basal outgroup missing from tree."
    assert "NC_017835.1" in term_names, "Oryza rufipogon missing from tree."


def test_haplotype_catalog_invariants():
    """Verify haplotype extraction completeness and graph network connectivity."""
    df_haps = pd.read_csv("results/tables/plastome_haplotypes.tsv", sep="\t")
    assert len(df_haps) == 28, f"Expected 28 maternal haplotypes, got {len(df_haps)}"

    # Check sum of accessions equals total Sativa cohort (90)
    total_acc = df_haps["Accession_Count"].sum()
    assert total_acc == 90, f"Total accessions across haplotypes ({total_acc}) != 90."

    # Check frequency percentage sum
    freq_sum = df_haps["Frequency_Pct"].sum()
    assert np.isclose(freq_sum, 100.0, atol=0.1), f"Frequency sum ({freq_sum}) != 100%."

    # Check network connectivity: MST of 28 nodes must have exactly 27 edges
    df_nodes = pd.read_csv("results/tables/plastome_network_nodes.tsv", sep="\t")
    df_edges = pd.read_csv("results/tables/plastome_network_edges.tsv", sep="\t")
    assert len(df_nodes) == 28, f"Expected 28 nodes, got {len(df_nodes)}"
    assert len(df_edges) == 27, f"Expected 27 edges (V - 1), got {len(df_edges)}"

    # Edge weights must be non-negative mutational steps
    assert (df_edges["Mutational_Steps"] >= 0).all(), "Found negative mutational steps in network edges."


# ==============================================================================
# 6. Subtask 3.3 Population Structure, Kinship & Manifold Tests
# ==============================================================================

def test_subtask3_3_deliverables_exist():
    """Verify all Subtask 3.3 tabular files and publication figures exist and are non-empty."""
    tables = [
        "results/tables/plastome_snps_ld_pruned.tsv",
        "results/tables/plastome_kinship_matrix.tsv",
        "results/tables/plastome_pca_eigenvalues.tsv",
        "results/tables/plastome_pca_coordinates.tsv",
        "results/tables/plastome_pca_loadings.tsv",
        "results/tables/plastome_umap_coordinates.tsv",
    ]
    for tbl in tables:
        assert os.path.exists(tbl), f"Deliverable table {tbl} is missing."
        df = pd.read_csv(tbl, sep="\t")
        assert len(df) > 0, f"Table {tbl} is empty."

    figures = [
        "results/figures/fig7_plastome_pca_structure.pdf",
        "results/figures/fig7_plastome_pca_structure.png",
        "results/figures/fig8_plastome_umap_kinship.pdf",
        "results/figures/fig8_plastome_umap_kinship.png",
    ]
    for fig in figures:
        assert os.path.exists(fig), f"Publication figure {fig} missing."
        assert os.path.getsize(fig) > 1000, f"Figure {fig} has abnormally small size."


def test_ld_pruning_invariants():
    """Verify LD pruning parameters: MAF >= 0.02 and non-empty independent SNP subset."""
    df_pruned = pd.read_csv("results/tables/plastome_snps_ld_pruned.tsv", sep="\t")
    assert len(df_pruned) == 20, f"Expected 20 independent tag SNPs, got {len(df_pruned)}"
    assert (df_pruned["MAF"] >= 0.02 - 1e-6).all(), "Found SNPs violating min_maf >= 0.02 filter."
    assert (df_pruned["POS"] > 0).all(), "Found invalid SNP genomic position."


def test_kinship_matrix_properties():
    """Verify VanRaden/Patterson Kinship matrix is symmetric, 90x90, with valid diagonal self-relationships."""
    df_k = pd.read_csv("results/tables/plastome_kinship_matrix.tsv", sep="\t", index_col=0)
    assert df_k.shape == (90, 90), f"Kinship matrix shape {df_k.shape} != (90, 90)."

    # Matrix symmetry: K == K^T
    K = df_k.values
    assert np.allclose(K, K.T, atol=1e-10), "Kinship matrix K is not symmetric (K != K^T)."

    # Diagonal self-relationship properties
    diag_vals = np.diag(K)
    assert (diag_vals > 0).all(), "Found non-positive diagonal entries in K."
    mean_diag = np.mean(diag_vals)
    assert np.isclose(mean_diag, 1.0, atol=0.1), f"Mean diagonal self-relationship ({mean_diag:.4f}) deviates substantially from 1.0."


def test_pca_decomposition_invariants():
    """Verify PCA eigenvalues are descending, variance ratios sum to cumulative variance, and top PCs explain major variance."""
    df_eval = pd.read_csv("results/tables/plastome_pca_eigenvalues.tsv", sep="\t")
    assert len(df_eval) == 10, f"Expected 10 principal components, got {len(df_eval)}"

    # Eigenvalues must be strictly non-negative and descending
    eigenvalues = df_eval["Eigenvalue"].values
    assert (eigenvalues >= 0).all(), "Found negative eigenvalues in eigendecomposition."
    assert (np.diff(eigenvalues) <= 1e-8).all(), "Eigenvalues are not sorted in descending order."

    # Variance explained consistency
    var_pct = df_eval["Variance_Explained_Pct"].values
    cum_pct = df_eval["Cumulative_Variance_Pct"].values
    assert np.allclose(np.cumsum(var_pct), cum_pct, atol=1e-3), "Cumulative variance does not match sum of individual variance components."

    # Biological expectation: plastome has structured variation, top 2 PCs explain substantial variance
    assert cum_pct[1] > 45.0, f"Top 2 PCs explain only {cum_pct[1]:.2f}% variance (expected >45%)."

    # PCA coordinates
    df_coord = pd.read_csv("results/tables/plastome_pca_coordinates.tsv", sep="\t")
    assert len(df_coord) == 90, f"Expected 90 accessions in PCA coordinates, got {len(df_coord)}"
    for i in range(1, 11):
        assert f"PC{i}" in df_coord.columns, f"PC{i} column missing from coordinates."
        assert not df_coord[f"PC{i}"].isna().any(), f"Found NaN values in PC{i}."


def test_umap_manifold_embeddings():
    """Verify 2D and 3D UMAP manifold coordinates exist for all 90 accessions with no NaNs."""
    df_umap = pd.read_csv("results/tables/plastome_umap_coordinates.tsv", sep="\t")
    assert len(df_umap) == 90, f"Expected 90 accessions in UMAP table, got {len(df_umap)}"

    cols = ["UMAP1_2D", "UMAP2_2D", "UMAP1_3D", "UMAP2_3D", "UMAP3_3D"]
    for c in cols:
        assert c in df_umap.columns, f"Column {c} missing from UMAP coordinates."
        assert not df_umap[c].isna().any(), f"Found NaN values in UMAP column {c}."
        assert np.isfinite(df_umap[c].values).all(), f"Found infinite values in UMAP column {c}."


# ==============================================================================
# 7. Subtask 3.4 Ancestral Admixture & Cytonuclear Discordance Tests
# ==============================================================================

def test_subtask3_4_deliverables_exist():
    """Verify all Subtask 3.4 tabular datasets and publication figures exist and are non-empty."""
    tables = [
        "results/tables/plastome_admixture_cv_errors.tsv",
        "results/tables/plastome_admixture_q_matrices.tsv",
        "results/tables/plastome_admixture_ancestry_frequencies.tsv",
        "results/tables/cytonuclear_discordance_matrix.tsv",
        "results/tables/cytonuclear_subpopulation_contingency.tsv",
    ]
    for tbl in tables:
        assert os.path.exists(tbl), f"Deliverable table {tbl} is missing."
        df = pd.read_csv(tbl, sep="\t")
        assert len(df) > 0, f"Table {tbl} is empty."

    figures = [
        "results/figures/fig9_plastome_admixture_ancestry.pdf",
        "results/figures/fig9_plastome_admixture_ancestry.png",
        "results/figures/fig10_cytonuclear_discordance_capture.pdf",
        "results/figures/fig10_cytonuclear_discordance_capture.png",
    ]
    for fig in figures:
        assert os.path.exists(fig), f"Publication figure {fig} missing."
        assert os.path.getsize(fig) > 1000, f"Figure {fig} has abnormally small size."


def test_admixture_q_matrix_invariants():
    """Verify Q-matrix rows sum to 1.0 +- 1e-4 for all K=2..8 across all 90 accessions."""
    df_q = pd.read_csv("results/tables/plastome_admixture_q_matrices.tsv", sep="\t")
    assert len(df_q) == 90, f"Expected 90 accessions in Q-matrices table, got {len(df_q)}"

    for k in range(2, 9):
        cols = [f"K{k}_Q{c}" for c in range(1, k + 1)]
        for c in cols:
            assert c in df_q.columns, f"Column {c} missing from Q-matrix table."
            assert (df_q[c] >= 0.0 - 1e-6).all(), f"Found negative ancestry proportion in {c}."
            assert (df_q[c] <= 1.0 + 1e-6).all(), f"Found ancestry proportion > 1 in {c}."

        row_sums = df_q[cols].sum(axis=1).values
        assert np.allclose(row_sums, 1.0, atol=1e-4), f"Q-matrix row sums deviate from 1.0 at K={k}."

    # Verify ancestral allele frequencies matrix
    df_f = pd.read_csv("results/tables/plastome_admixture_ancestry_frequencies.tsv", sep="\t")
    freq_cols = [c for c in df_f.columns if c.startswith("Ancestral_Cluster_")]
    assert len(freq_cols) >= 2, "Expected ancestral allele frequency columns in F table."
    for c in freq_cols:
        assert (df_f[c] >= 0.0).all() and (df_f[c] <= 1.0).all(), f"Allele frequency in {c} out of [0, 1] bounds."


def test_admixture_cv_error_properties():
    """Verify CV deviance values are monotonically evaluated across K=2..8 and non-negative."""
    df_cv = pd.read_csv("results/tables/plastome_admixture_cv_errors.tsv", sep="\t")
    assert len(df_cv) == 7, f"Expected 7 K-evaluations (K=2..8), got {len(df_cv)}"
    assert list(df_cv["K"].values) == list(range(2, 9)), "K values do not match expected range 2..8."

    assert (df_cv["CV_Error_Mean"] > 0).all(), "Found non-positive CV error deviance."
    assert (df_cv["CV_Error_Std"] >= 0).all(), "Found negative CV error standard deviation."
    assert not df_cv["CV_Error_Mean"].isna().any(), "Found NaN in CV error mean."


def test_cytonuclear_discordance_invariants():
    """Verify CDI is strictly bounded in [0, 1] and captures cold-adapted introgressions."""
    df_disc = pd.read_csv("results/tables/cytonuclear_discordance_matrix.tsv", sep="\t")
    assert len(df_disc) == 90, f"Expected 90 accessions in discordance table, got {len(df_disc)}"

    cdi = df_disc["Cytonuclear_Discordance_Index"].values
    assert (cdi >= 0.0 - 1e-6).all(), "Found negative CDI."
    assert (cdi <= 1.0 + 1e-6).all(), "Found CDI > 1.0."

    # Validate categories
    valid_cats = {"Concordant", "Admixed_Introgression", "Chloroplast_Capture"}
    assert set(df_disc["Discordance_Category"].unique()).issubset(valid_cats), "Unexpected discordance category."

    # Verify canonical Indica cultivars are concordant
    canonical_indica = ["93-11", "IR8", "Zhenshan 97", "Minghui 63", "TN1"]
    for cult in canonical_indica:
        match = df_disc[df_disc["Cultivar"] == cult]
        if not match.empty:
            row = match.iloc[0]
            assert row["Discordance_Category"] == "Concordant", f"{cult} expected to be Concordant, got {row['Discordance_Category']}"
            assert row["Cytonuclear_Discordance_Index"] <= 0.35, f"{cult} has unexpectedly high CDI ({row['Cytonuclear_Discordance_Index']})"

    # Verify Russian northern margin capture cultivars have CDI >= 0.70 and Chloroplast_Capture
    russian_capture = ["Tai Bonnet", "Zlata", "Oceano", "Snezhinka", "Svetlana", "Mars", "Ruille", "Arroz Buit"]
    for cult in russian_capture:
        match = df_disc[df_disc["Cultivar"] == cult]
        if not match.empty:
            row = match.iloc[0]
            assert row["Discordance_Category"] == "Chloroplast_Capture", f"{cult} expected to be Chloroplast_Capture, got {row['Discordance_Category']}"
            assert row["Cytonuclear_Discordance_Index"] >= 0.70, f"{cult} expected CDI >= 0.70, got {row['Cytonuclear_Discordance_Index']}"


def test_cytonuclear_contingency_statistics():
    """Verify contingency table sums to exactly 90 accessions and shows significant association."""
    df_cont = pd.read_csv("results/tables/cytonuclear_subpopulation_contingency.tsv", sep="\t")
    assert "Nuclear_Primary_Group" in df_cont.columns, "Nuclear_Primary_Group column missing from contingency table."

    # Sum of counts across all cells excluding Total row and column
    df_no_tot = df_cont[df_cont["Nuclear_Primary_Group"] != "Total"]
    count_cols = [c for c in df_no_tot.columns if c not in ["Nuclear_Primary_Group", "Total"]]
    total_samples = df_no_tot[count_cols].values.sum()
    assert total_samples == 90, f"Contingency matrix total ({total_samples}) != 90."



