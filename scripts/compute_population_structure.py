"""
Population Structure, Kinship & Manifold Projection Pipeline
============================================================
Subtask 3.3: Global Rice Plastome Population Genomics

Executes linkage disequilibrium (LD) pruning, genomic relationship matrix
(GRM / Kinship K) computation, Principal Component Analysis (PCA), and
2D/3D UMAP manifold projections across global rice germplasm (n = 90).

Key Methodological & Statistical Invariants:
1. LD-Pruning: Sliding-window r^2 < 0.2 pruning to remove redundant co-segregating
   alleles along the non-recombining plastome molecule.
2. Kinship Matrix: VanRaden (2008) / Patterson et al. (2006) normalized covariance:
      K = (1 / M) * W @ W.T
   where W_ik = (G_ik - p_k) / sqrt(p_k * (1 - p_k)).
3. Eigendecomposition: Exact SVD on standardized genotype matrix W to extract
   eigenvalues, variance explained percentages, and principal component coordinates.
4. Manifold Projections: Uniform Manifold Approximation & Projection (UMAP)
   in 2D and 3D embedding spaces preserving local and global geometric relationships.
5. Publication Figures:
   - Figure 7: PCA Population Structure Composite (Scree plot, PC1 vs PC2, PC1 vs PC3).
   - Figure 8: 2D & 3D UMAP Projections + Hierarchical Kinship Heatmap.
"""

import os
import gzip
import time
import warnings
from typing import List, Dict, Tuple, Any, Optional

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import umap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import seaborn as sns
from adjustText import adjust_text

warnings.filterwarnings("ignore", category=UserWarning)


# ==============================================================================
# 1. Pipeline Engine Class
# ==============================================================================

class PlastomePopulationStructureEngine:
    """Engine for LD-pruning, Kinship K, PCA, and UMAP manifold inference."""

    def __init__(
        self,
        vcf_path: str = "data/processed/variants/rice_plastome_variants_polarized.vcf.gz",
        sativa_metadata_path: str = "results/tables/ncbi_oryza_sativa_chloroplast_genomes.tsv",
        crossref_path: str = "results/tables/plastome_nuclear_crossref.tsv",
        haplotypes_path: str = "results/tables/plastome_haplotypes.tsv",
        out_table_dir: str = "results/tables",
        out_fig_dir: str = "results/figures",
    ):
        self.vcf_path = vcf_path
        self.sativa_metadata_path = sativa_metadata_path
        self.crossref_path = crossref_path
        self.haplotypes_path = haplotypes_path
        self.out_table_dir = out_table_dir
        self.out_fig_dir = out_fig_dir

        self.samples: List[str] = []
        self.sample_metadata: pd.DataFrame = pd.DataFrame()
        self.raw_variants_df: pd.DataFrame = pd.DataFrame()
        self.raw_gt_matrix: np.ndarray = np.array([])  # (M, N)

        self.pruned_variants_df: pd.DataFrame = pd.DataFrame()
        self.pruned_gt_matrix: np.ndarray = np.array([])  # (M_pruned, N)
        self.kinship_matrix: np.ndarray = np.array([])  # (N, N)

        self.pca_model: Optional[PCA] = None
        self.pca_coords_df: pd.DataFrame = pd.DataFrame()
        self.eigenvalues_df: pd.DataFrame = pd.DataFrame()
        self.loadings_df: pd.DataFrame = pd.DataFrame()
        self.umap_coords_df: pd.DataFrame = pd.DataFrame()

    # ==========================================================================
    # 2. Data Ingestion & Metadata Harmonization
    # ==========================================================================

    def load_data(self) -> None:
        """Load sample metadata and polarized VCF genotype records."""
        print("--> Ingesting sample passport metadata and variant calls...")

        # 1. Load Sativa metadata
        df_sat = pd.read_csv(self.sativa_metadata_path, sep=r"\s{2,}|\t", engine="python")
        df_sat.columns = [c.strip() for c in df_sat.columns]
        for c in df_sat.columns:
            df_sat[c] = df_sat[c].astype(str).str.strip()
        sat_dict = df_sat.set_index("Accession").to_dict("index")

        # 2. Load crossref table
        df_xref = pd.read_csv(self.crossref_path, sep="\t")
        df_xref.columns = [c.strip() for c in df_xref.columns]
        for c in df_xref.columns:
            df_xref[c] = df_xref[c].astype(str).str.strip()
        xref_dict = df_xref.set_index("Plastome_Accession").to_dict("index")

        # 3. Read VCF sample headers
        with gzip.open(self.vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    self.samples = line.strip().split("\t")[9:]
                    break

        print(f"    Loaded {len(self.samples)} accessions from VCF header.")

        # Map metadata
        meta_list = []
        for s in self.samples:
            if s == "NC_001320.1":
                taxon = "Oryza sativa Japonica Group"
                cultivar = "Nipponbare"
                origin = "Japan"
                group = "Japonica"
                region = "East_Asia"
            else:
                rec_sat = sat_dict.get(s, {})
                rec_x = xref_dict.get(s, {})
                taxon = rec_sat.get("Taxon_Organism", "Unknown")
                cultivar = rec_sat.get("Cultivar", "Unknown")
                origin = rec_sat.get("Nation_Origin", "Unknown")

                if "Japonica" in taxon or "japonica" in taxon:
                    group = "Japonica"
                elif "Indica" in taxon or "indica" in taxon:
                    group = "Indica"
                elif "spontanea" in taxon:
                    group = "Spontanea"
                elif "aromatic" in taxon:
                    group = "Aromatic"
                else:
                    group = "Other_Sativa"

                if origin in ["Russia", "Europe"]:
                    region = "Northern_Margin_Russia"
                elif origin in ["China", "Japan", "South Korea", "Taiwan", "Taiwan / Japan"]:
                    region = "East_Asia"
                elif origin in ["India", "Pakistan / India", "Philippines (IRRI)", "Thailand"]:
                    region = "South_Southeast_Asia"
                elif origin in ["Saudi Arabia"]:
                    region = "Middle_East"
                else:
                    region = "Other_Unspecified"

            meta_list.append({
                "Sample_ID": s,
                "Cultivar": cultivar,
                "Taxon": taxon,
                "Varietal_Group": group,
                "Nation_Origin": origin,
                "Region": region,
                "Is_Structural_Inversion": s in ["CP054688.1", "CP056064.1"],
            })

        self.sample_metadata = pd.DataFrame(meta_list)

        # 4. Ingest VCF variants
        var_rows = []
        gt_rows = []
        with gzip.open(self.vcf_path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if parts[6] != "PASS":
                    continue

                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                info_tokens = dict(tok.split("=", 1) if "=" in tok else (tok, True) for tok in parts[7].split(";"))

                aa = info_tokens.get("AA", ref)
                domain = info_tokens.get("QUAD", "Unknown")
                eff = info_tokens.get("EFF", "intergenic_variant")
                gene = info_tokens.get("GENE", ".")
                cplx = info_tokens.get("COMPLEX", ".")

                gts = parts[9:]
                gt_nums = [0 if g == "0" else (1 if g == "1" else -1) for g in gts]

                var_rows.append({
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "AA": aa,
                    "DOMAIN": domain,
                    "EFF": eff,
                    "GENE": gene,
                    "COMPLEX": cplx,
                })
                gt_rows.append(gt_nums)

        self.raw_variants_df = pd.DataFrame(var_rows)
        self.raw_gt_matrix = np.array(gt_rows, dtype=np.int8)  # (M, N)
        print(f"    Loaded {len(self.raw_variants_df):,} PASS variants across {self.raw_gt_matrix.shape[1]} accessions.")

    # ==========================================================================
    # 3. Linkage Disequilibrium Pruning (r^2 < 0.2)
    # ==========================================================================

    def perform_ld_pruning(
        self,
        min_maf: float = 0.02,
        r2_threshold: float = 0.2,
        window_bp: int = 20000,
    ) -> pd.DataFrame:
        """Perform physical genomic sliding-window LD pruning with r^2 < 0.2 threshold."""
        print(f"--> Performing Linkage Disequilibrium (LD) pruning (MAF >= {min_maf}, window = {window_bp:,} bp, r^2 < {r2_threshold})...")
        t0 = time.time()

        # Compute allele frequency and MAF
        n_samples = self.raw_gt_matrix.shape[1]
        allele_counts = np.sum(self.raw_gt_matrix == 1, axis=1)
        freq = allele_counts / n_samples
        maf = np.minimum(freq, 1.0 - freq)

        self.raw_variants_df["DAF"] = freq
        self.raw_variants_df["MAF"] = maf

        # Filter by MAF >= min_maf
        maf_mask = maf >= min_maf
        maf_variants = self.raw_variants_df[maf_mask].reset_index(drop=True)
        maf_gt = self.raw_gt_matrix[maf_mask, :]
        positions = maf_variants["POS"].values
        n_maf = len(maf_variants)

        print(f"    Variants passing MAF >= {min_maf}: {n_maf:,} SNPs (from {len(self.raw_variants_df):,} raw).")

        # Physical sliding-window LD pruning with fast reversed lookback
        retained_indices = [0]
        for i in range(1, n_maf):
            pos_i = positions[i]
            g_i = maf_gt[i, :]

            is_independent = True
            for r in reversed(retained_indices):
                if pos_i - positions[r] > window_bp:
                    break
                g_r = maf_gt[r, :]
                corr = np.corrcoef(g_i, g_r)[0, 1]
                if not np.isnan(corr) and (corr ** 2) >= r2_threshold:
                    is_independent = False
                    break

            if is_independent:
                retained_indices.append(i)

        self.pruned_variants_df = maf_variants.iloc[retained_indices].reset_index(drop=True)
        self.pruned_gt_matrix = maf_gt[retained_indices, :]  # (M_pruned, N)

        print(f"    Retained {len(self.pruned_variants_df):,} independent tag SNPs (r^2 < {r2_threshold}) in {time.time() - t0:.2f}s.")
        return self.pruned_variants_df

    # ==========================================================================
    # 4. Genomic Relationship Matrix (Kinship K)
    # ==========================================================================

    def compute_kinship_matrix(self) -> np.ndarray:
        """Compute VanRaden (2008) / Patterson et al. (2006) normalized Kinship Matrix."""
        print("--> Computing Genomic Relationship Matrix (Kinship K)...")
        t0 = time.time()

        # Genotype matrix: accessions x SNPs
        X = self.pruned_gt_matrix.T.astype(np.float64)  # (N, M_pruned)
        n_samples, n_snps = X.shape

        # Allele frequencies
        p_k = np.mean(X == 1, axis=0)
        p_k = np.clip(p_k, 1e-4, 1.0 - 1e-4)

        # Standardized genotype matrix: W_ik = (G_ik - p_k) / sqrt(p_k * (1 - p_k))
        std_k = np.sqrt(p_k * (1.0 - p_k))
        W = (X - p_k) / std_k

        # Realized kinship matrix K = (1 / M) * W @ W.T
        K = (W @ W.T) / n_snps
        self.kinship_matrix = K

        print(f"    Computed {n_samples}x{n_samples} Kinship Matrix in {time.time() - t0:.3f}s.")
        print(f"    Mean diagonal self-relationship: {np.mean(np.diag(K)):.4f} (theoretical expectation = 1.0).")
        return K

    # ==========================================================================
    # 5. Principal Component Analysis (PCA)
    # ==========================================================================

    def perform_pca(self, n_components: int = 10) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Perform PCA on standardized genotype matrix and compute loadings."""
        print(f"--> Executing Principal Component Analysis (top {n_components} components)...")
        t0 = time.time()

        # Standardize matrix W
        X = self.pruned_gt_matrix.T.astype(np.float64)  # (N, M_pruned)
        p_k = np.mean(X == 1, axis=0)
        p_k = np.clip(p_k, 1e-4, 1.0 - 1e-4)
        std_k = np.sqrt(p_k * (1.0 - p_k))
        W = (X - p_k) / std_k

        pca = PCA(n_components=n_components, random_state=42)
        pcs = pca.fit_transform(W)
        self.pca_model = pca

        # 1. Eigenvalues and Variance Explained table
        eval_rows = []
        cum_var = 0.0
        for i, (eig_val, var_pct) in enumerate(zip(pca.explained_variance_, pca.explained_variance_ratio_)):
            cum_var += var_pct * 100.0
            eval_rows.append({
                "Component": f"PC{i + 1}",
                "Eigenvalue": float(eig_val),
                "Variance_Explained_Pct": round(var_pct * 100.0, 3),
                "Cumulative_Variance_Pct": round(cum_var, 3),
            })
        self.eigenvalues_df = pd.DataFrame(eval_rows)

        # 2. PC coordinates table
        coord_dict = {
            "Sample_ID": self.sample_metadata["Sample_ID"],
            "Cultivar": self.sample_metadata["Cultivar"],
            "Taxon": self.sample_metadata["Taxon"],
            "Varietal_Group": self.sample_metadata["Varietal_Group"],
            "Nation_Origin": self.sample_metadata["Nation_Origin"],
            "Region": self.sample_metadata["Region"],
        }
        for i in range(n_components):
            coord_dict[f"PC{i + 1}"] = pcs[:, i]

        self.pca_coords_df = pd.DataFrame(coord_dict)

        # 3. SNP Loadings on top 3 PCs
        loadings = pca.components_[:3, :].T  # (M_pruned, 3)
        load_dict = {
            "POS": self.pruned_variants_df["POS"],
            "REF": self.pruned_variants_df["REF"],
            "ALT": self.pruned_variants_df["ALT"],
            "DOMAIN": self.pruned_variants_df["DOMAIN"],
            "GENE": self.pruned_variants_df["GENE"],
            "COMPLEX": self.pruned_variants_df["COMPLEX"],
            "PC1_Loading": loadings[:, 0],
            "PC2_Loading": loadings[:, 1],
            "PC3_Loading": loadings[:, 2],
        }
        self.loadings_df = pd.DataFrame(load_dict)

        print(f"    PCA completed in {time.time() - t0:.2f}s.")
        print(f"    PC1 explains {self.eigenvalues_df.iloc[0]['Variance_Explained_Pct']}% variance.")
        print(f"    PC2 explains {self.eigenvalues_df.iloc[1]['Variance_Explained_Pct']}% variance.")
        print(f"    PC3 explains {self.eigenvalues_df.iloc[2]['Variance_Explained_Pct']}% variance.")
        print(f"    Top 3 PCs explain cumulative {self.eigenvalues_df.iloc[2]['Cumulative_Variance_Pct']}% variance.")

        return self.pca_coords_df, self.eigenvalues_df, self.loadings_df

    # ==========================================================================
    # 6. UMAP Manifold Embeddings (2D & 3D)
    # ==========================================================================

    def perform_umap(self) -> pd.DataFrame:
        """Compute 2D and 3D UMAP manifold embeddings."""
        print("--> Computing 2D and 3D UMAP manifold projections...")
        t0 = time.time()

        # Standardize matrix W
        X = self.pruned_gt_matrix.T.astype(np.float64)
        p_k = np.mean(X == 1, axis=0)
        p_k = np.clip(p_k, 1e-4, 1.0 - 1e-4)
        std_k = np.sqrt(p_k * (1.0 - p_k))
        W = (X - p_k) / std_k

        # 2D UMAP
        u2 = umap.UMAP(
            n_neighbors=15,
            min_dist=0.1,
            metric="euclidean",
            n_components=2,
            random_state=42,
            n_jobs=1,
        ).fit_transform(W)

        # 3D UMAP
        u3 = umap.UMAP(
            n_neighbors=15,
            min_dist=0.1,
            metric="euclidean",
            n_components=3,
            random_state=42,
            n_jobs=1,
        ).fit_transform(W)

        umap_dict = {
            "Sample_ID": self.sample_metadata["Sample_ID"],
            "Cultivar": self.sample_metadata["Cultivar"],
            "Varietal_Group": self.sample_metadata["Varietal_Group"],
            "Nation_Origin": self.sample_metadata["Nation_Origin"],
            "Region": self.sample_metadata["Region"],
            "UMAP1_2D": u2[:, 0],
            "UMAP2_2D": u2[:, 1],
            "UMAP1_3D": u3[:, 0],
            "UMAP2_3D": u3[:, 1],
            "UMAP3_3D": u3[:, 2],
        }
        self.umap_coords_df = pd.DataFrame(umap_dict)

        print(f"    UMAP manifold projections computed in {time.time() - t0:.2f}s.")
        return self.umap_coords_df

    # ==========================================================================
    # 7. Deliverable Exports
    # ==========================================================================

    def export_tables(self) -> None:
        """Export all tabular results in TSV and CSV formats."""
        os.makedirs(self.out_table_dir, exist_ok=True)

        # 1. Pruned SNPs table
        self.pruned_variants_df.to_csv(f"{self.out_table_dir}/plastome_snps_ld_pruned.tsv", sep="\t", index=False)
        self.pruned_variants_df.to_csv(f"{self.out_table_dir}/plastome_snps_ld_pruned.csv", index=False)

        # 2. Kinship Matrix
        df_k = pd.DataFrame(self.kinship_matrix, index=self.samples, columns=self.samples)
        df_k.to_csv(f"{self.out_table_dir}/plastome_kinship_matrix.tsv", sep="\t")
        df_k.to_csv(f"{self.out_table_dir}/plastome_kinship_matrix.csv")

        # 3. PCA tables
        self.eigenvalues_df.to_csv(f"{self.out_table_dir}/plastome_pca_eigenvalues.tsv", sep="\t", index=False)
        self.eigenvalues_df.to_csv(f"{self.out_table_dir}/plastome_pca_eigenvalues.csv", index=False)

        self.pca_coords_df.to_csv(f"{self.out_table_dir}/plastome_pca_coordinates.tsv", sep="\t", index=False)
        self.pca_coords_df.to_csv(f"{self.out_table_dir}/plastome_pca_coordinates.csv", index=False)

        self.loadings_df.to_csv(f"{self.out_table_dir}/plastome_pca_loadings.tsv", sep="\t", index=False)
        self.loadings_df.to_csv(f"{self.out_table_dir}/plastome_pca_loadings.csv", index=False)

        # 4. UMAP tables
        self.umap_coords_df.to_csv(f"{self.out_table_dir}/plastome_umap_coordinates.tsv", sep="\t", index=False)
        self.umap_coords_df.to_csv(f"{self.out_table_dir}/plastome_umap_coordinates.csv", index=False)

        print(f"--> Exported: {self.out_table_dir}/plastome_snps_ld_pruned.[tsv/csv] ({len(self.pruned_variants_df)} SNPs)")
        print(f"--> Exported: {self.out_table_dir}/plastome_kinship_matrix.[tsv/csv] (90x90 matrix)")
        print(f"--> Exported: {self.out_table_dir}/plastome_pca_eigenvalues.[tsv/csv]")
        print(f"--> Exported: {self.out_table_dir}/plastome_pca_coordinates.[tsv/csv]")
        print(f"--> Exported: {self.out_table_dir}/plastome_pca_loadings.[tsv/csv]")
        print(f"--> Exported: {self.out_table_dir}/plastome_umap_coordinates.[tsv/csv]")

    # ==========================================================================
    # 8. Publication-Grade Visualizations
    # ==========================================================================

    def plot_pca_figures(self) -> None:
        """Generate Figure 7: Publication-grade 4-panel PCA composite."""
        print("--> Generating publication-grade PCA figure (Figure 7)...")
        os.makedirs(self.out_fig_dir, exist_ok=True)

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # Color palettes
        group_palette = {
            "Japonica": "#2b5c8f",
            "Indica": "#d95f02",
            "Spontanea": "#e6ab02",
            "Aromatic": "#7570b3",
            "Other_Sativa": "#555555",
        }
        region_palette = {
            "Northern_Margin_Russia": "#1b9e77",
            "East_Asia": "#d95f02",
            "South_Southeast_Asia": "#7570b3",
            "Middle_East": "#e7298a",
            "Other_Unspecified": "#666666",
        }

        # ----------------------------------------------------------------------
        # Panel A: Scree Plot (Eigenvalues & Cumulative Variance)
        # ----------------------------------------------------------------------
        ax_a = axes[0, 0]
        pcs = np.arange(1, len(self.eigenvalues_df) + 1)
        var_pcts = self.eigenvalues_df["Variance_Explained_Pct"].values
        cum_pcts = self.eigenvalues_df["Cumulative_Variance_Pct"].values

        ax_a.bar(pcs, var_pcts, color="#2b5c8f", edgecolor="black", alpha=0.85, label="Individual Variance (%)")
        ax_a_twin = ax_a.twinx()
        ax_a_twin.plot(pcs, cum_pcts, color="#d95f02", marker="o", linewidth=2.0, label="Cumulative Variance (%)")

        ax_a.set_xlabel("Principal Component (PC)")
        ax_a.set_ylabel("Variance Explained (%)", color="#2b5c8f")
        ax_a_twin.set_ylabel("Cumulative Variance (%)", color="#d95f02")
        ax_a.set_title("A. Scree Plot of Eigenvalues & Cumulative Variance Explained", weight="bold")
        ax_a.set_xticks(pcs)
        ax_a.grid(axis="y", linestyle=":", alpha=0.6)

        # ----------------------------------------------------------------------
        # Panel B: PC1 vs PC2 colored by Varietal Group with 95% Confidence Ellipses
        # ----------------------------------------------------------------------
        ax_b = axes[0, 1]
        pc1_pct = self.eigenvalues_df.iloc[0]["Variance_Explained_Pct"]
        pc2_pct = self.eigenvalues_df.iloc[1]["Variance_Explained_Pct"]

        for grp, color in group_palette.items():
            sub = self.pca_coords_df[self.pca_coords_df["Varietal_Group"] == grp]
            if len(sub) == 0:
                continue
            ax_b.scatter(
                sub["PC1"],
                sub["PC2"],
                c=color,
                label=f"{grp} (n={len(sub)})",
                s=60,
                alpha=0.85,
                edgecolors="black",
                linewidth=0.8,
            )
            # Add confidence ellipse if n >= 4
            if len(sub) >= 4:
                cov = np.cov(sub["PC1"], sub["PC2"])
                lambda_, v = np.linalg.eig(cov)
                lambda_ = np.sqrt(np.maximum(lambda_, 0))
                ell = Ellipse(
                    xy=(np.mean(sub["PC1"]), np.mean(sub["PC2"])),
                    width=lambda_[0] * 2 * 2.0,
                    height=lambda_[1] * 2 * 2.0,
                    angle=np.rad2deg(np.arccos(v[0, 0])),
                    edgecolor=color,
                    facecolor=color,
                    alpha=0.15,
                    linestyle="--",
                    linewidth=1.2,
                )
                ax_b.add_patch(ell)

        ax_b.axhline(0, color="gray", linestyle=":", linewidth=0.8)
        ax_b.axvline(0, color="gray", linestyle=":", linewidth=0.8)
        ax_b.set_xlabel(f"PC1 ({pc1_pct:.1f}% Variance Explained)")
        ax_b.set_ylabel(f"PC2 ({pc2_pct:.1f}% Variance Explained)")
        ax_b.set_title("B. Population Structure: PC1 vs PC2 by Varietal Group", weight="bold")
        ax_b.legend(loc="upper left", frameon=True, fontsize=8)
        ax_b.grid(True, linestyle=":", alpha=0.5)

        # ----------------------------------------------------------------------
        # Panel C: PC1 vs PC2 colored by Geographic Provenance
        # ----------------------------------------------------------------------
        ax_c = axes[1, 0]
        for reg, color in region_palette.items():
            sub = self.pca_coords_df[self.pca_coords_df["Region"] == reg]
            if len(sub) == 0:
                continue
            ax_c.scatter(
                sub["PC1"],
                sub["PC2"],
                c=color,
                label=f"{reg.replace('_', ' ')} (n={len(sub)})",
                s=60,
                alpha=0.85,
                edgecolors="black",
                linewidth=0.8,
            )

        ax_c.axhline(0, color="gray", linestyle=":", linewidth=0.8)
        ax_c.axvline(0, color="gray", linestyle=":", linewidth=0.8)
        ax_c.set_xlabel(f"PC1 ({pc1_pct:.1f}% Variance Explained)")
        ax_c.set_ylabel(f"PC2 ({pc2_pct:.1f}% Variance Explained)")
        ax_c.set_title("C. Geographic Stratification: PC1 vs PC2 by Agro-Zone", weight="bold")
        ax_c.legend(loc="upper left", frameon=True, fontsize=8)
        ax_c.grid(True, linestyle=":", alpha=0.5)

        # ----------------------------------------------------------------------
        # Panel D: PC1 vs PC3 resolving Wild Progenitors & Subgroups
        # ----------------------------------------------------------------------
        ax_d = axes[1, 1]
        pc3_pct = self.eigenvalues_df.iloc[2]["Variance_Explained_Pct"]

        for grp, color in group_palette.items():
            sub = self.pca_coords_df[self.pca_coords_df["Varietal_Group"] == grp]
            if len(sub) == 0:
                continue
            ax_d.scatter(
                sub["PC1"],
                sub["PC3"],
                c=color,
                label=grp,
                s=60,
                alpha=0.85,
                edgecolors="black",
                linewidth=0.8,
            )

        # Annotate prominent anchor lines using adjust_text
        texts = []
        anchor_cultivars = ["Nipponbare", "Dular", "Zhenshan 97", "Minghui 63", "Saw-An-T", "Tai Bonnet"]
        for _, r in self.pca_coords_df.iterrows():
            if any(a.lower() in r["Cultivar"].lower() for a in anchor_cultivars):
                txt = ax_d.text(
                    r["PC1"],
                    r["PC3"],
                    r["Cultivar"],
                    fontsize=8,
                    weight="bold",
                    color="#111111",
                )
                texts.append(txt)

        if texts:
            adjust_text(texts, arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8), ax=ax_d)

        ax_d.axhline(0, color="gray", linestyle=":", linewidth=0.8)
        ax_d.axvline(0, color="gray", linestyle=":", linewidth=0.8)
        ax_d.set_xlabel(f"PC1 ({pc1_pct:.1f}% Variance Explained)")
        ax_d.set_ylabel(f"PC3 ({pc3_pct:.1f}% Variance Explained)")
        ax_d.set_title("D. Fine-Scale Differentiation: PC1 vs PC3 with Anchor Cultivars", weight="bold")
        ax_d.legend(loc="upper left", frameon=True, fontsize=8)
        ax_d.grid(True, linestyle=":", alpha=0.5)

        plt.tight_layout()
        fig.savefig(f"{self.out_fig_dir}/fig7_plastome_pca_structure.pdf")
        fig.savefig(f"{self.out_fig_dir}/fig7_plastome_pca_structure.png")
        plt.close(fig)
        print("    Saved Figure 7: fig7_plastome_pca_structure.[pdf/png]")

    def plot_umap_and_kinship(self) -> None:
        """Generate Figure 8: 2D/3D UMAP Manifold & Hierarchical Kinship Heatmap."""
        print("--> Generating publication-grade UMAP & Kinship figure (Figure 8)...")
        os.makedirs(self.out_fig_dir, exist_ok=True)

        fig = plt.figure(figsize=(16, 13))

        group_palette = {
            "Japonica": "#2b5c8f",
            "Indica": "#d95f02",
            "Spontanea": "#e6ab02",
            "Aromatic": "#7570b3",
            "Other_Sativa": "#555555",
        }
        region_palette = {
            "Northern_Margin_Russia": "#1b9e77",
            "East_Asia": "#d95f02",
            "South_Southeast_Asia": "#7570b3",
            "Middle_East": "#e7298a",
            "Other_Unspecified": "#666666",
        }

        # ----------------------------------------------------------------------
        # Panel A: 2D UMAP by Varietal Group
        # ----------------------------------------------------------------------
        ax_a = fig.add_subplot(2, 2, 1)
        for grp, color in group_palette.items():
            sub = self.umap_coords_df[self.umap_coords_df["Varietal_Group"] == grp]
            if len(sub) == 0:
                continue
            ax_a.scatter(
                sub["UMAP1_2D"],
                sub["UMAP2_2D"],
                c=color,
                label=f"{grp} (n={len(sub)})",
                s=70,
                alpha=0.88,
                edgecolors="black",
                linewidth=0.8,
            )

        # Label key anchor cultivars
        texts_a = []
        anchor_cultivars = ["Nipponbare", "Dular", "Zhenshan 97", "Minghui 63", "Saw-An-T", "Tai Bonnet"]
        for _, r in self.umap_coords_df.iterrows():
            if any(a.lower() in r["Cultivar"].lower() for a in anchor_cultivars):
                txt = ax_a.text(
                    r["UMAP1_2D"],
                    r["UMAP2_2D"],
                    r["Cultivar"],
                    fontsize=8,
                    weight="bold",
                    color="#111111",
                )
                texts_a.append(txt)

        if texts_a:
            adjust_text(texts_a, arrowprops=dict(arrowstyle="->", color="#333333", lw=0.8), ax=ax_a)

        ax_a.set_xlabel("UMAP 1")
        ax_a.set_ylabel("UMAP 2")
        ax_a.set_title("A. 2D UMAP Manifold Embedding by Varietal Group", weight="bold")
        ax_a.legend(loc="lower left", frameon=True, fontsize=8)
        ax_a.grid(True, linestyle=":", alpha=0.5)

        # ----------------------------------------------------------------------
        # Panel B: 2D UMAP by Geographic Provenance
        # ----------------------------------------------------------------------
        ax_b = fig.add_subplot(2, 2, 2)
        for reg, color in region_palette.items():
            sub = self.umap_coords_df[self.umap_coords_df["Region"] == reg]
            if len(sub) == 0:
                continue
            ax_b.scatter(
                sub["UMAP1_2D"],
                sub["UMAP2_2D"],
                c=color,
                label=f"{reg.replace('_', ' ')} (n={len(sub)})",
                s=70,
                alpha=0.88,
                edgecolors="black",
                linewidth=0.8,
            )

        ax_b.set_xlabel("UMAP 1")
        ax_b.set_ylabel("UMAP 2")
        ax_b.set_title("B. 2D UMAP Manifold Embedding by Geographic Provenance", weight="bold")
        ax_b.legend(loc="lower left", frameon=True, fontsize=8)
        ax_b.grid(True, linestyle=":", alpha=0.5)

        # ----------------------------------------------------------------------
        # Panel C: 3D UMAP Manifold Isometric Projection
        # ----------------------------------------------------------------------
        ax_c = fig.add_subplot(2, 2, 3, projection="3d")
        for grp, color in group_palette.items():
            sub = self.umap_coords_df[self.umap_coords_df["Varietal_Group"] == grp]
            if len(sub) == 0:
                continue
            ax_c.scatter(
                sub["UMAP1_3D"],
                sub["UMAP2_3D"],
                sub["UMAP3_3D"],
                c=color,
                label=grp,
                s=60,
                alpha=0.85,
                edgecolors="black",
                linewidth=0.7,
            )

        ax_c.set_xlabel("UMAP 1", fontsize=8, labelpad=2)
        ax_c.set_ylabel("UMAP 2", fontsize=8, labelpad=2)
        ax_c.set_zlabel("UMAP 3", fontsize=8, labelpad=2)
        ax_c.set_title("C. 3D UMAP Isometric Manifold Projection", weight="bold", pad=10)
        ax_c.view_init(elev=25, azim=45)
        ax_c.legend(loc="upper left", frameon=True, fontsize=7)

        # ----------------------------------------------------------------------
        # Panel D: Hierarchically Clustered Kinship Matrix Heatmap
        # ----------------------------------------------------------------------
        ax_d = fig.add_subplot(2, 2, 4)

        # Build colored row annotation
        group_colors = [group_palette.get(g, "#555555") for g in self.sample_metadata["Varietal_Group"]]
        
        # Sort samples by Varietal Group for clear block structure
        sort_idx = np.lexsort((self.sample_metadata["Region"].values, self.sample_metadata["Varietal_Group"].values))
        sorted_K = self.kinship_matrix[sort_idx, :][:, sort_idx]

        im = ax_d.imshow(sorted_K, cmap="mako", vmin=-0.5, vmax=2.5, aspect="auto")
        cbar = fig.colorbar(im, ax=ax_d, fraction=0.046, pad=0.04)
        cbar.set_label("Realized Genomic Kinship (K)")

        # Add divider lines between groups
        sorted_groups = self.sample_metadata.iloc[sort_idx]["Varietal_Group"].values
        boundaries = np.where(sorted_groups[:-1] != sorted_groups[1:])[0]
        for b in boundaries:
            ax_d.axhline(b + 0.5, color="white", linewidth=1.2)
            ax_d.axvline(b + 0.5, color="white", linewidth=1.2)

        ax_d.set_title("D. Realized Genomic Kinship Matrix (K, 90x90 Accessions)", weight="bold")
        ax_d.set_xlabel("Accessions (Organized by Varietal Group)")
        ax_d.set_ylabel("Accessions")
        ax_d.set_xticks([])
        ax_d.set_yticks([])

        plt.tight_layout()
        fig.savefig(f"{self.out_fig_dir}/fig8_plastome_umap_kinship.pdf")
        fig.savefig(f"{self.out_fig_dir}/fig8_plastome_umap_kinship.png")
        plt.close(fig)
        print("    Saved Figure 8: fig8_plastome_umap_kinship.[pdf/png]")

    # ==========================================================================
    # 9. Master Runner
    # ==========================================================================

    def run_pipeline(self) -> None:
        """Execute full Subtask 3.3 population structure analysis."""
        t0 = time.time()
        print("\n" + "=" * 80)
        print("Starting Subtask 3.3: Population Structure, Kinship & Manifold Engine")
        print("=" * 80 + "\n")

        self.load_data()
        self.perform_ld_pruning(min_maf=0.02, r2_threshold=0.2, window_bp=20000)
        self.compute_kinship_matrix()
        self.perform_pca(n_components=10)
        self.perform_umap()
        self.export_tables()

        self.plot_pca_figures()
        self.plot_umap_and_kinship()

        print("\n" + "=" * 80)
        print(f"Subtask 3.3 Pipeline Completed Successfully in {time.time() - t0:.2f} seconds!")
        print("=" * 80 + "\n")


if __name__ == "__main__":
    engine = PlastomePopulationStructureEngine()
    engine.run_pipeline()
