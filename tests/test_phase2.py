"""
Unit tests for Phase 2: Sequence Alignment, Annotation & Variant Calling
Verifies dataset completeness, cryptographic manifests, coordinate liftover integrity,
callable site denominator accounting, and polarized VCF 4.2 conformance.
"""

import os
import gzip
import hashlib
import pytest
import pandas as pd
from Bio import SeqIO


def compute_sha256(filepath: str) -> str:
    """Compute SHA-256 cryptographic digest of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def test_plastomes_standardization_manifest():
    """Verify cryptographic integrity of all standardized plastome fasta files."""
    manifest_path = "data/interim/plastomes_standardized/manifest.sha256"
    assert os.path.exists(manifest_path), "Standardized plastomes manifest missing"
    
    with open(manifest_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
        
    assert len(lines) >= 110, f"Expected at least 110 entries in manifest, got {len(lines)}"
    
    for line in lines:
        expected_hash, rel_path = line.split("  ")
        full_path = os.path.join("data/interim/plastomes_standardized", rel_path)
        assert os.path.exists(full_path), f"File {full_path} missing"
        actual_hash = compute_sha256(full_path)
        assert actual_hash == expected_hash, f"SHA-256 mismatch for {rel_path}"


def test_plastome_standardization_metrics():
    """Verify standardization summary table and sequence start homology."""
    tsv_path = "results/tables/plastome_standardization_summary.tsv"
    assert os.path.exists(tsv_path), "Plastome standardization summary TSV missing"
    
    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) == 109, f"Expected 109 accessions in summary table, got {len(df)}"
    
    # Check that all sequences have high 5' homology to reference start
    assert (df["Start_Homology_25bp"] >= 20).all(), "Some accessions failed start homology check"
    
    # Check that known inverted plastomes CP054688.1 and CP056064.1 were reoriented
    cp1 = df[df["Accession"] == "CP054688.1"]
    cp2 = df[df["Accession"] == "CP056064.1"]
    assert not cp1.empty and cp1.iloc[0]["Reversed"], "CP054688.1 was not properly reversed"
    assert not cp2.empty and cp2.iloc[0]["Reversed"], "CP056064.1 was not properly reversed"


def test_plastome_annotations():
    """Verify GFF3, BED tracks, and gene annotation tables."""
    gff_path = "data/interim/annotation/rice_chloroplast_NC_001320.1.gff3"
    assert os.path.exists(gff_path), "GFF3 annotation missing"
    with open(gff_path, "r", encoding="utf-8") as f:
        first_line = f.readline().strip()
    assert first_line == "##gff-version 3", f"Invalid GFF version: {first_line}"
    
    gene_tsv = "results/tables/plastome_gene_annotations.tsv"
    assert os.path.exists(gene_tsv), "Gene annotation TSV missing"
    df_genes = pd.read_csv(gene_tsv, sep="\t")
    assert len(df_genes) >= 100, f"Expected >= 100 CDS features, got {len(df_genes)}"
    
    required_cols = ["Gene_Symbol", "Start_1based", "End_1based", "Functional_Complex", "Translation"]
    for col in required_cols:
        assert col in df_genes.columns, f"Required column {col} missing in gene annotations"
        
    # Check key adaptive genes
    genes = set(df_genes["Gene_Symbol"].unique())
    assert "rbcL" in genes, "Rubisco large subunit (rbcL) missing"
    assert "psbA" in genes, "Photosystem II reaction center (psbA) missing"
    assert "matK" in genes, "Maturase K (matK) missing"
    assert any("ndh" in g.lower() for g in genes), "NADH dehydrogenase (ndh) family missing"
    
    # Check quadripartite structure table
    quad_tsv = "results/tables/plastome_quadripartite_structure.tsv"
    assert os.path.exists(quad_tsv), "Quadripartite structure TSV missing"
    df_quad = pd.read_csv(quad_tsv, sep="\t")
    assert set(df_quad["Domain"]) == {"LSC", "IRb", "SSC", "IRa"}
    assert df_quad["Length_bp"].sum() == 134525, "Quadripartite domains do not sum to 134,525 bp"


def test_callable_sites_mask():
    """Verify callable bed mask invariants and exact denominator accounting."""
    callable_bed = "data/interim/callable_mask/callable_sites.bed"
    masked_bed = "data/interim/callable_mask/masked_regions.bed"
    summary_tsv = "results/tables/callable_sites_summary.tsv"
    
    assert os.path.exists(callable_bed), "callable_sites.bed missing"
    assert os.path.exists(masked_bed), "masked_regions.bed missing"
    assert os.path.exists(summary_tsv), "callable_sites_summary.tsv missing"
    
    df_sum = pd.read_csv(summary_tsv, sep="\t")
    whole_row = df_sum[df_sum["Genomic_Partition"] == "Whole_Plastome"].iloc[0]
    
    assert whole_row["Total_Span_bp"] == 134525
    assert whole_row["L_callable_bp"] > 130000, f"Callable span too low: {whole_row['L_callable_bp']}"
    assert whole_row["L_callable_bp"] + whole_row["Masked_Span_bp"] == 134525, "Callable + Masked mismatch"
    
    # Verify non-empty per-gene callable table
    gene_call_tsv = "results/tables/gene_callable_lengths.tsv"
    assert os.path.exists(gene_call_tsv), "gene_callable_lengths.tsv missing"
    df_gene_call = pd.read_csv(gene_call_tsv, sep="\t")
    assert len(df_gene_call) >= 100
    assert (df_gene_call["L_callable_bp"] > 0).all(), "Some genes have zero callable sites"


def test_msa_alignments():
    """Verify multi-FASTA alignments and manifest cryptographic checksums."""
    manifest_path = "results/alignments/manifest.sha256"
    assert os.path.exists(manifest_path), "Alignments manifest missing"
    
    with open(manifest_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
    assert len(lines) >= 3, "Expected at least 3 alignment files in manifest"
    
    for line in lines:
        exp_hash, fname = line.split("  ")
        fpath = os.path.join("results/alignments", fname)
        assert os.path.exists(fpath), f"File {fpath} does not exist"
        act_hash = compute_sha256(fpath)
        assert act_hash == exp_hash, f"SHA-256 mismatch for {fname}"
        
    # Check all_plastomes_aligned.fasta
    all_fa = "results/alignments/all_plastomes_aligned.fasta"
    recs_all = list(SeqIO.parse(all_fa, "fasta"))
    assert len(recs_all) == 109, f"Expected 109 taxa, got {len(recs_all)}"
    aln_len = len(recs_all[0].seq)
    assert all(len(r.seq) == aln_len for r in recs_all), "Inconsistent alignment lengths across taxa"
    assert aln_len >= 134525, f"Alignment length {aln_len} less than reference span"
    
    # Check sativa_plastomes_aligned.fasta
    sat_fa = "results/alignments/sativa_plastomes_aligned.fasta"
    recs_sat = list(SeqIO.parse(sat_fa, "fasta"))
    assert len(recs_sat) == 90, f"Expected 90 sativa taxa, got {len(recs_sat)}"
    assert all(len(r.seq) == aln_len for r in recs_sat)
    
    # Check single-copy partitioned alignment
    sc_fa = "results/alignments/all_plastomes_single_copy.fasta"
    recs_sc = list(SeqIO.parse(sc_fa, "fasta"))
    assert len(recs_sc) == 109
    sc_len = len(recs_sc[0].seq)
    assert sc_len < aln_len, "Single-copy alignment should be shorter than full alignment (IRa masked)"
    assert all(len(r.seq) == sc_len for r in recs_sc)


def test_coordinate_liftover_index():
    """Verify coordinate liftover table maps 100% of canonical reference coordinates."""
    liftover_path = "results/tables/coordinate_liftover_index.tsv.gz"
    assert os.path.exists(liftover_path), "coordinate_liftover_index.tsv.gz missing"
    
    with gzip.open(liftover_path, "rt", encoding="utf-8") as f:
        df = pd.read_csv(f, sep="\t")
        
    ref_coords = set(df[df["Reference_Coordinate"] > 0]["Reference_Coordinate"])
    assert len(ref_coords) == 134525, f"Expected 134,525 unique reference coordinates, got {len(ref_coords)}"
    assert min(ref_coords) == 1
    assert max(ref_coords) == 134525
    assert set(df["Domain"].unique()) == {"LSC", "IRb", "SSC", "IRa"}


def test_polarized_vcf():
    """Verify VCF 4.2 format conformance, provenance headers, and ancestral polarization."""
    vcf_path = "data/processed/variants/rice_plastome_variants_polarized.vcf"
    vcf_gz_path = "data/processed/variants/rice_plastome_variants_polarized.vcf.gz"
    manifest_path = "data/processed/variants/manifest.sha256"
    
    assert os.path.exists(vcf_path), "VCF file missing"
    assert os.path.exists(vcf_gz_path), "Compressed VCF file missing"
    assert os.path.exists(manifest_path), "Variants manifest missing"
    
    with open(manifest_path, "r", encoding="utf-8") as f:
        for line in f:
            exp_hash, fname = line.strip().split("  ")
            act_hash = compute_sha256(os.path.join("data/processed/variants", fname))
            assert act_hash == exp_hash, f"SHA-256 mismatch for {fname}"
            
    # Check VCF headers and records
    header_lines = []
    record_lines = []
    with open(vcf_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                header_lines.append(line.strip())
            else:
                record_lines.append(line.strip())
                
    assert any("##fileformat=VCFv4.2" in h for h in header_lines)
    assert any("##INFO=<ID=AA" in h for h in header_lines), "AA (Ancestral Allele) INFO header missing"
    assert any("##INFO=<ID=DAF" in h for h in header_lines), "DAF (Derived Allele Frequency) INFO header missing"
    assert any("##INFO=<ID=EFF" in h for h in header_lines), "EFF (Functional Effect) INFO header missing"
    assert any("##agent_command" in h for h in header_lines), "agent_command provenance header missing"
    
    # Check column header
    col_header = [h for h in header_lines if h.startswith("#CHROM")][0].split("\t")
    assert col_header[:9] == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    # 90 Sativa accessions
    assert len(col_header[9:]) == 90
    
    assert len(record_lines) >= 5000, f"Expected >= 5,000 SNPs, got {len(record_lines)}"
    
    # Check sample records
    aa_polarized_count = 0
    for line in record_lines[:500]:
        fields = line.split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        info = fields[7]
        
        assert chrom == "NC_001320.1"
        assert 1 <= pos <= 134525
        assert len(ref) >= 1
        assert len(alt) >= 1
        
        if "AA=" in info and "AA=." not in info:
            aa_polarized_count += 1
            
    assert aa_polarized_count >= 400, "Insufficient polarized ancestral alleles in VCF sample"


def test_diagnostic_indels():
    """Verify diagnostic indels table and detection of the canonical ~50 bp length dimorphism."""
    indel_tsv = "results/tables/plastome_indels_diagnostic.tsv"
    assert os.path.exists(indel_tsv), "Diagnostic indels TSV missing"
    
    df_indels = pd.read_csv(indel_tsv, sep="\t")
    assert len(df_indels) >= 10, f"Expected >= 10 polymorphic indels, got {len(df_indels)}"
    
    # Check that large diagnostic indels (>= 20 bp) are identified
    large_indels = df_indels[df_indels["Is_Diagnostic_Marker"]]
    assert len(large_indels) >= 1, "No large diagnostic indels (>= 20 bp) identified"
    
    # Verify ~50 bp dimorphism presence (length between 40 and 65 bp)
    has_50bp = (large_indels["Length_bp"] >= 40) & (large_indels["Length_bp"] <= 65)
    assert has_50bp.any(), "Canonical ~50 bp length dimorphism not found among diagnostic indels"
