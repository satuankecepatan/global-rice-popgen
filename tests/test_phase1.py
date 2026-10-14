"""
Unit tests for Phase 1: Data Acquisition & Quality Curation
Verifies dataset completeness, sequence integrity, SHA-256 manifest cryptographic consistency,
and metadata schema adherence.
"""

import os
import hashlib
import pytest
import pandas as pd
from Bio import SeqIO

def compute_sha256(filepath: str) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def test_sativa_plastomes_manifest():
    manifest_path = "data/raw/plastomes/sativa/manifest.sha256"
    assert os.path.exists(manifest_path), "Sativa manifest.sha256 missing"
    
    with open(manifest_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
        
    assert len(lines) >= 90, f"Expected at least 90 manifest entries, got {len(lines)}"
    
    # Verify cryptographic integrity of all files in manifest
    for line in lines:
        expected_hash, fname = line.split("  ")
        fpath = os.path.join("data/raw/plastomes/sativa", fname)
        assert os.path.exists(fpath), f"File {fpath} does not exist"
        actual_hash = compute_sha256(fpath)
        assert actual_hash == expected_hash, f"SHA-256 mismatch for {fname}"

def test_sativa_sequence_properties():
    combined_fa = "data/raw/plastomes/sativa/all_sativa_plastomes.fasta"
    assert os.path.exists(combined_fa), "all_sativa_plastomes.fasta missing"
    
    records = list(SeqIO.parse(combined_fa, "fasta"))
    assert len(records) == 90, f"Expected 90 sativa plastomes, got {len(records)}"
    
    for rec in records:
        seq_len = len(rec.seq)
        assert 134_000 <= seq_len <= 135_500, f"Plastome {rec.id} length {seq_len} bp out of bounds"
        # Validate minimal ambiguous characters (< 0.1%)
        ambig_count = rec.seq.upper().count("N")
        assert ambig_count / seq_len < 0.01, f"Plastome {rec.id} has excessive ambiguous bases ({ambig_count})"

def test_wild_outgroups_manifest():
    manifest_path = "data/raw/plastomes/outgroups/manifest.sha256"
    assert os.path.exists(manifest_path), "Outgroups manifest.sha256 missing"
    
    with open(manifest_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
        
    assert len(lines) >= 19, f"Expected at least 19 manifest entries, got {len(lines)}"
    
    for line in lines:
        expected_hash, fname = line.split("  ")
        fpath = os.path.join("data/raw/plastomes/outgroups", fname)
        assert os.path.exists(fpath), f"File {fpath} does not exist"
        actual_hash = compute_sha256(fpath)
        assert actual_hash == expected_hash, f"SHA-256 mismatch for {fname}"

def test_wild_outgroups_metadata():
    tsv_path = "results/tables/ncbi_wild_outgroup_chloroplast_genomes.tsv"
    assert os.path.exists(tsv_path), "Wild outgroups TSV table missing"
    
    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) == 19, f"Expected 19 outgroup accessions, got {len(df)}"
    
    species = set(df["Species"].unique())
    assert "Oryza rufipogon" in species
    assert "Oryza nivara" in species
    assert "Oryza punctata" in species

def test_3k_passport_metadata():
    tsv_path = "results/tables/rice_3k_passport_metadata.tsv"
    assert os.path.exists(tsv_path), "3K passport metadata TSV missing"
    
    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) == 3024, f"Expected 3,024 accessions in passport, got {len(df)}"
    
    required_cols = [
        "Accession_ID_3K", "Variety_Name", "Country_Origin", "Subpopulation_Full",
        "ind1A", "ind1B", "ind2", "ind3", "aus", "temp", "trop1", "trop2", "aro"
    ]
    for col in required_cols:
        assert col in df.columns, f"Required column {col} missing in passport metadata"
        
    # Check that Q-matrix sums to approximately 1.0 (for non-NaN admixture entries)
    q_cols = ["ind1A", "ind1B", "ind2", "ind3", "aus", "temp", "trop1", "trop2", "aro"]
    q_sums = df[q_cols].dropna().sum(axis=1)
    assert ((q_sums >= 0.95) & (q_sums <= 1.05)).all(), "Admixture proportions do not sum to 1.0"

def test_3k_phenotype_metadata():
    tsv_path = "results/tables/rice_3k_phenotype_metadata.tsv"
    assert os.path.exists(tsv_path), "3K phenotype metadata TSV missing"
    
    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) == 2266, f"Expected 2,266 evaluated accessions in phenotype table, got {len(df)}"
    assert "Cultivar_IRGC_Name" in df.columns
    assert "Days_To_Heading_Flowering" in df.columns

def test_plastome_nuclear_crossref():
    tsv_path = "results/tables/plastome_nuclear_crossref.tsv"
    assert os.path.exists(tsv_path), "Plastome-nuclear cross-reference TSV missing"
    
    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) >= 89, f"Expected at least 89 accessions, got {len(df)}"
    assert "Plastome_Accession" in df.columns
    assert "Cultivar" in df.columns
    assert "Concordance_Status" in df.columns
