"""
Standardize Plastomes Pipeline
==============================
Performs orientation detection, circular linearization, and start-coordinate
standardization for all Oryza sativa and wild outgroup chloroplast genomes against
the canonical Nipponbare reference standard (NC_001320.1).

Ensures all sequences begin at the exact homologous 5'-terminal locus (coordinate 1)
and reside on the positive (+) forward strand, preventing circular permutation
artifacts during downstream multiple sequence alignment and variant calling.
"""

import os
import glob
import hashlib
from typing import Tuple, List, Dict
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def compute_sha256(filepath: str) -> str:
    """Compute SHA-256 cryptographic digest of a file."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def detect_strand(seq: str, ref_seq: str, kmer_len: int = 31, step: int = 5000) -> Tuple[str, bool]:
    """
    Detect whether sequence is on forward (+) or reverse-complement (-) strand
    relative to reference genome by multi-kmer voting across the chromosome.
    """
    fwd_hits = 0
    rc_hits = 0
    for offset in range(1000, len(ref_seq) - 1000, step):
        kmer = ref_seq[offset : offset + kmer_len]
        if kmer in seq:
            fwd_hits += 1
        rc_kmer = str(Seq(kmer).reverse_complement())
        if rc_kmer in seq:
            rc_hits += 1
            
    is_rc = rc_hits > fwd_hits
    strand = "-" if is_rc else "+"
    return strand, is_rc


def find_circular_cut(seq: str, ref_seq: str, probe_len: int = 18, probe_span: int = 150, probe_step: int = 5) -> Tuple[int, int]:
    """
    Find circular breakpoint corresponding to coordinate 0 of reference sequence
    using sliding k-mer consensus voting on doubled sequence (seq + seq).
    """
    seq_len = len(seq)
    double_seq = seq + seq
    votes: Dict[int, int] = {}
    
    total_probes = 0
    for offset in range(0, probe_span, probe_step):
        kmer = ref_seq[offset : offset + probe_len]
        total_probes += 1
        idx = double_seq.find(kmer)
        if idx != -1:
            candidate_cut = (idx - offset) % seq_len
            votes[candidate_cut] = votes.get(candidate_cut, 0) + 1
            
    if not votes:
        # Fallback to shorter k-mers if high divergence at boundary
        for offset in range(0, 100, 5):
            kmer = ref_seq[offset : offset + 14]
            total_probes += 1
            idx = double_seq.find(kmer)
            if idx != -1:
                candidate_cut = (idx - offset) % seq_len
                votes[candidate_cut] = votes.get(candidate_cut, 0) + 1
                
    if not votes:
        raise ValueError("Failed to locate homologous circular start position against reference")
        
    best_cut = max(votes.items(), key=lambda x: x[1])[0]
    best_votes = votes[best_cut]
    return best_cut, best_votes


def standardize_sequence(
    record: SeqRecord,
    ref_seq: str
) -> Tuple[SeqRecord, Dict[str, any]]:
    """Standardize a single plastome record: reorient and relinearize."""
    raw_seq = str(record.seq).upper()
    orig_len = len(raw_seq)
    
    # 1. Orientation check
    strand, is_rc = detect_strand(raw_seq, ref_seq)
    oriented_seq = str(Seq(raw_seq).reverse_complement()) if is_rc else raw_seq
    
    # 2. Circular linearization
    cut_pos, votes = find_circular_cut(oriented_seq, ref_seq)
    std_seq = (oriented_seq + oriented_seq)[cut_pos : cut_pos + orig_len]
    
    # Verify homology of the first 25 bp
    start_homology = sum(1 for b1, b2 in zip(std_seq[:25], ref_seq[:25]) if b1 == b2)
    assert start_homology >= 20, f"Low start homology ({start_homology}/25) for {record.id}"
    
    std_record = SeqRecord(
        Seq(std_seq),
        id=record.id,
        name=record.name,
        description=f"{record.description} [standardized coordinate=NC_001320.1 strand={strand} cut_pos={cut_pos}]"
    )
    
    meta = {
        "Accession": record.id,
        "Original_Length_bp": orig_len,
        "Strand": strand,
        "Reversed": is_rc,
        "Circular_Cut_Offset": cut_pos,
        "Consensus_Votes": votes,
        "Standardized_Length_bp": len(std_seq),
        "Start_25bp": std_seq[:25],
        "Start_Homology_25bp": start_homology,
        "End_25bp": std_seq[-25:]
    }
    return std_record, meta


def main() -> None:
    print("=" * 80)
    print("Plastome Circular Standardization & Orientation Pipeline")
    print("=" * 80)
    
    ref_path = "data/raw/plastomes/sativa/reference_nipponbare_NC_001320.1.fasta"
    assert os.path.exists(ref_path), f"Reference genome {ref_path} not found"
    ref_rec = SeqIO.read(ref_path, "fasta")
    ref_seq = str(ref_rec.seq).upper()
    print(f"Loaded reference: {ref_rec.id} ({len(ref_seq):,} bp)")
    
    # Output directories
    out_dir_sativa = "data/interim/plastomes_standardized/sativa"
    out_dir_outgroups = "data/interim/plastomes_standardized/outgroups"
    os.makedirs(out_dir_sativa, exist_ok=True)
    os.makedirs(out_dir_outgroups, exist_ok=True)
    os.makedirs("results/tables", exist_ok=True)
    
    cohorts = [
        ("sativa", "data/raw/plastomes/sativa", out_dir_sativa),
        ("outgroups", "data/raw/plastomes/outgroups", out_dir_outgroups)
    ]
    
    all_standardized_records: List[SeqRecord] = []
    sativa_standardized_records: List[SeqRecord] = []
    summary_rows: List[Dict[str, any]] = []
    
    for cohort_name, in_dir, out_dir in cohorts:
        files = sorted(glob.glob(f"{in_dir}/*.fasta"))
        # Exclude multi-fastas and duplicate references
        files = [f for f in files if "all_" not in f and "reference_" not in f]
        print(f"\nProcessing {len(files)} {cohort_name} plastome files...")
        
        for fpath in files:
            rec = SeqIO.read(fpath, "fasta")
            std_rec, meta = standardize_sequence(rec, ref_seq)
            meta["Cohort"] = cohort_name
            meta["Source_File"] = os.path.basename(fpath)
            summary_rows.append(meta)
            
            # Write individual standardized fasta
            out_file = os.path.join(out_dir, f"{rec.id}.fasta")
            SeqIO.write(std_rec, out_file, "fasta")
            
            all_standardized_records.append(std_rec)
            if cohort_name == "sativa":
                sativa_standardized_records.append(std_rec)
                
            print(f"  [OK] {rec.id}: strand={meta['Strand']}, cut={meta['Circular_Cut_Offset']}, len={meta['Standardized_Length_bp']:,} bp")
            
    # Always include the canonical reference in sativa records if not already in files
    ref_std_rec = SeqRecord(
        Seq(ref_seq),
        id=ref_rec.id,
        name=ref_rec.name,
        description=f"{ref_rec.description} [canonical standard NC_001320.1]"
    )
    if not any(r.id == ref_rec.id for r in sativa_standardized_records):
        sativa_standardized_records.insert(0, ref_std_rec)
        all_standardized_records.insert(0, ref_std_rec)
        summary_rows.insert(0, {
            "Accession": ref_rec.id,
            "Cohort": "sativa_reference",
            "Source_File": "reference_nipponbare_NC_001320.1.fasta",
            "Original_Length_bp": len(ref_seq),
            "Strand": "+",
            "Reversed": False,
            "Circular_Cut_Offset": 0,
            "Consensus_Votes": 30,
            "Standardized_Length_bp": len(ref_seq),
            "Start_25bp": ref_seq[:25],
            "Start_Homology_25bp": 25,
            "End_25bp": ref_seq[-25:]
        })
        
    # Write combined multi-FASTAs
    all_sativa_path = "data/interim/plastomes_standardized/all_sativa_standardized.fasta"
    all_plastomes_path = "data/interim/plastomes_standardized/all_plastomes_standardized.fasta"
    SeqIO.write(sativa_standardized_records, all_sativa_path, "fasta")
    SeqIO.write(all_standardized_records, all_plastomes_path, "fasta")
    print(f"\nWrote combined Sativa multi-FASTA: {all_sativa_path} ({len(sativa_standardized_records)} records)")
    print(f"Wrote combined All-Plastomes multi-FASTA: {all_plastomes_path} ({len(all_standardized_records)} records)")
    
    # Generate cryptographic SHA-256 manifest
    manifest_path = "data/interim/plastomes_standardized/manifest.sha256"
    manifest_lines = []
    
    for root, _, fnames in os.walk("data/interim/plastomes_standardized"):
        for fname in sorted(fnames):
            if fname == "manifest.sha256":
                continue
            fpath = os.path.join(root, fname)
            rel_path = os.path.relpath(fpath, "data/interim/plastomes_standardized").replace("\\", "/")
            digest = compute_sha256(fpath)
            manifest_lines.append(f"{digest}  {rel_path}\n")
            
    with open(manifest_path, "w", encoding="utf-8") as f:
        f.writelines(manifest_lines)
    print(f"Generated cryptographic manifest: {manifest_path} ({len(manifest_lines)} entries)")
    
    # Save master summary table
    df_summary = pd.DataFrame(summary_rows)
    tsv_path = "results/tables/plastome_standardization_summary.tsv"
    csv_path = "results/tables/plastome_standardization_summary.csv"
    df_summary.to_csv(tsv_path, sep="\t", index=False)
    df_summary.to_csv(csv_path, index=False)
    print(f"Saved standardization summary table: {tsv_path}")
    
    # Verification report
    inverted_count = sum(1 for r in summary_rows if r["Reversed"])
    permuted_count = sum(1 for r in summary_rows if r["Circular_Cut_Offset"] != 0)
    print("\n" + "=" * 80)
    print("Standardization Metrics Summary:")
    print(f"  Total plastomes standardized: {len(summary_rows)}")
    print(f"  Reverse-complement strand inverted: {inverted_count}")
    print(f"  Circularly shifted relinearized: {permuted_count}")
    print(f"  All sequences verified 5' homologous to NC_001320.1: 100%")
    print("=" * 80)


if __name__ == "__main__":
    main()
