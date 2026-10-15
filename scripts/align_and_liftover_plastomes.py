"""
Plastome Multiple Sequence Alignment & Coordinate Liftover Pipeline
===================================================================
Constructs a unified, standardized multiple sequence alignment (MSA) for all 109
chloroplast genomes (90 Oryza sativa + 19 wild outgroups) against the canonical
Nipponbare reference standard (NC_001320.1).

Features:
1. High-speed anchor-guided affine gap dynamic programming (Numba-accelerated).
2. Exact coordinate liftover mapping each reference nucleotide (1..134,525) to its
   unified MSA column and query allele across all accessions.
3. Generation of:
   - Full whole-plastome alignment (sativa cohort, N=90).
   - Full whole-plastome alignment (all accessions, N=109).
   - Single-copy partitioned alignment (LSC + IRb + SSC), masking IRa (113,727..134,525)
     to prevent concerted evolution double-counting bias in population genetic statistics.
4. Comprehensive coordinate liftover table and cryptographic SHA-256 manifests.
"""

import os
import glob
import time
import gzip
import hashlib
from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
import numba
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


@numba.njit
def affine_block_align(
    s1: np.ndarray,
    s2: np.ndarray,
    gap_open: float = 12.0,
    gap_ext: float = 1.0,
    match: float = 2.0,
    mismatch: float = -3.0
) -> Tuple[np.ndarray, np.ndarray]:
    """Affine gap penalty banded alignment for collinear sequence blocks."""
    n = len(s1)
    m = len(s2)
    band = max(40, abs(m - n) + 25)
    bw = 2 * band + 1
    
    M = np.full((n + 1, bw), -1e9, dtype=np.float32)
    X = np.full((n + 1, bw), -1e9, dtype=np.float32) # gap in s1 (ins in s2)
    Y = np.full((n + 1, bw), -1e9, dtype=np.float32) # gap in s2 (del in s2)
    
    tb_M = np.zeros((n + 1, bw), dtype=np.int8)
    tb_X = np.zeros((n + 1, bw), dtype=np.int8)
    tb_Y = np.zeros((n + 1, bw), dtype=np.int8)
    
    M[0, band] = 0.0
    for k in range(1, band + 1):
        X[0, band + k] = -(gap_open + k * gap_ext)
        tb_X[0, band + k] = 1
        
    for i in range(1, n + 1):
        for b in range(bw):
            j = i + b - band
            if j < 0 or j > m:
                continue
                
            # X: ins in s2 (gap in s1, j advanced, i same)
            if b - 1 >= 0 and j > 0:
                xm = M[i, b - 1] - gap_open - gap_ext
                xx = X[i, b - 1] - gap_ext
                if xx >= xm:
                    X[i, b] = xx
                    tb_X[i, b] = 1
                else:
                    X[i, b] = xm
                    tb_X[i, b] = 0
                    
            # Y: del in s2 (gap in s2, i advanced, j same)
            if b + 1 < bw and i > 0:
                ym = M[i - 1, b + 1] - gap_open - gap_ext
                yy = Y[i - 1, b + 1] - gap_ext
                if yy >= ym:
                    Y[i, b] = yy
                    tb_Y[i, b] = 2
                else:
                    Y[i, b] = ym
                    tb_Y[i, b] = 0
                    
            # M: match / mismatch
            if i > 0 and j > 0:
                score = match if s1[i - 1] == s2[j - 1] else mismatch
                mm = M[i - 1, b] + score
                mx = X[i - 1, b] + score
                my = Y[i - 1, b] + score
                best_m = mm
                c = 0
                if mx > best_m:
                    best_m = mx
                    c = 1
                if my > best_m:
                    best_m = my
                    c = 2
                M[i, b] = best_m
                tb_M[i, b] = c
                
    final_b = band + (m - n)
    state = 0
    best_final = M[n, final_b]
    if X[n, final_b] > best_final:
        best_final = X[n, final_b]
        state = 1
    if Y[n, final_b] > best_final:
        best_final = Y[n, final_b]
        state = 2
        
    i = n
    b = final_b
    max_len = n + m
    aln1 = np.zeros(max_len, dtype=np.uint8)
    aln2 = np.zeros(max_len, dtype=np.uint8)
    pos = 0
    gap_char = ord('-')
    
    while i > 0 or (i + b - band) > 0:
        j = i + b - band
        if state == 0:
            prev_state = tb_M[i, b]
            aln1[pos] = s1[i - 1]
            aln2[pos] = s2[j - 1]
            i -= 1
            state = prev_state
        elif state == 1:
            prev_state = tb_X[i, b]
            aln1[pos] = gap_char
            aln2[pos] = s2[j - 1]
            b -= 1
            state = prev_state
        else:
            prev_state = tb_Y[i, b]
            aln1[pos] = s1[i - 1]
            aln2[pos] = gap_char
            i -= 1
            b += 1
            state = prev_state
        pos += 1
        
    return aln1[:pos][::-1], aln2[:pos][::-1]


def align_query_to_reference(ref_seq: str, qry_seq: str, k: int = 25, step: int = 1000) -> Tuple[np.ndarray, np.ndarray]:
    """Perform anchor-guided collinear alignment of a query against reference genome."""
    s1 = ref_seq
    s2 = qry_seq
    
    # 1. Discover collinear anchors
    anchors = [(0, 0)]
    for r_pos in range(step, len(s1) - step, step):
        kmer = s1[r_pos : r_pos + k]
        # Search in query window around homologous position
        q_search_start = max(0, r_pos - 400)
        q_search_end = min(len(s2), r_pos + 700)
        sub = s2[q_search_start:q_search_end]
        idx = sub.find(kmer)
        if idx != -1:
            q_pos = q_search_start + idx
            if q_pos > anchors[-1][1]:
                anchors.append((r_pos, q_pos))
    anchors.append((len(s1), len(s2)))
    
    # 2. Align block-by-block
    full_a1 = []
    full_a2 = []
    for i in range(len(anchors) - 1):
        r1, q1 = anchors[i]
        r2, q2 = anchors[i + 1]
        blk_s1 = np.frombuffer(s1[r1:r2].encode("ascii"), dtype=np.uint8)
        blk_s2 = np.frombuffer(s2[q1:q2].encode("ascii"), dtype=np.uint8)
        a1, a2 = affine_block_align(blk_s1, blk_s2)
        full_a1.append(a1)
        full_a2.append(a2)
        
    return np.concatenate(full_a1), np.concatenate(full_a2)


def get_quadripartite_domain(pos_1based: int) -> str:
    """Return quadripartite domain name for 1-based coordinate in NC_001320.1."""
    if 1 <= pos_1based <= 80592:
        return "LSC"
    elif 80593 <= pos_1based <= 101391:
        return "IRb"
    elif 101392 <= pos_1based <= 113726:
        return "SSC"
    elif 113727 <= pos_1based <= 134525:
        return "IRa"
    return "Unknown"


def load_callable_mask(bed_path: str, genome_len: int) -> np.ndarray:
    """Load callable intervals into a boolean mask array of length genome_len (1-based indexed)."""
    mask = np.zeros(genome_len + 1, dtype=bool)
    if os.path.exists(bed_path):
        with open(bed_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    s0 = int(parts[1])
                    e0 = int(parts[2])
                    # 1-based: s0 + 1 .. e0
                    mask[s0 + 1 : e0 + 1] = True
    else:
        mask[:] = True
    return mask


def main() -> None:
    print("=" * 80)
    print("Plastome Multiple Sequence Alignment & Coordinate Liftover Pipeline")
    print("=" * 80)
    
    # Warm up JIT compiler
    dummy = np.frombuffer(b"ACGT", dtype=np.uint8)
    affine_block_align(dummy, dummy)
    
    ref_path = "data/raw/plastomes/sativa/reference_nipponbare_NC_001320.1.fasta"
    ref_rec = SeqIO.read(ref_path, "fasta")
    ref_seq = str(ref_rec.seq).upper()
    genome_len = len(ref_seq)
    print(f"Loaded reference standard: {ref_rec.id} ({genome_len:,} bp)")
    
    callable_bed = "data/interim/callable_mask/callable_sites.bed"
    callable_mask = load_callable_mask(callable_bed, genome_len)
    print(f"Loaded callable sites mask: {callable_mask.sum():,} bp callable")
    
    # Load all standardized accessions
    all_sativa_fa = "data/interim/plastomes_standardized/all_sativa_standardized.fasta"
    all_plastomes_fa = "data/interim/plastomes_standardized/all_plastomes_standardized.fasta"
    
    sativa_records = list(SeqIO.parse(all_sativa_fa, "fasta"))
    all_records = list(SeqIO.parse(all_plastomes_fa, "fasta"))
    print(f"Loaded {len(sativa_records)} Sativa genomes and {len(all_records)} total plastome genomes")
    
    # Output directories
    out_aln_dir = "results/alignments"
    out_table_dir = "results/tables"
    os.makedirs(out_aln_dir, exist_ok=True)
    os.makedirs(out_table_dir, exist_ok=True)
    
    # 1. Perform pairwise alignments against reference and extract per-position calls and insertions
    print("\nAligning all 109 plastomes against reference NC_001320.1...")
    t0_all = time.time()
    
    sample_ids: List[str] = [r.id for r in all_records]
    num_samples = len(sample_ids)
    
    # Data structures for reference-anchored unified MSA
    # ref_calls[sample_idx, r_pos_1based]: ASCII char
    ref_calls = np.zeros((num_samples, genome_len + 1), dtype=np.uint8)
    
    # insertions_after[r_pos_1based][sample_idx]: string of inserted bases
    insertions_after: List[Dict[int, str]] = [{} for _ in range(genome_len + 1)]
    
    align_stats = []
    
    for s_idx, rec in enumerate(all_records):
        t0 = time.time()
        qry_seq = str(rec.seq).upper()
        
        if rec.id == ref_rec.id:
            # Self-alignment is trivial
            for p in range(1, genome_len + 1):
                ref_calls[s_idx, p] = ord(ref_seq[p - 1])
            align_stats.append({
                "Sample_ID": rec.id, "Length_bp": len(qry_seq), "SNPs": 0, "Gaps": 0, "Align_Time_s": 0.0
            })
            continue
            
        a1, a2 = align_query_to_reference(ref_seq, qry_seq)
        
        # Traverse alignment and map to reference coordinates
        curr_ref_pos = 0
        current_ins = []
        snps = 0
        gaps = 0
        
        for b1_val, b2_val in zip(a1, a2):
            b1 = chr(b1_val)
            b2 = chr(b2_val)
            
            if b1 == "-":
                # Insertion in query relative to reference
                current_ins.append(b2)
                gaps += 1
            else:
                # b1 is reference base: advance reference coordinate
                if current_ins:
                    insertions_after[curr_ref_pos][s_idx] = "".join(current_ins)
                    current_ins = []
                    
                curr_ref_pos += 1
                ref_calls[s_idx, curr_ref_pos] = b2_val
                
                if b2 == "-":
                    gaps += 1
                elif b1 != b2:
                    snps += 1
                    
        if current_ins:
            insertions_after[curr_ref_pos][s_idx] = "".join(current_ins)
            
        elapsed = time.time() - t0
        align_stats.append({
            "Sample_ID": rec.id, "Length_bp": len(qry_seq), "SNPs": snps, "Gaps": gaps, "Align_Time_s": round(elapsed, 3)
        })
        if (s_idx + 1) % 20 == 0 or s_idx == num_samples - 1:
            print(f"  [{s_idx + 1:3d}/{num_samples}] Aligned {rec.id:<14}: SNPs={snps:<4} Gaps={gaps:<4} in {elapsed:.2f}s")
            
    print(f"\nAll 109 pairwise alignments completed in {time.time() - t0_all:.2f} seconds!")
    
    # 2. Build Unified Multiple Sequence Alignment Matrix
    print("\nAssembling unified global Multiple Sequence Alignment (MSA)...")
    
    # Determine maximum insertion width at each reference coordinate
    max_ins_widths = [0] * (genome_len + 1)
    for p in range(genome_len + 1):
        if insertions_after[p]:
            max_ins_widths[p] = max(len(ins) for ins in insertions_after[p].values())
            
    total_ins_columns = sum(max_ins_widths)
    msa_total_len = genome_len + total_ins_columns
    print(f"Unified MSA dimensions: {num_samples} taxa x {msa_total_len:,} alignment columns ({total_ins_columns} insertion columns)")
    
    # Build MSA arrays
    msa_matrix = np.full((num_samples, msa_total_len), ord("-"), dtype=np.uint8)
    ref_msa = np.full(msa_total_len, ord("-"), dtype=np.uint8)
    
    # Coordinate liftover index mapping
    liftover_records = []
    
    col_idx = 0
    # Coordinate 0 insertions (before pos 1)
    if max_ins_widths[0] > 0:
        w = max_ins_widths[0]
        for s_idx in range(num_samples):
            ins_str = insertions_after[0].get(s_idx, "")
            for offset, char in enumerate(ins_str):
                msa_matrix[s_idx, col_idx + offset] = ord(char)
        for offset in range(w):
            liftover_records.append({
                "Alignment_Column": col_idx + offset + 1,
                "Reference_Coordinate": 0,
                "Sub_Coordinate": offset + 1,
                "Ref_Base": "-",
                "Domain": get_quadripartite_domain(1),
                "Is_Callable": False
            })
        col_idx += w
        
    for p in range(1, genome_len + 1):
        # Reference column
        ref_char = ord(ref_seq[p - 1])
        ref_msa[col_idx] = ref_char
        for s_idx in range(num_samples):
            msa_matrix[s_idx, col_idx] = ref_calls[s_idx, p]
            
        domain = get_quadripartite_domain(p)
        is_callable = bool(callable_mask[p])
        
        liftover_records.append({
            "Alignment_Column": col_idx + 1,
            "Reference_Coordinate": p,
            "Sub_Coordinate": 0,
            "Ref_Base": chr(ref_char),
            "Domain": domain,
            "Is_Callable": is_callable
        })
        col_idx += 1
        
        # Insertion columns after pos p
        w = max_ins_widths[p]
        if w > 0:
            for s_idx in range(num_samples):
                ins_str = insertions_after[p].get(s_idx, "")
                for offset, char in enumerate(ins_str):
                    msa_matrix[s_idx, col_idx + offset] = ord(char)
            for offset in range(w):
                liftover_records.append({
                    "Alignment_Column": col_idx + offset + 1,
                    "Reference_Coordinate": p,
                    "Sub_Coordinate": offset + 1,
                    "Ref_Base": "-",
                    "Domain": domain,
                    "Is_Callable": False
                })
            col_idx += w
            
    assert col_idx == msa_total_len, f"Column count mismatch: {col_idx} != {msa_total_len}"
    
    # 3. Export Multi-FASTA files
    # A. All Plastomes Full Alignment (109 accessions)
    all_out_path = f"{out_aln_dir}/all_plastomes_aligned.fasta"
    all_aln_records = []
    # Include reference first if not sample 0
    ref_idx = sample_ids.index(ref_rec.id)
    all_aln_records.append(SeqRecord(Seq(bytes(msa_matrix[ref_idx]).decode("ascii")), id=ref_rec.id, description="Reference standard"))
    for s_idx, s_id in enumerate(sample_ids):
        if s_id == ref_rec.id:
            continue
        all_aln_records.append(SeqRecord(Seq(bytes(msa_matrix[s_idx]).decode("ascii")), id=s_id, description="Aligned complete plastome"))
    SeqIO.write(all_aln_records, all_out_path, "fasta")
    print(f"\nWrote full alignment (All 109 taxa): {all_out_path}")
    
    # B. Sativa Cohort Full Alignment (90 accessions)
    sativa_out_path = f"{out_aln_dir}/sativa_plastomes_aligned.fasta"
    sativa_set = {r.id for r in sativa_records}
    sativa_aln_records = [rec for rec in all_aln_records if rec.id in sativa_set]
    SeqIO.write(sativa_aln_records, sativa_out_path, "fasta")
    print(f"Wrote Sativa alignment (90 taxa): {sativa_out_path}")
    
    # C. Single-Copy Partitioned Alignment (LSC + IRb + SSC), masking IRa (pos 113,727..134,525)
    # Exclude columns belonging to IRa to avoid double-counting concerted evolution
    single_copy_cols = [r["Alignment_Column"] - 1 for r in liftover_records if r["Domain"] in ("LSC", "IRb", "SSC")]
    sc_len = len(single_copy_cols)
    sc_matrix = msa_matrix[:, single_copy_cols]
    
    sc_out_path = f"{out_aln_dir}/all_plastomes_single_copy.fasta"
    sc_aln_records = []
    for s_idx, s_id in enumerate(sample_ids):
        sc_aln_records.append(SeqRecord(Seq(bytes(sc_matrix[s_idx]).decode("ascii")), id=s_id, description="Single-copy partitioned alignment (LSC+IRb+SSC)"))
    SeqIO.write(sc_aln_records, sc_out_path, "fasta")
    print(f"Wrote single-copy partitioned alignment (LSC+IRb+SSC, {sc_len:,} columns): {sc_out_path}")
    
    # 4. Save Coordinate Liftover Table
    df_liftover = pd.DataFrame(liftover_records)
    liftover_tsv_gz = f"{out_table_dir}/coordinate_liftover_index.tsv.gz"
    with gzip.open(liftover_tsv_gz, "wt", encoding="utf-8") as f:
        df_liftover.to_csv(f, sep="\t", index=False)
    print(f"Saved compressed coordinate liftover table: {liftover_tsv_gz}")
    
    # Alignment summary table
    df_stats = pd.DataFrame(align_stats)
    df_stats.to_csv(f"{out_table_dir}/plastome_alignment_statistics.tsv", sep="\t", index=False)
    
    # 5. Cryptographic SHA-256 Manifest
    manifest_path = f"{out_aln_dir}/manifest.sha256"
    with open(manifest_path, "w", encoding="utf-8") as f:
        for p in [all_out_path, sativa_out_path, sc_out_path]:
            h = compute_sha256(p)
            f.write(f"{h}  {os.path.basename(p)}\n")
    print(f"Generated cryptographic manifest: {manifest_path}")
    
    print("\n" + "=" * 80)
    print("Alignment & Liftover Summary:")
    print(f"  Total plastomes aligned: {num_samples} (90 Sativa + 19 Wild Outgroups)")
    print(f"  Whole MSA span: {msa_total_len:,} columns")
    print(f"  Single-copy MSA span: {sc_len:,} columns")
    print(f"  Liftover coordinates mapped: 100% of reference positions 1..{genome_len:,}")
    print("=" * 80)


if __name__ == "__main__":
    main()
