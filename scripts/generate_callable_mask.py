"""
Callable Sites Mask Generator
=============================
Implements Denominator Integrity (AGENTS.md Rule 3) for chloroplast population genomics.
Identifies and masks:
1. Homopolymer tracts (>= 8 bp identical nucleotide runs) prone to polymerase stutter / slippage.
2. Hypervariable simple sequence repeats (SSRs / microsatellites >= 12 bp).
3. Ambiguous characters (N) and terminal edge boundaries (10 bp).

Outputs 0-based BED intervals of callable and masked regions, computing exact callable
denominators (L_callable) across whole genome, quadripartite compartments (LSC, IRb, SSC, IRa),
single-copy partitioned domains, and individual CDS genes.
"""

import os
import re
from typing import List, Tuple, Dict
import pandas as pd
from Bio import SeqIO


def find_homopolymers(seq: str, min_len: int = 8) -> List[Tuple[int, int, str]]:
    """Locate homopolymer runs of identical nucleotides >= min_len."""
    pattern = re.compile(r"([A-Z])\1{" + str(min_len - 1) + r",}")
    masks = []
    for match in pattern.finditer(seq):
        start = match.start()
        end = match.end()
        base = match.group(1)
        length = end - start
        masks.append((start, end, f"homopolymer_{base}_{length}bp"))
    return masks


def find_ssrs(seq: str, min_total_len: int = 12) -> List[Tuple[int, int, str]]:
    """Find simple sequence repeats (di-, tri-, tetra-nucleotide tandem repeats)."""
    # 2-4 bp repeat units repeated at least min_total_len bp
    ssr_pattern = re.compile(r"(([A-Z]{2,4}))\1{2,}")
    masks = []
    for match in ssr_pattern.finditer(seq):
        start = match.start()
        end = match.end()
        if end - start >= min_total_len:
            unit = match.group(2)
            masks.append((start, end, f"SSR_unit_{unit}_{end-start}bp"))
    return masks


def merge_intervals(intervals: List[Tuple[int, int, str]]) -> List[Tuple[int, int, str]]:
    """Merge overlapping or adjacent intervals into non-redundant intervals."""
    if not intervals:
        return []
    sorted_int = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [sorted_int[0]]
    for curr in sorted_int[1:]:
        prev_s, prev_e, prev_desc = merged[-1]
        curr_s, curr_e, curr_desc = curr
        if curr_s <= prev_e:
            # Overlap
            new_e = max(prev_e, curr_e)
            new_desc = prev_desc if curr_desc in prev_desc else f"{prev_desc};{curr_desc}"
            merged[-1] = (prev_s, new_e, new_desc)
        else:
            merged.append(curr)
    return merged


def invert_intervals(merged_masks: List[Tuple[int, int, str]], total_len: int) -> List[Tuple[int, int]]:
    """Invert masked intervals to yield complementary callable intervals."""
    callable_intervals = []
    current_pos = 0
    for s, e, _ in merged_masks:
        if s > current_pos:
            callable_intervals.append((current_pos, s))
        current_pos = max(current_pos, e)
    if current_pos < total_len:
        callable_intervals.append((current_pos, total_len))
    return callable_intervals


def compute_interval_overlap(int1: Tuple[int, int], int2: Tuple[int, int]) -> int:
    """Compute number of overlapping bases between two 0-based half-open intervals."""
    s = max(int1[0], int2[0])
    e = min(int1[1], int2[1])
    return max(0, e - s)


def main() -> None:
    print("=" * 80)
    print("Callable Sites BED Mask Generation Pipeline (NC_001320.1)")
    print("=" * 80)
    
    ref_path = "data/raw/plastomes/sativa/reference_nipponbare_NC_001320.1.fasta"
    ref = SeqIO.read(ref_path, "fasta")
    seq = str(ref.seq).upper()
    genome_len = len(seq)
    print(f"Analyzing reference {ref.id} ({genome_len:,} bp)...")
    
    out_mask_dir = "data/interim/callable_mask"
    out_table_dir = "results/tables"
    os.makedirs(out_mask_dir, exist_ok=True)
    os.makedirs(out_table_dir, exist_ok=True)
    
    # 1. Identify raw mask intervals
    raw_masks: List[Tuple[int, int, str]] = []
    
    # A. Terminal boundaries (first 10 bp and last 10 bp)
    raw_masks.append((0, 10, "terminal_edge_5prime"))
    raw_masks.append((genome_len - 10, genome_len, "terminal_edge_3prime"))
    
    # B. Ambiguous bases (N)
    for m in re.finditer(r"[Nn]+", seq):
        raw_masks.append((m.start(), m.end(), "ambiguous_N"))
        
    # C. Homopolymer tracts (>= 8 bp)
    hp_masks = find_homopolymers(seq, min_len=8)
    print(f"Found {len(hp_masks)} homopolymer tracts (>= 8 bp)")
    raw_masks.extend(hp_masks)
    
    # D. Simple sequence repeats (SSRs >= 12 bp)
    ssr_masks = find_ssrs(seq, min_total_len=12)
    print(f"Found {len(ssr_masks)} simple sequence repeats (SSRs >= 12 bp)")
    raw_masks.extend(ssr_masks)
    
    # 2. Merge into non-overlapping masked intervals
    merged_masks = merge_intervals(raw_masks)
    total_masked_bp = sum(e - s for s, e, _ in merged_masks)
    print(f"Merged into {len(merged_masks)} non-overlapping masked intervals: {total_masked_bp:,} bp masked ({total_masked_bp/genome_len*100:.2f}%)")
    
    # 3. Derive callable intervals (complement of masked regions)
    callable_intervals = invert_intervals(merged_masks, genome_len)
    total_callable_bp = sum(e - s for s, e in callable_intervals)
    print(f"Calculated {len(callable_intervals)} callable intervals: {total_callable_bp:,} bp callable ({total_callable_bp/genome_len*100:.2f}%)")
    assert total_masked_bp + total_callable_bp == genome_len, "Masked + Callable does not sum to genome length!"
    
    # 4. Write BED files
    masked_bed_path = f"{out_mask_dir}/masked_regions.bed"
    with open(masked_bed_path, "w", encoding="utf-8") as f:
        for s, e, desc in merged_masks:
            f.write(f"{ref.id}\t{s}\t{e}\t{desc}\t0\t+\n")
            
    callable_bed_path = f"{out_mask_dir}/callable_sites.bed"
    with open(callable_bed_path, "w", encoding="utf-8") as f:
        for s, e in callable_intervals:
            f.write(f"{ref.id}\t{s}\t{e}\tcallable_site\t1000\t+\n")
    print(f"Wrote masked intervals BED: {masked_bed_path}")
    print(f"Wrote callable intervals BED: {callable_bed_path}")
    
    # 5. Compute callable denominators across genomic partitions
    compartments = [
        ("Whole_Plastome", (0, genome_len)),
        ("LSC", (0, 80592)),
        ("IRb", (80592, 101391)),
        ("SSC", (101391, 113726)),
        ("IRa", (113726, 134525)),
        ("Single_Copy_Partition", (0, 113726)) # LSC + IRb + SSC (excluding duplicated IRa)
    ]
    
    summary_rows: List[Dict[str, any]] = []
    for comp_name, comp_span in compartments:
        comp_len = comp_span[1] - comp_span[0]
        callable_bp = sum(compute_interval_overlap(comp_span, c_int) for c_int in callable_intervals)
        masked_bp = comp_len - callable_bp
        summary_rows.append({
            "Genomic_Partition": comp_name,
            "Start_0based": comp_span[0],
            "End_0based": comp_span[1],
            "Total_Span_bp": comp_len,
            "L_callable_bp": callable_bp,
            "Masked_Span_bp": masked_bp,
            "Callable_Percentage": round((callable_bp / comp_len) * 100, 3)
        })
        
    df_comp = pd.DataFrame(summary_rows)
    tsv_path = f"{out_table_dir}/callable_sites_summary.tsv"
    df_comp.to_csv(tsv_path, sep="\t", index=False)
    print(f"Saved callable sites summary table: {tsv_path}")
    
    # 6. Compute per-gene callable lengths for CDS features
    cds_bed_path = "data/interim/annotation/cds.bed"
    if os.path.exists(cds_bed_path):
        cds_rows = []
        with open(cds_bed_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    chrom, s0, e0, gene = parts[0], int(parts[1]), int(parts[2]), parts[3]
                    gene_span = (s0, e0)
                    gene_len = e0 - s0
                    c_bp = sum(compute_interval_overlap(gene_span, c_int) for c_int in callable_intervals)
                    cds_rows.append({
                        "Gene_Symbol": gene,
                        "Start_0based": s0,
                        "End_0based": e0,
                        "Gene_Length_bp": gene_len,
                        "L_callable_bp": c_bp,
                        "Callable_Percentage": round((c_bp / gene_len) * 100, 2)
                    })
        df_gene_callable = pd.DataFrame(cds_rows)
        gene_tsv = f"{out_table_dir}/gene_callable_lengths.tsv"
        df_gene_callable.to_csv(gene_tsv, sep="\t", index=False)
        print(f"Saved per-gene callable lengths table: {gene_tsv} ({len(df_gene_callable)} CDS features)")
        
    print("\n" + "=" * 80)
    print("Callable Sites Accounting Summary:")
    for row in summary_rows:
        print(f"  {row['Genomic_Partition']:<22}: Total={row['Total_Span_bp']:,} bp | L_callable={row['L_callable_bp']:,} bp ({row['Callable_Percentage']}%)")
    print("=" * 80)


if __name__ == "__main__":
    main()
