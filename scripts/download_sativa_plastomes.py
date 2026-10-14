"""
Download complete chloroplast genomes for Oryza sativa (curated cohort + Nipponbare reference)
from NCBI Entrez, validate sequence integrity, generate concatenated multi-FASTA,
and construct SHA-256 cryptographic manifest.
"""

import os
import sys
import hashlib
import time
from typing import List, Dict, Tuple
import pandas as pd
from Bio import Entrez, SeqIO
import io

def compute_sha256(filepath: str) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()

def download_sativa_plastomes(
    metadata_tsv: str = "results/tables/ncbi_oryza_sativa_chloroplast_genomes.tsv",
    outdir: str = "data/raw/plastomes/sativa",
    email: str = "popgen_researcher@example.com",
    reference_acc: str = "NC_001320.1"
):
    Entrez.email = email
    Entrez.tool = "GlobalRicePopgen"
    os.makedirs(outdir, exist_ok=True)
    
    # Read curated accessions and strip padded whitespace
    df = pd.read_csv(metadata_tsv, sep=r"\t", engine="python")
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    
    # Overwrite the CSV cleanly with proper quoting for compatibility
    clean_csv = "results/tables/ncbi_oryza_sativa_chloroplast_genomes.csv"
    df.to_csv(clean_csv, index=False)
    
    accessions: List[str] = df["Accession"].dropna().unique().tolist()
    
    # Ensure reference Nipponbare NC_001320.1 is included
    if reference_acc not in accessions:
        accessions.append(reference_acc)
        
    print(f"Total Oryza sativa accessions to download: {len(accessions)}")
    
    batch_size = 10
    downloaded_files: List[Tuple[str, str, int]] = []
    all_records = []
    
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        print(f"Fetching batch {i + 1} - {i + len(batch)} of {len(accessions)} ({', '.join(batch)})...")
        
        success = False
        for attempt in range(3):
            try:
                h = Entrez.efetch(db="nucleotide", id=",".join(batch), rettype="fasta", retmode="text")
                fasta_data = h.read()
                h.close()
                success = True
                break
            except Exception as e:
                print(f"  Attempt {attempt + 1} failed: {e}. Retrying in 2s...")
                time.sleep(2)
                
        if not success:
            raise RuntimeError(f"Failed to fetch batch {batch} after 3 attempts.")
            
        # Parse individual fasta records
        records = list(SeqIO.parse(io.StringIO(fasta_data), "fasta"))
        for rec in records:
            acc = rec.id
            clean_acc = acc.replace("/", "_")
            out_file = os.path.join(outdir, f"{clean_acc}.fasta")
            
            with open(out_file, "w", encoding="utf-8") as f:
                SeqIO.write(rec, f, "fasta")
                
            all_records.append(rec)
            downloaded_files.append((clean_acc, out_file, len(rec.seq)))
            
        time.sleep(0.4)
        
    # Write combined multi-FASTA
    combined_fasta = os.path.join(outdir, "all_sativa_plastomes.fasta")
    print(f"\nWriting combined multi-FASTA with {len(all_records)} sequences: {combined_fasta}")
    with open(combined_fasta, "w", encoding="utf-8") as f:
        SeqIO.write(all_records, f, "fasta")
        
    # If reference accession was downloaded, ensure dedicated alias
    for rec in all_records:
        if reference_acc in rec.id:
            ref_alias = os.path.join(outdir, f"reference_nipponbare_{reference_acc}.fasta")
            with open(ref_alias, "w", encoding="utf-8") as f:
                SeqIO.write(rec, f, "fasta")
            downloaded_files.append((f"reference_nipponbare_{reference_acc}", ref_alias, len(rec.seq)))
            break
            
    # Include combined fasta in manifest
    downloaded_files.append(("all_sativa_plastomes", combined_fasta, sum(len(r.seq) for r in all_records)))
    
    # Generate cryptographic SHA-256 manifest
    manifest_path = os.path.join(outdir, "manifest.sha256")
    print(f"Generating SHA-256 manifest: {manifest_path}")
    manifest_lines = []
    
    for _, fpath, seq_len in sorted(downloaded_files, key=lambda x: os.path.basename(x[1])):
        fname = os.path.basename(fpath)
        sha256_hash = compute_sha256(fpath)
        manifest_lines.append(f"{sha256_hash}  {fname}")
        
    with open(manifest_path, "w", encoding="utf-8") as f:
        f.write("\n".join(manifest_lines) + "\n")
        
    print(f"\nSuccessfully downloaded and validated {len(all_records)} Oryza sativa plastomes.")
    print(f"Manifest written with {len(manifest_lines)} file entries.")
    
    # Sequence length summary
    lengths = [len(r.seq) for r in all_records]
    print(f"Sequence length range: {min(lengths)} bp - {max(lengths)} bp (mean: {sum(lengths)/len(lengths):.1f} bp)")
    
    return all_records

if __name__ == "__main__":
    download_sativa_plastomes()
