"""
Search NCBI Nucleotide database for strictly Oryza sativa complete chloroplast genomes,
extract metadata (cultivar name, isolate name, nation/geographic origin, length, subspecies),
and generate reproducible summary tables and reports.
"""

import os
import sys
import argparse
from typing import List, Dict, Any
import pandas as pd
from Bio import Entrez, SeqIO
import io
import time

def parse_args():
    parser = argparse.ArgumentParser(description="Query NCBI for Oryza sativa chloroplast genomes.")
    parser.add_argument("--email", default="popgen_researcher@example.com", help="NCBI Entrez email")
    parser.add_argument("--outdir", default="results/tables", help="Output directory for generated tables")
    return parser.parse_args()

def fetch_and_curate(email: str, outdir: str):
    Entrez.email = email
    Entrez.tool = "GlobalRicePopgen"
    
    os.makedirs(outdir, exist_ok=True)
    
    # Strict query
    strict_query = '"Oryza sativa"[Organism] AND "chloroplast complete genome"'
    print(f"Executing strict Entrez query: {strict_query}")
    h = Entrez.esearch(db="nucleotide", term=strict_query, retmax=500)
    res = Entrez.read(h)
    h.close()
    
    id_list = res["IdList"]
    print(f"Found {len(id_list)} records matching strict organism and title phrase.")
    
    # Title variation search to capture complete chloroplast genomes that use minor phrasing differences
    broad_query = '("Oryza sativa"[Organism]) AND (chloroplast[Title] AND "complete genome"[Title])'
    h2 = Entrez.esearch(db="nucleotide", term=broad_query, retmax=500)
    res2 = Entrez.read(h2)
    h2.close()
    
    all_uids = list(dict.fromkeys(id_list + res2["IdList"]))
    print(f"Total UIDs retrieved across strict & title variations: {len(all_uids)}")
    
    # Batch fetch GenBank records
    batch_size = 30
    records: List[Dict[str, Any]] = []
    
    for i in range(0, len(all_uids), batch_size):
        batch = all_uids[i:i + batch_size]
        print(f"Fetching batch {i + 1} - {i + len(batch)} of {len(all_uids)} from NCBI...")
        fetch_h = Entrez.efetch(db="nucleotide", id=",".join(batch), rettype="gb", retmode="text")
        raw_text = fetch_h.read()
        fetch_h.close()
        
        for rec in SeqIO.parse(io.StringIO(raw_text), "genbank"):
            acc = rec.id
            desc = rec.description
            length = len(rec.seq)
            org = rec.annotations.get("organism", "")
            date = rec.annotations.get("date", "")
            
            source_feats = [f for f in rec.features if f.type == "source"]
            q = source_feats[0].qualifiers if source_feats else {}
            
            cultivar_raw = q.get("cultivar", [None])[0]
            isolate_raw = q.get("isolate", [None])[0]
            country_raw = q.get("country", [None])[0] or q.get("geo_loc_name", [None])[0]
            voucher = q.get("specimen_voucher", [None])[0]
            strain = q.get("strain", [None])[0]
            note = q.get("note", [None])[0]
            
            refs = rec.annotations.get("references", [])
            ref_titles = " | ".join([r.title for r in refs if r.title])
            ref_journals = " | ".join([r.journal for r in refs if r.journal])
            
            # Curate cultivar
            cultivar_curated = cultivar_raw if cultivar_raw else ""
            if not cultivar_curated:
                if "PA64S" in str(isolate_raw) or "PA64S" in desc:
                    cultivar_curated = "PA64S"
                elif "Jingeng818" in str(isolate_raw):
                    cultivar_curated = "Jingeng 818"
                elif "Nipponbare" in desc:
                    cultivar_curated = "Nipponbare"
                elif "93-11" in desc:
                    cultivar_curated = "93-11"
                elif "HSAGSDYD1802" in str(isolate_raw):
                    cultivar_curated = "HSAGSDYD1802 (Yellow mutant)"
                elif "Hassawi" in desc or (note and "Hassawi" in note):
                    cultivar_curated = "Hassawi hybrid (IR1112 x Hassawi)" if "hybrid" in desc else "Hassawi"
                elif "Taichung 65" in desc:
                    cultivar_curated = "Taichung 65"
                elif "spontanea" in org or "nivara" in desc.lower():
                    cultivar_curated = "Wild form (O. sativa f. spontanea / O. nivara)"
                elif note:
                    cultivar_curated = f"Line {note}"
                else:
                    cultivar_curated = "Not specified"
            
            # Curate isolate
            isolate_curated = isolate_raw if isolate_raw else ""
            if not isolate_curated:
                if voucher:
                    isolate_curated = f"Voucher: {voucher}"
                elif strain:
                    isolate_curated = f"Strain: {strain}"
                elif "PA64S" in desc:
                    isolate_curated = "PA64S"
                elif "93-11" in desc:
                    isolate_curated = "93-11"
                elif "W0038" in desc or "W0038" in ref_journals:
                    isolate_curated = "W0038"
                else:
                    isolate_curated = "Not specified"
                    
            # Curate nation origin
            origin_nation = ""
            if cultivar_curated == "Nipponbare":
                origin_nation = "Japan"
            elif cultivar_curated in ["93-11", "PA64S", "Zhenshan 97", "Minghui 63", "Shuhui498", "SR86", "Jingeng 818"]:
                origin_nation = "China"
            elif cultivar_curated in ["TN1", "Taichung 65"]:
                origin_nation = "Taiwan" if cultivar_curated == "TN1" else "Taiwan / Japan"
            elif cultivar_curated == "IR8":
                origin_nation = "Philippines (IRRI)"
            elif cultivar_curated in ["N22", "RP Bio-226"]:
                origin_nation = "India"
            elif cultivar_curated == "basmati":
                origin_nation = "Pakistan / India"
            elif "Hassawi" in cultivar_curated:
                origin_nation = "Saudi Arabia"
            elif cultivar_curated == "Pin Gaew 56":
                origin_nation = "Thailand"
            elif cultivar_curated in ["Milyang 23", "Yukara"]:
                origin_nation = "South Korea" if cultivar_curated == "Milyang 23" else "Japan"
            elif country_raw:
                origin_nation = country_raw.split(":")[0].strip()
            elif "Russia" in ref_journals:
                origin_nation = "Russia"
            elif "China" in ref_journals or "Qingdao" in ref_journals or "Kunming" in ref_journals:
                origin_nation = "China"
            elif "Korea" in ref_journals or "Seoul National" in ref_journals:
                origin_nation = "South Korea"
            else:
                origin_nation = "Unspecified"
                
            records.append({
                "Accession": acc,
                "Length_bp": length,
                "Taxon_Organism": org,
                "Cultivar": cultivar_curated,
                "Isolate": isolate_curated,
                "GenBank_Country_Field": country_raw if country_raw else "Not reported in source",
                "Nation_Origin": origin_nation,
                "Definition": desc,
                "Submission_Date": date
            })
        time.sleep(0.4)
        
    df = pd.DataFrame(records).drop_duplicates(subset=["Accession"])
    
    out_csv = os.path.join(outdir, "ncbi_oryza_sativa_chloroplast_genomes.csv")
    out_tsv = os.path.join(outdir, "ncbi_oryza_sativa_chloroplast_genomes.tsv")
    
    df.to_csv(out_csv, index=False)
    df.to_csv(out_tsv, sep="\t", index=False)
    
    print(f"Successfully processed {len(df)} records.")
    print(f"Outputs written to:\n - {out_csv}\n - {out_tsv}")
    return df

if __name__ == "__main__":
    args = parse_args()
    fetch_and_curate(args.email, args.outdir)
