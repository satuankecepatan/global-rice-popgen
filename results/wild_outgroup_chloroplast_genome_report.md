# Wild Outgroup Chloroplast Genome Curation & Retrieval Report

**Project:** Global Rice Population Genomics (`global-rice-popgen`)  
**Target Purpose:** Ancestral Allele Polarization & Maternal Phylogeny Outgroup Rooting  
**Total Outgroup Accessions:** 19  
**Generated Tables:**
- TSV: [`results/tables/ncbi_wild_outgroup_chloroplast_genomes.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/ncbi_wild_outgroup_chloroplast_genomes.tsv)
- CSV: [`results/tables/ncbi_wild_outgroup_chloroplast_genomes.csv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/ncbi_wild_outgroup_chloroplast_genomes.csv)
- Multi-FASTA: [`data/raw/plastomes/outgroups/all_wild_outgroup_plastomes.fasta`](file:///d:/Magister/Proyek/global-rice-popgen/data/raw/plastomes/outgroups/all_wild_outgroup_plastomes.fasta)
- Manifest: [`data/raw/plastomes/outgroups/manifest.sha256`](file:///d:/Magister/Proyek/global-rice-popgen/data/raw/plastomes/outgroups/manifest.sha256)

---

## 1. Executive Summary & Taxonomic Representation

Accurate estimation of the Site Frequency Spectrum (SFS), Fay & Wu's $H$, McDonald-Kreitman tests, and rooted maximum-likelihood genealogies requires unadulterated ancestral allele polarization. In Asian rice (*Oryza sativa*), using a single outgroup or an overly distant outgroup introduces alignment artifacts and homoplasy.

To establish a robust tripartite outgroup system, complete plastome assemblies were curated from three strategic evolutionary lineages:

1. ***Oryza rufipogon* ($N = 10$):**
   The perennial wild ancestor of Asian cultivated rice. Characterized by high outcrossing, broad geographic distribution across South/Southeast Asia and Southern China, and extensive genetic diversity. Crucially includes accession `KF562709.1` (DongXiang wild rice, the northernmost wild rice harboring cold-tolerance alleles pertinent to Hypothesis $H_2$), and NCBI RefSeq reference `NC_017835.1`.

2. ***Oryza nivara* ($N = 4$):**
   The annual wild ancestor of Asian cultivated rice, adapted to drought-prone, seasonally dry savannas and marshes. Represents the immediate wild sister group to *indica* / *aus* gene pools.

3. ***Oryza punctata* ($N = 5$):**
   A diploid wild rice possessing the BB nuclear genome (~400 Mb) and divergent plastome (~134.6–134.8 kb). Functions as the definitive ancestral outgroup to root the AA-genome clade (*O. sativa* + *O. rufipogon* + *O. nivara*) without confounding introgression. Includes NCBI RefSeq standard `NC_027676.1`.

---

## 2. Curated Accessions Summary

| Accession | Species | Evolutionary Role | Length (bp) | Sample / Voucher | Geographic Provenance | Submission Date |
| :--- | :--- | :--- | :---: | :--- | :--- | :---: |
| `PZ137539.1` | *Oryza rufipogon* | Perennial progenitor | 134,549 | ICROPS wild rice | Philippines | 16-MAR-2026 |
| `NC_017835.1` | *Oryza rufipogon* | Perennial progenitor | 134,544 | Wild accession | Viet Nam | 03-APR-2023 |
| `MW711322.1` | *Oryza rufipogon* | Perennial progenitor | 134,921 | Or136 | China | 21-APR-2021 |
| `MW711321.1` | *Oryza rufipogon* | Perennial progenitor | 134,963 | Or135 | China | 21-APR-2021 |
| `KF359902.1` | *Oryza rufipogon* | Perennial progenitor | 134,587 | Wild accession | China | 01-MAY-2017 |
| `KM103372.1` | *Oryza rufipogon* | Perennial progenitor | 134,510 | W1588 | Not reported in source | 03-AUG-2015 |
| `KF428978.1` | *Oryza rufipogon* | Perennial progenitor | 134,557 | Wild accession | Australia | 23-OCT-2014 |
| `KF562709.1` | *Oryza rufipogon* | Perennial progenitor | 134,537 | DongXiang (Cold-tolerant wild rice) | China | 06-NOV-2013 |
| `JN005833.1` | *Oryza rufipogon* | Perennial progenitor | 134,557 | AusTRCF 309313; AC01-1002323 | Australia | 24-FEB-2012 |
| `JN005832.1` | *Oryza rufipogon* | Perennial progenitor | 134,544 | Wild accession | Viet Nam | 24-FEB-2012 |
| `OL912836.1` | *Oryza nivara* | Annual progenitor | 134,555 | Wild accession | China | 17-OCT-2022 |
| `MZ450798.1` | *Oryza nivara* | Annual progenitor | 134,502 | Wild accession | China | 23-FEB-2022 |
| `LC666910.1` | *Oryza nivara* | Annual progenitor | 134,555 | CN_HLJ_W_S163 | China | 25-DEC-2021 |
| `LC642244.1` | *Oryza nivara* | Annual progenitor | 134,502 | CX19 | China | 21-JUL-2021 |
| `LC907692.1` | *Oryza punctata* | BB-genome ancestral outgroup | 134,844 | W1024 | Ghana | 14-JUL-2026 |
| `NC_027676.1` | *Oryza punctata* | BB-genome ancestral outgroup | 134,604 | GEN1233_2 | Not reported in source | 03-APR-2023 |
| `MT726932.1` | *Oryza punctata* | BB-genome ancestral outgroup | 134,829 | Wild accession | China | 15-SEP-2020 |
| `KF359908.1` | *Oryza punctata* | BB-genome ancestral outgroup | 134,582 | Wild accession | China | 01-MAY-2017 |
| `KM103375.1` | *Oryza punctata* | BB-genome ancestral outgroup | 134,604 | GEN1233_2 | Not reported in source | 03-AUG-2015 |

---

## 3. Data Integrity & Reproducibility

- **Sequence Length Range:** 134,502 bp – 134,963 bp (Mean: 134623.5 bp)
- **Cryptographic Checksum:** All sequences verified and recorded with SHA-256 hashes in `data/raw/plastomes/outgroups/manifest.sha256`.
- **Data Immutability:** Raw sequences are housed under `data/raw/` with zero in-place modification permitted pursuant to `AGENTS.md` Rule 2.

