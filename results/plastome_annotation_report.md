# Rice Plastome Structural & Functional Annotation Report

> **Reference Genome:** *Oryza sativa* Japonica Group canonical standard RefSeq **NC_001320.1** (134,525 bp)

---

## 1. Executive Summary & Genomic Architecture

The *Oryza sativa* chloroplast genome is a 134,525 bp circular double-stranded DNA molecule exhibiting the canonical quadripartite organization characteristic of angiosperms. The genome comprises a **Large Single Copy (LSC)** region, a **Small Single Copy (SSC)** region, and two **Inverted Repeats (IRb and IRa)** that undergo flip-flop recombination.

| Genomic Compartment | 1-based Coordinates | Span (bp) | Proportion (%) | Coding Features |
| :--- | :--- | :--- | :--- | :--- |
| **Large Single Copy (LSC)** | 1 – 80,592 | 80,592 bp | 59.91% | Photosystems I/II, Cytochrome b6/f, Rubisco (*rbcL*), RNA polymerase |
| **Inverted Repeat B (IRb)** | 80,593 – 101,391 | 20,799 bp | 15.46% | Ribosomal RNAs (*rrn16, 23, 4.5, 5*), *rps12, rps7, ndhB* |
| **Small Single Copy (SSC)** | 101,392 – 113,726 | 12,335 bp | 9.17% | *ndh* family (*ndhF, ndhD, ndhE, ndhG, ndhI, ndhA, ndhH*), *rps15* |
| **Inverted Repeat A (IRa)** | 113,727 – 134,525 | 20,799 bp | 15.46% | Inverted duplicate copy of IRb |
| **Whole Plastome** | **1 – 134,525** | **134,525 bp** | **100.00%** | **106 CDS, 37 tRNA, 8 rRNA, 78 IGS** |

---

## 2. Functional Feature Breakdown

- **Protein-Coding Genes (CDS):** 106 CDS features spanning 155,871 bp (115.87% coding density).
- **Transfer RNAs (tRNA):** 37 tRNA genes spanning 9,954 bp.
- **Ribosomal RNAs (rRNA):** 8 rRNA genes (4 unique ribosomal operons duplicated in IRb and IRa) spanning 9,182 bp.
- **Intergenic Spacers (IGS):** 78 defined spacers spanning 26,191 bp (19.47% non-coding fraction).
- **Total Functional Density:** 130.09% of the genome encodes functional RNA or protein elements.

### Distribution Across Photosynthetic & Metabolic Complexes

| Functional Complex | Gene Count | Representative Genes | Biological Role |
| :--- | :--- | :--- | :--- |
| **Other / Conserved hypothetical** | 25 | hypothetical protein, ORF42, acetyl-CoA carboxylase beta subunit, photosystem I subunit IX | Core photosynthetic / metabolic apparatus |
| **Ribosomal protein (SSU)** | 15 | rps12, rps16, rps2, rps14, rps4, rps18 | Core photosynthetic / metabolic apparatus |
| **Photosystem II** | 15 | psbA, psbK, psbI, psbD, psbC, psbZ | Core photosynthetic / metabolic apparatus |
| **NADH dehydrogenase** | 12 | ndhJ, ndhK, ndhC, ndhB, ndhF, ndhD | Core photosynthetic / metabolic apparatus |
| **Ribosomal protein (LSU)** | 11 | rpl33, rpl20, rpl36, rpl14, rpl16, rpl22 | Core photosynthetic / metabolic apparatus |
| **Cytochrome b6/f complex** | 6 | petN, petA, petL, petG, petB, petD | Core photosynthetic / metabolic apparatus |
| **ATP synthase** | 6 | atpI, atpH, atpF, atpA, atpE, atpB | Core photosynthetic / metabolic apparatus |
| **RNA polymerase** | 4 | rpoB, rpoC1, rpoC2, rpoA | Core photosynthetic / metabolic apparatus |
| **Photosystem I** | 4 | psaB, psaA, psaI, psaC | Core photosynthetic / metabolic apparatus |
| **Conserved chloroplast ORF (ycf)** | 2 | ycf3, ycf4 | Core photosynthetic / metabolic apparatus |
| **Maturase K (splicing)** | 1 | matK | Core photosynthetic / metabolic apparatus |
| **Rubisco large subunit** | 1 | rbcL | Core photosynthetic / metabolic apparatus |
| **Chloroplast envelope membrane protein** | 1 | cemA | Core photosynthetic / metabolic apparatus |
| **ATP-dependent Clp protease** | 1 | clpP | Core photosynthetic / metabolic apparatus |
| **Translation initiation factor 1** | 1 | infA | Core photosynthetic / metabolic apparatus |
| **Cytochrome c biogenesis** | 1 | ccsA | Core photosynthetic / metabolic apparatus |

---

## 3. Key Adaptive Loci & Diagnostic Diagnostic Markers

1. **Rubisco Large Subunit (*rbcL*, LSC 53,959 – 55,392 bp):**  
   Primary catalytic engine of carbon fixation, key candidate for thermal adaptation and carbon assimilation efficiency.
2. **NADH Dehydrogenase Complex (*ndhA*–*K*, SSC & LSC):**  
   Critical for cyclic electron flow around Photosystem I, mediating drought, photoinhibition, and temperature stress resilience.
3. **Canonical ~50 bp Length Dimorphism (*rpl23-trnI* / ORF100 Spacer):**  
   Diagnostic indel situated between *rpl23* and *trnI* in the inverted repeat, distinguishing ancestral *indica* (deletion form) from *japonica* (insertion form).
4. **Photosystem II Core (*psbA*, LSC 335 – 1,396 bp):**  
   Encodes the D1 reaction center protein, turnover target under intense irradiance and cold-induced photoinhibition.

---

## 4. Deliverables Manifest

- **GFF3 Track:** [`data/interim/annotation/rice_chloroplast_NC_001320.1.gff3`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/annotation/rice_chloroplast_NC_001320.1.gff3)
- **BED Tracks:**
  - CDS: [`data/interim/annotation/cds.bed`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/annotation/cds.bed)
  - tRNA: [`data/interim/annotation/trna.bed`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/annotation/trna.bed)
  - rRNA: [`data/interim/annotation/rrna.bed`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/annotation/rrna.bed)
  - IGS: [`data/interim/annotation/igs.bed`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/annotation/igs.bed)
  - Quadripartite: [`data/interim/annotation/quadripartite.bed`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/annotation/quadripartite.bed)
- **TSV Tables:**
  - Master Gene Table: [`results/tables/plastome_gene_annotations.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_gene_annotations.tsv)
  - Quadripartite Architecture: [`results/tables/plastome_quadripartite_structure.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_quadripartite_structure.tsv)
  - Intergenic Spacers: [`results/tables/plastome_igs_annotations.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_igs_annotations.tsv)
