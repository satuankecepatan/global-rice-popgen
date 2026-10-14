# Global Rice Population Genomics (`global-rice-popgen`)

> **Investigating Maternal Plastid Phylogeography, Cytonuclear Discordance, and Evolutionary Adaptation Across Global Rice (*Oryza sativa*) Germplasm**

---

## 1. Project Overview

Rice (*Oryza sativa* L.) is the foundational staple feeding over half the global human population. While extensive research has focused on nuclear genome diversity and domestication history, the maternal plastid genome (chloroplast DNA, ~134.5 kb) and its coevolutionary interplay with the nuclear genome remain critical yet under-characterized dimensions of global rice germplasm.

This project integrates computational population genomics, coalescent theory, and phylogenetic methods to dissect the maternal evolutionary lineages, geographic dispersal patterns, and environmental adaptations of Asian rice. Utilizing complete chloroplast genome assemblies alongside nuclear genomic datasets, we analyze wild ancestors, historical landraces, elite Green Revolution cultivars, and edge-range accessions across the globe.

---

## 2. Scientific Hypotheses

* **$H_1$ (Plastid Origin & Cytonuclear Discordance):**  
  Ancient independent chloroplast capture events and persistent asymmetric hybridization between *indica* and *japonica* gene pools have produced pervasive cytonuclear discordance. Distinct maternal plastome haplotypes (characterized by diagnostic ~50 bp length dimorphisms and specific indels) will reveal maternal donor lineages that diverge from predominant nuclear ancestry in key elite and hybrid cultivars.

* **$H_2$ (Adaptive Plastome Evolution at Geographic Margins):**  
  Plastid-encoded photosynthetic and metabolic core complexes (e.g., *rbcL*, *psbA*, *matK*, and the *ndh* gene family) harbor non-neutral nucleotide substitutions that correlate with adaptation to extreme environmental gradients—specifically cold tolerance at the northern range margin (e.g., high-latitude European/Russian cohorts) versus heat and drought resilience in tropical adapted landraces (e.g., *aus* / *N22*).

---

## 3. Proposed Outcomes

1. **Curated Global Rice Plastome Catalog:**  
   A high-integrity, contamination-free dataset of complete *Oryza sativa* and outgroup chloroplast genomes, cross-referenced with verified taxonomic classifications, geographic coordinates, and agronomic metadata.
2. **High-Resolution Maternal Phylogeographic Atlas:**  
   A robust genealogical tree and ancestral haplotype network mapping the maternal radiation of *japonica*, *indica*, *aus*, *aromatic*, and *spontanea* populations across global agroecological zones.
3. **Quantified Cytonuclear Discordance Matrix:**  
   Systematic identification and quantification of nuclear-cytoplasmic mismatch events in prominent breeding lines and hybrid parents (e.g., *PA64S*, *Milyang 23*).
4. **Selection & Nucleotide Diversity Signatures:**  
   Per-base callable site-adjusted nucleotide diversity ($\pi$, $\theta_W$), neutrality statistics (Tajima's $D$, Fay & Wu's $H$), and differentiation metrics (Hudson's $F_{ST}$, $d_{XY}$) across plastome functional domains and geographical populations.
5. **Open Reproducible Computational Pipeline:**  
   A standardized, container-ready workflow governed by strict data immutability, automated checksums, and deterministic random seeds.
6. **Peer-Reviewed Scientific Publication:**  
   A comprehensive manuscript prepared for submission to a leading genomics or plant biology journal (e.g., *Molecular Biology and Evolution*, *New Phytologist*, or *The Plant Journal*), accompanied by open-access data and code deposition.

---

## 4. Project Roadmap & Task Tracker

### Phase 1: Data Acquisition & Quality Curation
- [x] Establish isolated Python virtual environment (`.venv`) and pinned dependencies (`requirements.txt`)
- [x] Implement programmatic NCBI Entrez query pipeline for plastome retrieval ([`scripts/search_ncbi_chloroplast.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/search_ncbi_chloroplast.py))
- [x] Curate complete *Oryza sativa* chloroplast accessions ($N = 89$) and generate reference metadata tables ([`results/rice_chloroplast_genome_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/rice_chloroplast_genome_report.md))
- [ ] Download raw FASTA sequences and verify sequence integrity with SHA-256 manifests
- [ ] Retrieve wild outgroup plastomes (*Oryza rufipogon*, *Oryza nivara*, *Oryza punctata*) for ancestral allele polarization

### Phase 2: Sequence Alignment, Annotation & Variant Calling
- [ ] Perform whole-plastome multiple sequence alignment (e.g., MAFFT with manual refinement of inverted repeats)
- [ ] Define standardized coordinate system and coordinate-liftover against reference *Nipponbare* (`NC_001320.1`)
- [ ] Extract protein-coding genes (CDS), ribosomal RNAs (rRNA), transfer RNAs (tRNA), and intergenic spacers (IGS)
- [ ] Construct callable site bed mask excluding homopolymer tracts and unalignable hypervariable repeats
- [ ] Call and polarize SNPs and diagnostic indels into polarized VCF/BCF format

### Phase 3: Population Genetics & Phylogeographic Reconstruction
- [ ] Reconstruct maximum-likelihood and Bayesian plastid phylogenetic trees (IQ-TREE / BEAST)
- [ ] Build median-joining haplotype networks (TCS / PopART) resolving maternal sub-lineages
- [ ] Calculate unbiased diversity statistics ($\hat{\pi}$, $\hat{\theta}_W$, Tajima's $D$) using exact callable denominators ($L_{\text{callable}}$)
- [ ] Compute pairwise population differentiation (Hudson's $\hat{F}_{ST}$, $d_{XY}$) across varietal groups and regions
- [ ] Evaluate cytonuclear concordance by integrating paired nuclear genotype data (e.g., 3K Rice Genomes Project)

### Phase 4: Molecular Evolution & Environmental Adaptation
- [ ] Conduct gene-by-gene $d_N/d_S$ ($\omega$) selection tests (PAML / codeml and HyPhy RELAX/BUSTED)
- [ ] Test for cold-climate adaptive convergence in northern latitude accessions (Krasnodar cohort) vs. tropical lineages
- [ ] Screen for functional consequences of nonsynonymous mutations in photosystem and ATP synthase complexes
- [ ] Perform delete-$m_b$ weighted block jackknife for robust error estimation across genomic intervals

### Phase 5: Synthesis, Manuscript Preparation & Publication
- [ ] Compile publication-grade figures and vector graphics (SciencePlots styling, Nature/Cell formats)
- [ ] Draft manuscript text: Introduction, Methodology, Results, Discussion, and Data Availability
- [ ] Internal peer-review and critical validation against scientific governance invariants
- [ ] Archive datasets on Zenodo / NCBI BioProject and release pipeline on GitHub
- [ ] Deposit preprint on bioRxiv
- [ ] Journal submission, respond to reviewer revisions, and reach final publication

---

## 5. Repository Structure

```text
global-rice-popgen/
├── .agents/                      # Scientific governance rules & operational skills
│   ├── RULES/AGENTS.md           # Mandatory .venv, data immutability, and rigor standards
│   └── SKILLS/                   # Bioinformatics, Genomics, and Statistics specifications
├── .venv/                        # Project virtual environment (Python 3.11+)
├── data/
│   ├── raw/                      # Read-only raw sequence data and manifests (immutable)
│   ├── interim/                  # Intermediate alignments and callable bed masks
│   └── processed/                # Curated alignments, polarized VCFs, and tree files
├── requirements.txt              # Core scientific and genomics Python dependencies
├── results/
│   ├── figures/                  # Publication-ready vector plots and diagrams
│   ├── tables/                   # Curated metadata TSV/CSV tables
│   └── rice_chloroplast_genome_report.md  # NCBI curation report (89 curated plastomes)
├── scripts/
│   └── search_ncbi_chloroplast.py         # Automated NCBI Entrez search and curation tool
└── README.md                     # Project summary, hypotheses, outcomes, and task tracker
```

---

## 6. Quick Start & Execution

All workflows run strictly inside the project's dedicated virtual environment:

```powershell
# 1. Activate project virtual environment
.\.venv\Scripts\Activate.ps1

# 2. Verify environment and installed genomics libraries
python -c "import allel, Bio, pandas, numpy; print('Environment operational.')"

# 3. Query and update curated chloroplast genome metadata
python scripts\search_ncbi_chloroplast.py --output results\tables\ncbi_oryza_sativa_chloroplast_genomes.tsv
```

---

## 7. Citation & Scientific Governance

This repository adheres to reproducible research standards:
* **Estimator Rigor:** Ratio of averages for population metrics; zero minor allele frequency (MAF) pruning prior to SFS or demographic analysis.
* **Denominator Integrity:** Summary statistics strictly normalized by true callable invariant + variant sites ($L_{\text{callable}}$).
* **Provenance:** Every artifact is tracked with parameters, Git commit hashes, and random seeds.
