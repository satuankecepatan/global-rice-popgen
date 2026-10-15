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

* **$H_3$ (Genetic Architecture of Adaptive Traits & Cytonuclear Epistasis):**  
  Phenotypic divergence across global rice varieties (e.g., cold tolerance at range margins, drought resilience, and phenological timing) is driven by significant nuclear quantitative trait loci (QTLs) identifiable via mixed-model GWAS after rigorously controlling for hierarchical subpopulation structure and kinship, with evidence of cytonuclear epistasis between nuclear-encoded plastid-targeted proteins and maternal plastome variants.

---

## 3. Proposed Outcomes

1. **Curated Global Rice Plastome Catalog:**  
   A high-integrity, contamination-free dataset of complete *Oryza sativa* and outgroup chloroplast genomes, cross-referenced with verified taxonomic classifications, geographic coordinates, and agronomic metadata.
2. **Global Population Structure & Admixture Atlas:**  
   Comprehensive ancestral decomposition of global rice germplasm via Principal Component Analysis (PCA), cross-validated ADMIXTURE ancestry barplots ($K = 2\dots8$), and UMAP embeddings linking nuclear sub-taxa (*indica*, *japonica*, *aus*, *aromatic*, *spontanea*) to geographic agro-zones.
3. **High-Resolution Maternal Phylogeographic Atlas:**  
   A robust genealogical tree and ancestral haplotype network mapping the maternal radiation of rice populations, resolving maternal sub-lineages and regional divergence centers.
4. **Quantified Cytonuclear Discordance Matrix:**  
   Systematic identification and quantification of nuclear-cytoplasmic mismatch events in prominent breeding lines and hybrid parents (e.g., *PA64S*, *Milyang 23*).
5. **Genome-Wide Association Study (GWAS) & Locus Mapping:**  
   Identification of high-confidence causal loci for adaptive agronomic traits using mixed linear models (MLM) with kinship and population stratification correction, visualized via publication-grade Manhattan and QQ plots ($\lambda_{GC}$ calibrated).
6. **Selection & Nucleotide Diversity Signatures:**  
   Per-base callable site-adjusted nucleotide diversity ($\pi$, $\theta_W$), neutrality statistics (Tajima's $D$, Fay & Wu's $H$), and differentiation metrics (Hudson's $F_{ST}$, $d_{XY}$) across plastome functional domains and geographical populations.
7. **Open Reproducible Computational Pipeline:**  
   A standardized, container-ready workflow governed by strict data immutability, automated checksums, and deterministic random seeds.
8. **Peer-Reviewed Scientific Publication:**  
   A comprehensive manuscript prepared for submission to a leading genomics or plant biology journal (e.g., *Molecular Biology and Evolution*, *New Phytologist*, or *The Plant Journal*), accompanied by open-access data and code deposition.

---

## 4. Project Roadmap & Task Tracker

### Phase 1: Data Acquisition & Quality Curation
- [x] Establish isolated Python virtual environment (`.venv`) and pinned dependencies (`requirements.txt`)
- [x] Implement programmatic NCBI Entrez query pipeline for plastome retrieval ([`scripts/search_ncbi_chloroplast.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/search_ncbi_chloroplast.py))
- [x] Curate complete *Oryza sativa* chloroplast accessions ($N = 89$) and generate reference metadata tables ([`results/rice_chloroplast_genome_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/rice_chloroplast_genome_report.md))
- [x] Download raw FASTA sequences and verify sequence integrity with SHA-256 manifests ([`scripts/download_sativa_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/download_sativa_plastomes.py), [`data/raw/plastomes/sativa/manifest.sha256`](file:///d:/Magister/Proyek/global-rice-popgen/data/raw/plastomes/sativa/manifest.sha256))
- [x] Retrieve wild outgroup plastomes (*Oryza rufipogon*, *Oryza nivara*, *Oryza punctata*) for ancestral allele polarization ([`scripts/curate_and_download_outgroups.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/curate_and_download_outgroups.py), [`results/wild_outgroup_chloroplast_genome_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/wild_outgroup_chloroplast_genome_report.md))
- [x] Curate paired nuclear SNP datasets (3K Rice Genomes Project, $N = 3,024$) and passport/phenotypic metadata tables ([`scripts/curate_nuclear_3k_metadata.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/curate_nuclear_3k_metadata.py), [`results/rice_3k_nuclear_curation_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/rice_3k_nuclear_curation_report.md))

### Phase 2: Sequence Alignment, Annotation & Variant Calling
- [x] Standardize plastome circular orientation and linear breakpoints against *Nipponbare* (`NC_001320.1`) ([`scripts/standardize_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/standardize_plastomes.py), [`results/tables/plastome_standardization_summary.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_standardization_summary.tsv))
- [x] Extract protein-coding genes (CDS), ribosomal RNAs (rRNA), transfer RNAs (tRNA), and intergenic spacers (IGS) ([`scripts/extract_plastome_annotations.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/extract_plastome_annotations.py), [`results/plastome_annotation_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/plastome_annotation_report.md))
- [x] Construct callable site bed mask excluding homopolymer tracts and simple sequence repeats ([`scripts/generate_callable_mask.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/generate_callable_mask.py), [`data/interim/callable_mask/callable_sites.bed`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/callable_mask/callable_sites.bed))
- [x] Perform whole-plastome multiple sequence alignment and single-copy partitioning ([`scripts/align_and_liftover_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/align_and_liftover_plastomes.py), [`results/alignments/manifest.sha256`](file:///d:/Magister/Proyek/global-rice-popgen/results/alignments/manifest.sha256))
- [x] Define standardized coordinate system and coordinate-liftover table against reference *Nipponbare* ([`results/tables/coordinate_liftover_index.tsv.gz`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/coordinate_liftover_index.tsv.gz))
- [x] Call and polarize SNPs and diagnostic indels into polarized VCF/BCF format ([`scripts/call_and_polarize_variants.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/call_and_polarize_variants.py), [`data/processed/variants/rice_plastome_variants_polarized.vcf.gz`](file:///d:/Magister/Proyek/global-rice-popgen/data/processed/variants/rice_plastome_variants_polarized.vcf.gz), [`results/plastome_alignment_variant_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/plastome_alignment_variant_report.md))

### Phase 3: Population Structure, Phylogeography & Cytonuclear Dynamics
- [ ] Perform LD-pruning of biallelic SNPs ($r^2 < 0.2$) for unbiased population structure inference
- [ ] Compute genomic relationship matrix (GRM / Kinship $K$) and run Principal Component Analysis (PCA)
- [ ] Model ancestral admixture proportions (ADMIXTURE / fastSTRUCTURE, $K = 2\dots8$) and generate stacked barplots
- [ ] Generate 2D/3D UMAP manifold embeddings colored by varietal groups and geographic provenance
- [ ] Reconstruct maximum-likelihood and Bayesian plastid phylogenetic trees (IQ-TREE / BEAST)
- [ ] Build median-joining haplotype networks (TCS / PopART) resolving maternal sub-lineages
- [ ] Calculate unbiased diversity statistics ($\hat{\pi}$, $\hat{\theta}_W$, Tajima's $D$) using exact callable denominators ($L_{\text{callable}}$)
- [ ] Compute pairwise population differentiation (Hudson's $\hat{F}_{ST}$, $d_{XY}$) across varietal groups and regions
- [ ] Quantify cytonuclear discordance by projecting nuclear Q-matrices against plastome clades

### Phase 4: Genome-Wide Association Studies (GWAS) & Adaptive Evolution
- [ ] Harmonize phenotypic datasets (e.g., cold tolerance indices, heading date, drought resilience metrics)
- [ ] Execute mixed linear model (MLM / FaST-LMM) association mapping incorporating fixed structure covariates ($Q$/PCs) and random kinship ($K$)
- [ ] Generate publication-grade Manhattan and Quantile-Quantile (QQ) plots with `scienceplots` and `adjustText`
- [ ] Evaluate genomic inflation factor ($\lambda_{GC}$) and establish significance thresholds (Bonferroni & Benjamini-Hochberg FDR)
- [ ] Fine-map lead candidate loci and assess functional overlap with plastid-nuclear interaction networks
- [ ] Conduct gene-by-gene $d_N/d_S$ ($\omega$) selection tests (PAML / codeml and HyPhy RELAX/BUSTED)
- [ ] Perform delete-$m_b$ weighted block jackknife for robust error estimation across genomic intervals

### Phase 5: Synthesis, Manuscript Preparation & Publication
- [ ] Compile publication-grade composite figures (PCA + ADMIXTURE panels, Manhattan/QQ plots, and phylogenetic trees)
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
│   │   ├── plastomes/
│   │   │   ├── sativa/           # 90 complete O. sativa plastomes + manifest.sha256
│   │   │   └── outgroups/        # 19 wild outgroup plastomes + manifest.sha256
│   │   └── nuclear_3k/           # 3K RGP download manifest + manifest.sha256
│   ├── interim/                  # Intermediate alignments and callable bed masks
│   └── processed/                # Curated alignments, polarized VCFs, and tree files
├── requirements.txt              # Core scientific and genomics Python dependencies
├── results/
│   ├── figures/                  # Publication-ready vector plots and diagrams
│   ├── tables/                   # Curated metadata TSV/CSV tables (sativa, outgroups, 3K passport, phenotypes, crossref)
│   ├── rice_chloroplast_genome_report.md    # NCBI curation report (89 curated plastomes)
│   ├── wild_outgroup_chloroplast_genome_report.md # Wild outgroups report (19 outgroups)
│   └── rice_3k_nuclear_curation_report.md   # 3K RGP nuclear genomics & phenotype curation report
├── scripts/
│   ├── search_ncbi_chloroplast.py           # Automated NCBI Entrez search and curation tool
│   ├── download_sativa_plastomes.py         # Batch download & SHA-256 verification of sativa plastomes
│   ├── curate_and_download_outgroups.py     # Wild outgroup retrieval, curation & manifest generator
│   └── curate_nuclear_3k_metadata.py       # 3K RGP passport, Q-matrix, phenotype harmonization & crossref
├── tests/
│   └── test_phase1.py                       # Automated verification suite for Phase 1 data integrity
└── README.md                     # Project summary, hypotheses, outcomes, and task tracker
```

---

## 6. Quick Start & Execution

All workflows run strictly inside the project's dedicated virtual environment:

```powershell
# 1. Activate project virtual environment
.\.venv\Scripts\Activate.ps1

# 2. Verify environment and run Phase 1 automated tests
pytest tests\test_phase1.py -v

# 3. Execute data acquisition pipelines
python scripts\download_sativa_plastomes.py
python scripts\curate_and_download_outgroups.py
python scripts\curate_nuclear_3k_metadata.py
```

---

## 7. Citation & Scientific Governance

This repository adheres to reproducible research standards:
* **Estimator Rigor:** Ratio of averages for population metrics; zero minor allele frequency (MAF) pruning prior to SFS or demographic analysis.
* **Denominator Integrity:** Summary statistics strictly normalized by true callable invariant + variant sites ($L_{\text{callable}}$).
* **Provenance:** Every artifact is tracked with parameters, Git commit hashes, and random seeds.
