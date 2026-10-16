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
* [x] Implement programmatic NCBI Entrez query pipeline for plastome retrieval ([`scripts/search_ncbi_chloroplast.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/search_ncbi_chloroplast.py))
* [x] Curate complete *Oryza sativa* chloroplast accessions ($N = 89$) and generate reference metadata tables ([`results/rice_chloroplast_genome_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/rice_chloroplast_genome_report.md))
* [x] Download raw FASTA sequences and verify sequence integrity with SHA-256 manifests ([`scripts/download_sativa_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/download_sativa_plastomes.py), [`data/raw/plastomes/sativa/manifest.sha256`](file:///d:/Magister/Proyek/global-rice-popgen/data/raw/plastomes/sativa/manifest.sha256))
* [x] Retrieve wild outgroup plastomes (*Oryza rufipogon*, *Oryza nivara*, *Oryza punctata*) for ancestral allele polarization ([`scripts/curate_and_download_outgroups.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/curate_and_download_outgroups.py), [`results/wild_outgroup_chloroplast_genome_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/wild_outgroup_chloroplast_genome_report.md))
* [x] Curate paired nuclear SNP datasets (3K Rice Genomes Project, $N = 3,024$) and passport/phenotypic metadata tables ([`scripts/curate_nuclear_3k_metadata.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/curate_nuclear_3k_metadata.py), [`results/rice_3k_nuclear_curation_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/rice_3k_nuclear_curation_report.md))

### Phase 2: Sequence Alignment, Annotation & Variant Calling

- [x] Standardize plastome circular orientation and linear breakpoints against *Nipponbare* (`NC_001320.1`) ([`scripts/standardize_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/standardize_plastomes.py), [`results/tables/plastome_standardization_summary.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_standardization_summary.tsv))
* [x] Extract protein-coding genes (CDS), ribosomal RNAs (rRNA), transfer RNAs (tRNA), and intergenic spacers (IGS) ([`scripts/extract_plastome_annotations.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/extract_plastome_annotations.py), [`results/plastome_annotation_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/plastome_annotation_report.md))
* [x] Construct callable site bed mask excluding homopolymer tracts and simple sequence repeats ([`scripts/generate_callable_mask.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/generate_callable_mask.py), [`data/interim/callable_mask/callable_sites.bed`](file:///d:/Magister/Proyek/global-rice-popgen/data/interim/callable_mask/callable_sites.bed))
* [x] Perform whole-plastome multiple sequence alignment and single-copy partitioning ([`scripts/align_and_liftover_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/align_and_liftover_plastomes.py), [`results/alignments/manifest.sha256`](file:///d:/Magister/Proyek/global-rice-popgen/results/alignments/manifest.sha256))
* [x] Define standardized coordinate system and coordinate-liftover table against reference *Nipponbare* ([`results/tables/coordinate_liftover_index.tsv.gz`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/coordinate_liftover_index.tsv.gz))
* [x] Call and polarize SNPs and diagnostic indels into polarized VCF/BCF format ([`scripts/call_and_polarize_variants.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/call_and_polarize_variants.py), [`data/processed/variants/rice_plastome_variants_polarized.vcf.gz`](file:///d:/Magister/Proyek/global-rice-popgen/data/processed/variants/rice_plastome_variants_polarized.vcf.gz), [`results/plastome_alignment_variant_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/plastome_alignment_variant_report.md))

### Phase 3: Population Structure, Phylogeography & Cytonuclear Dynamics

- [x] Calculate unbiased diversity statistics ($\hat{\pi}$, $\hat{\theta}_W$, Tajima's $D$, Fay & Wu's $H$) using exact callable denominators ($L_{\text{callable}}$) ([`scripts/compute_diversity_and_divergence.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/compute_diversity_and_divergence.py), [`results/tables/plastome_diversity_neutrality_summary.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_diversity_neutrality_summary.tsv), [`results/plastome_diversity_divergence_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/plastome_diversity_divergence_report.md))
* [x] Compute pairwise population differentiation (Hudson's $\hat{F}_{ST}$, $d_{XY}$, $d_A$) across varietal groups and regions with delete-$m_b$ weighted block jackknife ([`results/tables/plastome_pairwise_divergence_fst.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_pairwise_divergence_fst.tsv), [`results/figures/fig1_plastome_diversity_neutrality.pdf`](file:///d:/Magister/Proyek/global-rice-popgen/results/figures/fig1_plastome_diversity_neutrality.pdf), [`results/figures/fig3_pairwise_differentiation_heatmaps.pdf`](file:///d:/Magister/Proyek/global-rice-popgen/results/figures/fig3_pairwise_differentiation_heatmaps.pdf))
* [ ] Perform LD-pruning of biallelic SNPs ($r^2 < 0.2$) for unbiased population structure inference
* [ ] Compute genomic relationship matrix (GRM / Kinship $K$) and run Principal Component Analysis (PCA)
* [ ] Model ancestral admixture proportions (ADMIXTURE / fastSTRUCTURE, $K = 2\dots8$) and generate stacked barplots
* [ ] Generate 2D/3D UMAP manifold embeddings colored by varietal groups and geographic provenance
* [x] Reconstruct maternal plastid phylogenetic trees rooted on basal *O. punctata* with 100 bootstrap replicates ([`scripts/reconstruct_phylogeny_and_network.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/reconstruct_phylogeny_and_network.py), [`results/tables/rice_maternal_plastome_tree.nwk`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/rice_maternal_plastome_tree.nwk), [`results/figures/fig5_maternal_plastid_phylogeny.pdf`](file:///d:/Magister/Proyek/global-rice-popgen/results/figures/fig5_maternal_plastid_phylogeny.pdf))
* [x] Build median-joining haplotype networks resolving 28 maternal lineages and mutational steps ([`results/tables/plastome_haplotypes.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_haplotypes.tsv), [`results/figures/fig6_maternal_haplotype_network.pdf`](file:///d:/Magister/Proyek/global-rice-popgen/results/figures/fig6_maternal_haplotype_network.pdf), [`results/plastome_phylogeny_haplotype_network_report.md`](file:///d:/Magister/Proyek/global-rice-popgen/results/plastome_phylogeny_haplotype_network_report.md))
* [ ] Quantify cytonuclear discordance by projecting nuclear Q-matrices against plastome clades

### Phase 4: Genome-Wide Association Studies (GWAS) & Adaptive Evolution

- [ ] Harmonize phenotypic datasets (e.g., cold tolerance indices, heading date, drought resilience metrics)
* [ ] Execute mixed linear model (MLM / FaST-LMM) association mapping incorporating fixed structure covariates ($Q$/PCs) and random kinship ($K$)
* [ ] Generate publication-grade Manhattan and Quantile-Quantile (QQ) plots with `scienceplots` and `adjustText`
* [ ] Evaluate genomic inflation factor ($\lambda_{GC}$) and establish significance thresholds (Bonferroni & Benjamini-Hochberg FDR)
* [ ] Fine-map lead candidate loci and assess functional overlap with plastid-nuclear interaction networks
* [ ] Conduct gene-by-gene $d_N/d_S$ ($\omega$) selection tests (PAML / codeml and HyPhy RELAX/BUSTED)
* [ ] Perform delete-$m_b$ weighted block jackknife for robust error estimation across genomic intervals

### Phase 5: Synthesis, Manuscript Preparation & Publication

- [ ] Compile publication-grade composite figures (PCA + ADMIXTURE panels, Manhattan/QQ plots, and phylogenetic trees)
* [ ] Draft manuscript text: Introduction, Methodology, Results, Discussion, and Data Availability
* [ ] Internal peer-review and critical validation against scientific governance invariants
* [ ] Archive datasets on Zenodo / NCBI BioProject and release pipeline on GitHub
* [ ] Deposit preprint on bioRxiv
* [ ] Journal submission, respond to reviewer revisions, and reach final publication

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
│   ├── interim/                  # Intermediate standardized plastomes & callable masks
│   │   ├── plastomes_standardized/ # Standardized 5'-linear forward-strand fastas
│   │   └── callable_mask/        # Callable site BED masks & exclusion tracks
│   └── processed/                # Curated alignments, polarized VCFs, and tree sequences
│       ├── variants/             # Polarized VCF 4.2 datasets + diagnostic indels
│       └── phylogeny/            # Maternal plastid trees (.nwk, .nex, .xml)
├── requirements.txt              # Core scientific and genomics Python dependencies
├── results/
│   ├── alignments/               # Whole-plastome and single-copy partitioned MSAs
│   ├── figures/                  # Publication-ready vector plots (PDF/PNG)
│   ├── tables/                   # Curated metadata TSVs, diversity summaries, & liftover index
│   ├── rice_chloroplast_genome_report.md    # NCBI curation report (89 curated plastomes)
│   ├── wild_outgroup_chloroplast_genome_report.md # Wild outgroups report (19 outgroups)
│   ├── rice_3k_nuclear_curation_report.md   # 3K RGP nuclear genomics & phenotype report
│   ├── plastome_annotation_report.md        # Feature extraction & coordinate report
│   ├── plastome_alignment_variant_report.md # MSA liftover & variant calling report
│   ├── plastome_diversity_divergence_report.md # Population diversity & divergence report
│   └── plastome_phylogeny_haplotype_network_report.md # Phylogeny & haplotype network report
├── scripts/
│   ├── search_ncbi_chloroplast.py           # Automated NCBI Entrez search and curation tool
│   ├── download_sativa_plastomes.py         # Batch download & SHA-256 verification of sativa plastomes
│   ├── curate_and_download_outgroups.py     # Wild outgroup retrieval, curation & manifest generator
│   ├── curate_nuclear_3k_metadata.py       # 3K RGP passport, Q-matrix, phenotype harmonization & crossref
│   ├── standardize_plastomes.py             # Circular linearization & 5'-strand orientation
│   ├── extract_plastome_annotations.py      # CDS, rRNA, tRNA, and IGS feature extraction
│   ├── generate_callable_mask.py            # Homopolymer, SSR, and callable bed mask generator
│   ├── align_and_liftover_plastomes.py      # Numba affine-gap alignment & coordinate liftover
│   ├── call_and_polarize_variants.py        # Outgroup polarization, VCF 4.2 caller & indel marker
│   ├── compute_diversity_and_divergence.py  # pi, theta_W, Tajima's D, Hudson's Fst & block jackknife
│   └── reconstruct_phylogeny_and_network.py # K80/JC69 NJ trees, 100 bootstraps & median-joining network
├── tests/
│   ├── test_phase1.py                       # Automated verification for Phase 1 data integrity
│   ├── test_phase2.py                       # Automated tests for Phase 2 alignment, masking & VCF
│   └── test_phase3.py                       # Analytical invariant tests for diversity, Fst & trees
├── tools/
│   └── mafft-win/                           # MAFFT multiple sequence aligner distribution
└── README.md                                # Project summary, hypotheses, outcomes, methods & citations
```

---

## 6. Quick Start & Execution

All workflows run strictly inside the project's dedicated virtual environment:

```powershell
# 1. Activate project virtual environment
.\.venv\Scripts\Activate.ps1

# 2. Run complete automated test suite across all completed phases
pytest tests\ -v

# 3. Phase 1: Execute data acquisition pipelines
python scripts\download_sativa_plastomes.py
python scripts\curate_and_download_outgroups.py
python scripts\curate_nuclear_3k_metadata.py

# 4. Phase 2: Alignment, masking & polarized variant calling
python scripts\standardize_plastomes.py
python scripts\extract_plastome_annotations.py
python scripts\generate_callable_mask.py
python scripts\align_and_liftover_plastomes.py
python scripts\call_and_polarize_variants.py

# 5. Phase 3: Population diversity, divergence, phylogeny & haplotype networks
python scripts\compute_diversity_and_divergence.py
python scripts\reconstruct_phylogeny_and_network.py
python scripts\compute_population_structure.py
```

---

## 7. CLI Programs & Computational Tooling (To Cite for Reproducibility)

To ensure full computational reproducibility in accordance with FAIR (Findable, Accessible, Interoperable, and Reusable) principles, all external command-line programs, bioinformatics suites, and numerical libraries utilized across this project are cataloged below with their operational roles, implementation entry points, and primary literature citations.

| CLI Program / Tool | Pipeline Role & Operational Purpose | Script / Module Link | Primary Literature Citation |
| :--- | :--- | :--- | :--- |
| **MAFFT** (v7.x) | High-speed multiple sequence alignment (MSA) of full chloroplast genomes and collinear single-copy partitions using progressive fast Fourier transform algorithms. | [`tools/mafft-win`](file:///d:/Magister/Proyek/global-rice-popgen/tools/mafft-win) | Katoh & Standley (2013) *Mol Biol Evol* [[DOI](https://doi.org/10.1093/molbev/mst010)] |
| **SAMtools** (v1.x) | Reference sequence indexing (`faidx`), BAM/CRAM coordinate sorting, duplicate marking (`markdup`), and binary alignment processing. | [`.agents/SKILLS/Bioinformatics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Bioinformatics.md) | Danecek et al. (2021) *GigaScience* [[DOI](https://doi.org/10.1093/gigascience/giab008)] |
| **BCFtools** (v1.x) | Polarized VCF/BCF stream manipulation, CSI/TBI indexing, multi-allelic variant decomposition, and tabix-based interval querying. | [`scripts/call_and_polarize_variants.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/call_and_polarize_variants.py) | Danecek et al. (2021) *GigaScience* [[DOI](https://doi.org/10.1093/gigascience/giab008)] |
| **BWA-MEM2 / BWA** | High-performance Burrow-Wheeler Transform (BWT) / FM-index sequence read alignment against the *Nipponbare* reference genome. | [`.agents/SKILLS/Bioinformatics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Bioinformatics.md) | Vasimuddin et al. (2019) *IEEE IPDPS* [[DOI](https://doi.org/10.1109/IPDPS.2019.00041)] |
| **GATK** (v4.x) | Genome Analysis Toolkit for sequence dictionary generation (`CreateSequenceDictionary`), joint genotyping, and genomic interval definitions. | [`.agents/SKILLS/Bioinformatics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Bioinformatics.md) | Poplin et al. (2018) *Nat Biotechnol* [[DOI](https://doi.org/10.1038/nbt.4235)] |
| **NCBI Entrez Direct** | Automated batch query, programmatic taxonomy resolution, and XML/FASTA retrieval of chloroplast accessions. | [`scripts/search_ncbi_chloroplast.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/search_ncbi_chloroplast.py) | Kans (2020) *NCBI Handbook* [[NLM](https://www.ncbi.nlm.nih.gov/books/NBK179288/)] |
| **mosdepth** | High-performance per-base and windowed sequencing depth calculation and callable site coverage profiling. | [`.agents/SKILLS/Bioinformatics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Bioinformatics.md) | Pedersen & Quinlan (2018) *Bioinformatics* [[DOI](https://doi.org/10.1093/bioinformatics/btx699)] |
| **PLINK** (v1.9 / v2.0) | High-performance genotype matrix manipulation, pairwise linkage disequilibrium pruning ($r^2 < 0.2$), and Principal Component Analysis. | [`.agents/SKILLS/Genomics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Genomics.md) | Purcell et al. (2007) *Am J Hum Genet* [[DOI](https://doi.org/10.1086/519795)] |
| **ADMIXTURE** | Maximum likelihood estimation of individual genetic ancestry proportions across $K = 2\dots8$ ancestral populations. | [`.agents/SKILLS/Genomics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Genomics.md) | Alexander, Novembre, & Lange (2009) *Genome Res* [[DOI](https://doi.org/10.1101/gr.094052.109)] |
| **fastSTRUCTURE** | Variational Bayesian inference of population structure and admixture proportions under the Pritchard model. | [`.agents/SKILLS/Genomics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Genomics.md) | Raj, Stephens, & Pritchard (2014) *Genetics* [[DOI](https://doi.org/10.1534/genetics.114.164350)] |
| **PAML (codeml)** | Maximum likelihood codon substitution modeling and gene-by-gene nonsynonymous-to-synonymous rate ratio ($d_N/d_S$, $\omega$) selection scans. | [`.agents/SKILLS/Genomics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Genomics.md) | Yang (2007) *Mol Biol Evol* [[DOI](https://doi.org/10.1093/molbev/msm088)] |
| **HyPhy** | Hypothesis testing in phylogenetics using codon models (RELAX for selection relaxation/intensification and BUSTED for branch-site positive selection). | [`.agents/SKILLS/Genomics.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/SKILLS/Genomics.md) | Kosakovsky Pond et al. (2020) *Mol Biol Evol* [[DOI](https://doi.org/10.1093/molbev/msaa159)] |
| **Snakemake** | Deterministic workflow engine managing directed acyclic graph (DAG) job dependency scheduling and computational reproducibility. | [`.agents/RULES/AGENTS.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/RULES/AGENTS.md) | Mölder et al. (2021) *F1000Research* [[DOI](https://doi.org/10.12688/f1000research.29032.2)] |
| **pytest** | Automated test runner executing analytical invariant verification, coalescent recovery tests, and mathematical boundary audits. | [`tests/`](file:///d:/Magister/Proyek/global-rice-popgen/tests) | Krekel et al. (2004) *pytest.org* |
| **BioPython** | Core sequence manipulation, GenBank parsing, Phylo tree structures, and distance matrix construction. | [`requirements.txt`](file:///d:/Magister/Proyek/global-rice-popgen/requirements.txt) | Cock et al. (2009) *Bioinformatics* [[DOI](https://doi.org/10.1093/bioinformatics/btp163)] |
| **scikit-allel** | Population genetics algorithms: folded/unfolded SFS, windowed diversity, pairwise $F_{ST}$, Patterson's $D$, and PCA. | [`requirements.txt`](file:///d:/Magister/Proyek/global-rice-popgen/requirements.txt) | Miles et al. (2020) [[Zenodo](https://doi.org/10.5281/zenodo.3759363)] |
| **NetworkX** | Graph-theoretical algorithms for minimum spanning trees, graph distance matrices, and median-joining haplotype networks. | [`scripts/reconstruct_phylogeny_and_network.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/reconstruct_phylogeny_and_network.py) | Hagberg, Swart, & Chult (2008) *SciPy Proceedings* |

---

## 8. Algorithms, Mathematical Estimators & Theoretical Models (To Cite for Reproducibility)

Every statistic and transformation in this pipeline is governed by strict mathematical estimators rather than ad-hoc heuristics. Downstream publications should cite the following mathematical formulations:

### 8.1. Sequence Standardization & Collinear Dynamic Programming

* **Multi-kmer Strand & Linearization Voting:**  
  Plastome assemblies in public repositories exhibit arbitrary circular origin breakpoints and strand inversion. We implement multi-kmer voting ($k=31\text{ bp}$, step $=5,000\text{ bp}$) against the canonical *Nipponbare* RefSeq (`NC_001320.1`) to identify forward/reverse-complement orientation and perform exact circular permutation so all sequences linearize at coordinate 1.  
  *Script:* [`scripts/standardize_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/standardize_plastomes.py)
* **Banded Affine-Gap Dynamic Programming:**  
  Collinear inter-accession blocks are aligned via Needleman-Wunsch with Gotoh's affine gap penalty formulation ($g_{\text{open}} = 12.0$, $g_{\text{ext}} = 1.0$, $s_{\text{match}} = 2.0$, $s_{\text{mismatch}} = -3.0$) accelerated via Numba JIT compilation:
  $$M_{i,j} = \max\begin{cases} M_{i-1,j-1} + S(s1_i, s2_j) \\ X_{i-1,j-1} + S(s1_i, s2_j) \\ Y_{i-1,j-1} + S(s1_i, s2_j) \end{cases}$$
  $$X_{i,j} = \max\begin{cases} M_{i-1,j} - g_{\text{open}} \\ X_{i-1,j} - g_{\text{ext}} \end{cases}, \qquad Y_{i,j} = \max\begin{cases} M_{i,j-1} - g_{\text{open}} \\ Y_{i,j-1} - g_{\text{ext}} \end{cases}$$
  *Citation:* Gotoh (1982) *J Mol Biol* [[DOI](https://doi.org/10.1016/0022-2836(82)90398-9)]  
  *Script:* [`scripts/align_and_liftover_plastomes.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/align_and_liftover_plastomes.py)

### 8.2. Denominator Integrity & Variant Polarization

* **Callable Sites Masking Algorithm (AGENTS.md Rule 3):**  
  To eliminate sequencing slippage and false-positive calls, we construct exact 0-based BED masks filtering:
  1. Homopolymers: identical nucleotide runs $\ge 8\text{ bp}$.
  2. Simple sequence repeats (SSRs): di-, tri-, and tetra-nucleotide tandem repeats $\ge 12\text{ bp}$.
  3. Terminal edge boundaries: $10\text{ bp}$ flanking sequence margins.  
  All per-base diversity metrics are normalized by the exact count of callable bases ($L_{\text{callable}} = 133,526\text{ bp}$ whole plastome; $112,814\text{ bp}$ single-copy partition).  
  *Script:* [`scripts/generate_callable_mask.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/generate_callable_mask.py)
* **Single-Copy Partitioning (IRa Masking):**  
  Chloroplast genomes feature two large Inverted Repeats (IRb and IRa, ~20.8 kb each) maintained by active high-fidelity gene conversion. Analyzing both copies artificially inflates sample sizes and skews neutral coalescent expectations. We mask IRa ($113,727\dots134,525\text{ bp}$), retaining the Single-Copy Partition ($\text{LSC} + \text{IRb} + \text{SSC}$, $118,545\text{ bp}$).
* **Ancestral Allele Polarization:**  
  Variants across the 90 *O. sativa* accessions are polarized using parsimony against basal wild outgroups (*Oryza punctata*, *O. rufipogon*, *O. nivara*). Only sites where the ancestral state is unambiguously resolved are populated with the derived allele frequency (DAF) and derived counts ($c_s$) for the unfolded Site Frequency Spectrum (SFS).  
  *Script:* [`scripts/call_and_polarize_variants.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/call_and_polarize_variants.py)
* **Diagnostic Indel Length Dimorphism:**  
  Identification of the canonical ~50 bp indel length dimorphism (insertion in *Japonica*, deletion in *Indica*) serving as a diagnostic plastome marker.

### 8.3. Population Genetic Diversity & Neutrality Estimators

All summary statistics are calculated as linear combinations of the unfolded Site Frequency Spectrum $\boldsymbol{\xi} = (\xi_1, \dots, \xi_{n-1})$ normalized by $L_{\text{callable}}$:

* **Watterson's Estimator $\hat{\theta}_W$:**
  $$\hat{\theta}_W = \frac{S}{a_1 L_{\text{callable}}}, \qquad S = \sum_{i=1}^{n-1} \xi_i, \qquad a_1 = \sum_{i=1}^{n-1} \frac{1}{i}$$
  *Citation:* Watterson (1975) *Theor Popul Biol* [[DOI](https://doi.org/10.1016/0040-5809(75)90020-9)]
* **Pairwise Nucleotide Diversity $\hat{\pi}$:**
  $$\hat{\pi} = \frac{1}{L_{\text{callable}}} \sum_{i=1}^{n-1} \frac{2i(n-i)}{n(n-1)} \xi_i$$
  *Citation:* Nei & Li (1979) *Proc Natl Acad Sci USA* [[DOI](https://doi.org/10.1073/pnas.76.10.5269)]
* **Fay & Wu's $\hat{\theta}_H$ & Normalized $H$ Statistic:**
  $$\hat{\theta}_H = \frac{1}{L_{\text{callable}}} \sum_{i=1}^{n-1} \frac{2i^2}{n(n-1)} \xi_i, \qquad H = \hat{\pi} - \hat{\theta}_H$$
  *Citations:* Fay & Wu (2000) *Genetics* [[DOI](https://doi.org/10.1093/genetics/155.3.1405)]; Zeng et al. (2006) *Genetics* [[DOI](https://doi.org/10.1534/genetics.105.052605)]
* **Tajima's $D$ Neutrality Statistic:**
  $$D = \frac{\hat{\pi}_{\text{total}} - \hat{\theta}_{W,\text{total}}}{\sqrt{e_1 S + e_2 S(S-1)}}$$
  where $e_1 = c_1/a_1$, $e_2 = c_2/(a_1^2 + a_2)$, with exact analytical coefficients $b_1, b_2, c_1, c_2$ as formulated in Tajima (1989).  
  *Citation:* Tajima (1989) *Genetics* [[DOI](https://doi.org/10.1093/genetics/123.3.585)]  
  *Script:* [`scripts/compute_diversity_and_divergence.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/compute_diversity_and_divergence.py)

### 8.4. Population Differentiation & Resampling Under Autocorrelation

* **Hudson's Pairwise $F_{ST}$ (Ratio of Averages):**  
  To avoid the severe upward bias of Nei's $G_{ST}$ or Wright's raw $F_{ST}$ at low-diversity loci, we compute Hudson's estimator formulated strictly as the **ratio of averages**:
  $$\hat{F}_{ST}^{\text{Hudson}} = \frac{\sum_{k=1}^M \hat{N}_k}{\sum_{k=1}^M \hat{D}_k}$$
  $$\hat{N}_k = (\hat{p}_{1,k} - \hat{p}_{2,k})^2 - \frac{\hat{p}_{1,k}(1 - \hat{p}_{1,k})}{n_1 - 1} - \frac{\hat{p}_{2,k}(1 - \hat{p}_{2,k})}{n_2 - 1}$$
  $$\hat{D}_k = \hat{p}_{1,k}(1 - \hat{p}_{2,k}) + \hat{p}_{2,k}(1 - \hat{p}_{1,k})$$
  *Citations:* Hudson, Slatkin, & Maddison (1992) *Genetics* [[DOI](https://doi.org/10.1093/genetics/132.2.583)]; Bhatia et al. (2013) *Genome Res* [[DOI](https://doi.org/10.1101/gr.154831.113)]
* **Nei's Absolute ($d_{XY}$) and Net ($d_A$) Nucleotide Divergence:**
  $$d_{XY} = \frac{1}{L_{\text{callable}}} \sum_{k=1}^M \left[\hat{p}_{1,k}(1 - \hat{p}_{2,k}) + \hat{p}_{2,k}(1 - \hat{p}_{1,k})\right], \qquad d_A = d_{XY} - \frac{\hat{\pi}_1 + \hat{\pi}_2}{2}$$
  *Citation:* Nei & Li (1979) *Proc Natl Acad Sci USA* [[DOI](https://doi.org/10.1073/pnas.76.10.5269)]
* **Delete-$m_b$ Weighted Block Jackknife:**  
  Because physical linkage generates spatial autocorrelation across variant sites, standard i.i.d. bootstrap or naive normal approximations yield anti-conservative standard errors. We apply the delete-$m_b$ weighted block jackknife across $B = 24$ contiguous blocks ($5\text{ kb}$ window size) with unequal site weights $m_b$:
  $$\hat{\theta}_{\text{jack}} = B\hat{\theta} - \sum_{b=1}^B \left(1 - \frac{m_b}{M}\right)\hat{\theta}_{(-b)}, \qquad \widehat{\text{Var}}_{\text{jack}}(\hat{\theta}) = \frac{1}{B} \sum_{b=1}^B \frac{1}{h_b - 1} \left(\tilde{\tau}_b - \hat{\theta}_{\text{jack}}\right)^2$$
  where $h_b = M / m_b$, $\tilde{\tau}_b = h_b \hat{\theta} - (h_b - 1)\hat{\theta}_{(-b)}$, yielding asymptotically normal standard errors $\text{SE}_{\text{jack}} = \sqrt{\widehat{\text{Var}}_{\text{jack}}}$, $Z$-scores, and two-sided $p$-values.  
  *Citation:* Busing, Meijer, & Leeden (1999) *Stat Comput* [[DOI](https://doi.org/10.1023/A:1008900403334)]  
  *Script:* [`scripts/compute_diversity_and_divergence.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/compute_diversity_and_divergence.py)

### 8.5. Phylogenetics & Graph Topology Reconstruction

* **Kimura 2-Parameter (K80) & Jukes-Cantor (JC69) Evolutionary Distances:**  
  Pairwise genetic distances between chloroplast genomes are computed with pairwise gap deletion:
  $$d_{\text{K80}} = -\frac{1}{2}\ln(1 - 2P - Q) - \frac{1}{4}\ln(1 - 2Q)$$
  where $P$ is transition frequency and $Q$ is transversion frequency.  
  *Citations:* Kimura (1980) *J Mol Evol* [[DOI](https://doi.org/10.1007/BF01731581)]; Jukes & Cantor (1969) *Mamm Protein Metab* [[DOI](https://doi.org/10.1016/B978-1-4832-3211-9.50009-7)]
* **Neighbor-Joining (NJ) Algorithm with Non-Parametric Bootstrap:**  
  Phylogenetic trees are constructed using the Neighbor-Joining clustering algorithm and evaluated across 100 non-parametric bootstrap pseudo-replicates resampling polymorphic alignment columns. The tree is rooted on the basal *Oryza punctata* outgroup clade (BB genome).  
  *Citations:* Saitou & Nei (1987) *Mol Biol Evol* [[DOI](https://doi.org/10.1093/oxfordjournals.molbev.a040454)]; Felsenstein (1985) *Evolution* [[DOI](https://doi.org/10.1111/j.1558-5646.1985.tb00420.x)]
* **Median-Joining & Minimum Spanning Haplotype Networks:**  
  Maternal haplotype relationships are modeled as a Minimum Spanning Tree / Median-Joining graph topology where nodes represent unique maternal haplotypes (sized proportionally to accession frequency), edge lengths represent exact integer mutational steps, and color overlays denote varietal groups.  
  *Citation:* Bandelt, Forster, & Röhl (1999) *Mol Biol Evol* [[DOI](https://doi.org/10.1093/oxfordjournals.molbev.a026036)]  
  *Script:* [`scripts/reconstruct_phylogeny_and_network.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/reconstruct_phylogeny_and_network.py)

### 8.6. Linkage Disequilibrium Pruning, Kinship & Manifold Projections

* **Physical Sliding-Window LD Pruning ($r^2 < 0.2$):**  
  To eliminate collinear redundancies along the non-recombining plastome while preserving independent tag markers across compartments, variants with $\text{MAF} \ge 0.02$ are pruned within physical windows ($W_{\text{bp}} = 20,000\text{ bp}$) when pairwise squared Pearson correlation satisfies $r^2 \ge 0.20$. This collapses large structural rearrangements (such as the 12.5 kb SSC structural isomer) into independent tag variants.  
  *Script:* [`scripts/compute_population_structure.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/compute_population_structure.py)
* **Realized Genomic Relationship Matrix (Kinship $K$):**  
  Normalized covariance matrix computed via the VanRaden (2008) / Patterson et al. (2006) estimator:
  $$\mathbf{K} = \frac{1}{M} \mathbf{W}\mathbf{W}^T, \qquad W_{ik} = \frac{G_{ik} - p_k}{\sqrt{p_k(1 - p_k)}}$$
  where $\mathbf{K} \in \mathbb{R}^{N \times N}$, $\mathbf{K} = \mathbf{K}^T$, and $\mathbb{E}[\text{diag}(\mathbf{K})] = 1.0$.  
  *Citations:* VanRaden (2008) *J Dairy Sci* [[DOI](https://doi.org/10.3168/jds.2007-0980)]; Patterson, Price, & Reich (2006) *PLOS Genet* [[DOI](https://doi.org/10.1371/journal.pgen.0020190)]
* **Principal Component Analysis (PCA) & Eigendecomposition:**  
  Exact SVD decomposition of standardized genotypes extracting eigenvalues, variance explained percentages, and principal coordinates ($C = U\Sigma$).  
* **Non-Linear Manifold Projections (UMAP):**  
  Preservation of both local neighborhoods and global topological relationships in 2D and 3D manifold embeddings.  
  *Citation:* McInnes, Healy, & Melville (2018) *arXiv preprint arXiv:1802.03426*  
  *Script:* [`scripts/compute_population_structure.py`](file:///d:/Magister/Proyek/global-rice-popgen/scripts/compute_population_structure.py)

### 8.7. Association Mapping, Population Structure & Multiple Testing

* **Linear Mixed Models (LMM / MLM for GWAS):**  
  To control for cryptic relatedness and hierarchical population structure when mapping adaptive traits, phenotypes are modeled with fixed population covariates ($Q$ matrix or Principal Components) and random polygenic kinship ($K$):
  $$\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{S}\boldsymbol{\alpha} + \mathbf{Z}\mathbf{u} + \boldsymbol{\epsilon}, \qquad \mathbf{u} \sim \mathcal{N}(\mathbf{0}, \sigma_g^2 \mathbf{K}), \quad \boldsymbol{\epsilon} \sim \mathcal{N}(\mathbf{0}, \sigma_e^2 \mathbf{I})$$
  *Citations:* Kang et al. (2010) *Nat Genet* [[DOI](https://doi.org/10.1038/ng.548)]; Lippert et al. (2011) *Nat Methods* [[DOI](https://doi.org/10.1038/nmeth.1681)]
* **Genomic Inflation Factor ($\lambda_{GC}$) and False Discovery Rate (FDR):**  
  Calibration of genome-wide test statistics against unmodeled inflation:
  $$\lambda_{GC} = \frac{\text{median}(T_1^2, \dots, T_M^2)}{0.4549364}$$
  Multiple testing error control is enforced via the Benjamini-Hochberg step-up procedure at FDR $q < 0.05$.  
  *Citations:* Devlin & Roeder (1999) *Biometrics* [[DOI](https://doi.org/10.1111/j.0006-341x.1999.00997.x)]; Benjamini & Hochberg (1995) *J R Stat Soc B* [[DOI](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)]

---

## 9. Scientific Governance & Reproducibility Invariants

All analyses and software executions within this project must strictly comply with the governance rules defined in [`.agents/RULES/AGENTS.md`](file:///d:/Magister/Proyek/global-rice-popgen/.agents/RULES/AGENTS.md):

1. **Mandatory `.venv` Execution:** Every script and test runner executes exclusively inside the project's local virtual environment (`.venv/bin/python`, `pytest`) to avoid ABI drift and environment discrepancies.
2. **Strict Data Immutability:** Raw data (`data/raw/`) is read-only (chmod 444) and verified using SHA-256 cryptographic manifests. No pipeline step modifies raw inputs in-place.
3. **Denominator Integrity:** Nucleotide diversity ($\pi, \theta_W$) and divergence ($d_{XY}$) strictly account for all monomorphic and polymorphic callable positions ($L_{\text{callable}}$) identified in the BED mask. Dividing by genomic span is strictly prohibited.
4. **Ratio of Averages:** Windowed and genome-wide $F_{ST}$ estimates are calculated as the ratio of average numerators to average denominators ($\sum \hat{N} / \sum \hat{D}$) rather than the average of ratios.
5. **No MAF Pruning Prior to SFS/Neutrality:** Allele frequencies are never arbitrarily truncated ($\text{MAF} > 0.05$) prior to Site Frequency Spectrum, Tajima's $D$, or demographic history inference.
6. **PRNG Seeding & Determinism:** All stochastic methods (bootstrap pseudo-replicates, MCMC samplers, PCA) employ explicit 64-bit random seeds recorded in configuration manifests.

---

## 10. References & Literature Citations

### 10.1. Primary Methodological & Software References

1. **Alexander, D. H., Novembre, J., & Lange, K.** (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664. <https://doi.org/10.1101/gr.094052.109>
2. **Bandelt, H. J., Forster, P., & Röhl, A.** (1999). Median-joining networks for inferring intraspecific phylogenies. *Molecular Biology and Evolution*, 16(1), 37–48. <https://doi.org/10.1093/oxfordjournals.molbev.a026036>
3. **Benjamini, Y., & Hochberg, Y.** (1995). Controlling the false discovery rate: A practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society: Series B (Methodological)*, 57(1), 289–300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>
4. **Bhatia, G., Patterson, N., Sankararaman, S., & Price, A. L.** (2013). Estimating and interpreting $F_{ST}$: The impact of rare variants. *Genome Research*, 23(9), 1514–1521. <https://doi.org/10.1101/gr.154831.113>
5. **Busing, F. M. T. A., Meijer, E., & Van Der Leeden, R.** (1999). Delete-$m$ jackknife for unequal $m$. *Statistics and Computing*, 9(1), 3–8. <https://doi.org/10.1023/A:1008900403334>
6. **Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., ... & de Hoon, M. J. L.** (2009). Biopython: Freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. <https://doi.org/10.1093/bioinformatics/btp163>
7. **Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H.** (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. <https://doi.org/10.1093/gigascience/giab008>
8. **Devlin, B., & Roeder, K.** (1999). Genomic control for association studies. *Biometrics*, 55(4), 997–1004. <https://doi.org/10.1111/j.0006-341x.1999.00997.x>
9. **Fay, J. C., & Wu, C. I.** (2000). Hitchhiking under positive Darwinian selection. *Genetics*, 155(3), 1405–1413. <https://doi.org/10.1093/genetics/155.3.1405>
10. **Felsenstein, J.** (1985). Confidence limits on phylogenies: An approach using the bootstrap. *Evolution*, 39(4), 783–791. <https://doi.org/10.1111/j.1558-5646.1985.tb00420.x>
11. **Gotoh, O.** (1982). An improved algorithm for matching biological sequences. *Journal of Molecular Biology*, 162(3), 705–708. <https://doi.org/10.1016/0022-2836(82)90398-9>
12. **Hudson, R. R., Slatkin, M., & Maddison, W. P.** (1992). Estimation of levels of gene flow from DNA sequence data. *Genetics*, 132(2), 583–589. <https://doi.org/10.1093/genetics/132.2.583>
13. **Jukes, T. H., & Cantor, C. R.** (1969). Evolution of protein molecules. In H. N. Munro (Ed.), *Mammalian Protein Metabolism* (pp. 21–132). Academic Press. <https://doi.org/10.1016/B978-1-4832-3211-9.50009-7>
14. **Kang, H. M., Sul, J. H., Service, S. K., Zaitlen, N. A., Kong, S. Y., Freimer, N. B., ... & Eskin, E.** (2010). Variance component model to account for sample structure in genome-wide association studies. *Nature Genetics*, 42(4), 348–354. <https://doi.org/10.1038/ng.548>
15. **Katoh, K., & Standley, D. M.** (2013). MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. *Molecular Biology and Evolution*, 30(4), 772–780. <https://doi.org/10.1093/molbev/mst010>
16. **Kimura, M.** (1980). A simple method for estimating evolutionary rates of base substitutions through comparative studies of nucleotide sequences. *Journal of Molecular Evolution*, 16(2), 111–120. <https://doi.org/10.1007/BF01731581>
17. **Kosakovsky Pond, S. L., Poon, A. F. Y., Velazquez, R., Weaver, S., Hepler, N. L., Murrell, B., ... & Frost, S. D. W.** (2020). HyPhy 2.5—A customizable platform for evolutionary hypothesis testing using phylogenies. *Molecular Biology and Evolution*, 37(1), 295–299. <https://doi.org/10.1093/molbev/msaa159>
18. **Li, H.** (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv preprint arXiv:1303.3997*.
19. **Lippert, C., Listgarten, J., Liu, Y., Kadie, C. M., Davidson, R. I., & Heckerman, D.** (2011). FaST linear mixed models for genome-wide association studies. *Nature Methods*, 8(10), 833–835. <https://doi.org/10.1038/nmeth.1681>
20. **McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A.** (2010). The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Research*, 20(9), 1297–1303. <https://doi.org/10.1101/gr.107524.110>
21. **Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., ... & Köster, J.** (2021). Sustainable data analysis with Snakemake. *F1000Research*, 10, 33. <https://doi.org/10.12688/f1000research.29032.2>
22. **Nei, M., & Li, W. H.** (1979). Mathematical model for studying genetic variation in terms of restriction endonucleases. *Proceedings of the National Academy of Sciences of the United States of America*, 76(10), 5269–5273. <https://doi.org/10.1073/pnas.76.10.5269>
23. **Pedersen, B. S., & Quinlan, A. R.** (2018). Mosdepth: Quick coverage calculation for genomes and exomes. *Bioinformatics*, 34(5), 867–868. <https://doi.org/10.1093/bioinformatics/btx699>
24. **Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A., Bender, D., ... & Sham, P. C.** (2007). PLINK: A tool set for whole-genome association and population-based linkage analyses. *The American Journal of Human Genetics*, 81(3), 559–575. <https://doi.org/10.1086/519795>
25. **Raj, A., Stephens, M., & Pritchard, J. K.** (2014). fastSTRUCTURE: Variational inference of population structure in large SNP data sets. *Genetics*, 197(2), 573–589. <https://doi.org/10.1534/genetics.114.164350>
26. **Saitou, N., & Nei, M.** (1987). The neighbor-joining method: A new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406–425. <https://doi.org/10.1093/oxfordjournals.molbev.a040454>
27. **Tajima, F.** (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics*, 123(3), 585–595. <https://doi.org/10.1093/genetics/123.3.585>
28. **Vasimuddin, M., Misra, S., Li, H., & Aluru, S.** (2019). Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)*, 314–324. <https://doi.org/10.1109/IPDPS.2019.00041>
29. **Watterson, G. A.** (1975). On the number of segregating sites in genetical models without recombination. *Theoretical Population Biology*, 7(2), 256–276. <https://doi.org/10.1016/0040-5809(75)90020-9>
30. **Yang, Z.** (2007). PAML 4: Phylogenetic analysis by maximum likelihood. *Molecular Biology and Evolution*, 24(8), 1586–1591. <https://doi.org/10.1093/molbev/msm088>
31. **Zeng, K., Fu, Y. X., Shi, S., & Wu, C. I.** (2006). Statistical tests for detecting positive selection by utilizing high-frequency variants. *Genetics*, 174(3), 1431–1439. <https://doi.org/10.1534/genetics.105.052605>

### 10.2. BibTeX Bibliography (Copy-Ready for LaTeX / Overleaf / Zotero)

```bibtex
@article{katoh2013mafft,
  author    = {Katoh, Kazutaka and Standley, Daron M.},
  title     = {{MAFFT} multiple sequence alignment software version 7: improvements in performance and usability},
  journal   = {Molecular Biology and Evolution},
  volume    = {30},
  number    = {4},
  pages     = {772--780},
  year      = {2013},
  doi       = {10.1093/molbev/mst010}
}

@article{danecek2021samtools,
  author    = {Danecek, Petr and Bonfield, James K. and Liddle, Jennifer and Marshall, John and Ohan, Valeriu and Pollard, Martin O. and Whitwham, Andrew and Keane, Thomas and McCarthy, Shane A. and Davies, Robert M. and Li, Heng},
  title     = {Twelve years of {SAMtools} and {BCFtools}},
  journal   = {GigaScience},
  volume    = {10},
  number    = {2},
  pages     = {giab008},
  year      = {2021},
  doi       = {10.1093/gigascience/giab008}
}

@article{bhatia2013fst,
  author    = {Bhatia, Gaurav and Patterson, Nick and Sankararaman, Sriram and Price, Alkes L.},
  title     = {Estimating and interpreting ${F}_{ST}$: The impact of rare variants},
  journal   = {Genome Research},
  volume    = {23},
  number    = {9},
  pages     = {1514--1521},
  year      = {2013},
  doi       = {10.1101/gr.154831.113}
}

@article{tajima1989statistical,
  author    = {Tajima, Fumio},
  title     = {Statistical method for testing the neutral mutation hypothesis by {DNA} polymorphism},
  journal   = {Genetics},
  volume    = {123},
  number    = {3},
  pages     = {585--595},
  year      = {1989},
  doi       = {10.1093/genetics/123.3.585}
}

@article{watterson1975segregating,
  author    = {Watterson, G. A.},
  title     = {On the number of segregating sites in genetical models without recombination},
  journal   = {Theoretical Population Biology},
  volume    = {7},
  number    = {2},
  pages     = {256--276},
  year      = {1975},
  doi       = {10.1016/0040-5809(75)90020-9}
}

@article{nei1979mathematical,
  author    = {Nei, Masatoshi and Li, Wen-Hsiung},
  title     = {Mathematical model for studying genetic variation in terms of restriction endonucleases},
  journal   = {Proceedings of the National Academy of Sciences USA},
  volume    = {76},
  number    = {10},
  pages     = {5269--5273},
  year      = {1979},
  doi       = {10.1073/pnas.76.10.5269}
}

@article{fay2000hitchhiking,
  author    = {Fay, Justin C. and Wu, Chung-I},
  title     = {Hitchhiking under positive {Darwinian} selection},
  journal   = {Genetics},
  volume    = {155},
  number    = {3},
  pages     = {1405--1413},
  year      = {2000},
  doi       = {10.1093/genetics/155.3.1405}
}

@article{busing1999delete,
  author    = {Busing, Frank M. T. A. and Meijer, Erik and Van Der Leeden, Rien},
  title     = {Delete-$m$ jackknife for unequal $m$},
  journal   = {Statistics and Computing},
  volume    = {9},
  number    = {1},
  pages     = {3--8},
  year      = {1999},
  doi       = {10.1023/A:1008900403334}
}

@article{saitou1987neighbor,
  author    = {Saitou, Naruya and Nei, Masatoshi},
  title     = {The neighbor-joining method: a new method for reconstructing phylogenetic trees},
  journal   = {Molecular Biology and Evolution},
  volume    = {4},
  number    = {4},
  pages     = {406--425},
  year      = {1987},
  doi       = {10.1093/oxfordjournals.molbev.a040454}
}

@article{bandelt1999median,
  author    = {Bandelt, Hans-J{\"u}rgen and Forster, Peter and R{\"o}hl, Arne},
  title     = {Median-joining networks for inferring intraspecific phylogenies},
  journal   = {Molecular Biology and Evolution},
  volume    = {16},
  number    = {1},
  pages     = {37--48},
  year      = {1999},
  doi       = {10.1093/oxfordjournals.molbev.a026036}
}

@article{purcell2007plink,
  author    = {Purcell, Shaun and Neale, Benjamin and Todd-Brown, Kathe and Thomas, Lori and Ferreira, Manuel A. R. and Bender, David and Maller, Julian and Sklar, Pamela and de Bakker, Paul I. W. and Daly, Mark J. and Sham, Pak C.},
  title     = {{PLINK}: A tool set for whole-genome association and population-based linkage analyses},
  journal   = {The American Journal of Human Genetics},
  volume    = {81},
  number    = {3},
  pages     = {559--575},
  year      = {2007},
  doi       = {10.1086/519795}
}

@article{alexander2009admixture,
  author    = {Alexander, David H. and Novembre, John and Lange, Kenneth},
  title     = {Fast model-based estimation of ancestry in unrelated individuals},
  journal   = {Genome Research},
  volume    = {19},
  number    = {9},
  pages     = {1655--1664},
  year      = {2009},
  doi       = {10.1101/gr.094052.109}
}

@article{molder2021snakemake,
  author    = {M{\"o}lder, Felix and Jablonski, Kim Philipp and Letcher, Brice and Hall, Michael B. and Tomkins-Tinch, Christopher H. and Sochat, Vanessa and Forster, Jan and Lee, Sean and Twardziok, Sven O. and Kanitz, Alexander and Wilm, Andreas and Holtgrewe, Manuel and Rahmann, Sven and Nahnsen, Sven and K{\"o}ster, Johannes},
  title     = {Sustainable data analysis with {Snakemake}},
  journal   = {F1000Research},
  volume    = {10},
  pages     = {33},
  year      = {2021},
  doi       = {10.12688/f1000research.29032.2}
}
```
