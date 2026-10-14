# 3K Rice Genomes Project (3K RGP) Nuclear Genomics & Phenotypic Curation Report

**Project:** Global Rice Population Genomics (`global-rice-popgen`)  
**Cohort Scope:** 3,024 Accessions across 89 Countries (*Oryza sativa* L.)  
**Foundational Sources:** Wang et al. 2018 *Nature* (PMCID: `PMC6784863`), IRRI SNP-Seek, AWS Registry of Open Data  
**Generated Deliverables:**
- Passport Table: [`results/tables/rice_3k_passport_metadata.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/rice_3k_passport_metadata.tsv) | [`CSV`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/rice_3k_passport_metadata.csv)
- Phenotypic Table: [`results/tables/rice_3k_phenotype_metadata.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/rice_3k_phenotype_metadata.tsv) | [`CSV`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/rice_3k_phenotype_metadata.csv)
- Plastome-Nuclear Cross-Reference: [`results/tables/plastome_nuclear_crossref.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_nuclear_crossref.tsv) | [`CSV`](file:///d:/Magister/Proyek/global-rice-popgen/results/tables/plastome_nuclear_crossref.csv)
- Remote Download Manifest: [`data/raw/nuclear_3k/download_manifest_3k_rgp.tsv`](file:///d:/Magister/Proyek/global-rice-popgen/data/raw/nuclear_3k/download_manifest_3k_rgp.tsv)
- Cryptographic Manifest: [`data/raw/nuclear_3k/manifest.sha256`](file:///d:/Magister/Proyek/global-rice-popgen/data/raw/nuclear_3k/manifest.sha256)

---

## 1. Executive Summary & Scientific Purpose

Testing the project hypotheses ($H_1$: Plastid origin & cytonuclear discordance, $H_2$: Adaptive plastome evolution at range margins, and $H_3$: Cytonuclear epistasis & GWAS QTL mapping) requires connecting the maternal plastome phylogeny to a large-scale nuclear population genomic dataset. Asian rice (*Oryza sativa*) exhibits profound cytonuclear discordance resulting from ancient chloroplast capture and asymmetric inter-subspecies introgression between *indica* and *japonica* gene pools.

The **3,000 Rice Genomes Project (3K RGP)** provides the world standard reference cohort. This curation integrates the complete passport registry, the 9-component population admixture decomposition ($K=9$), the multi-season phenotypic evaluations, and establishes explicit links to the 90 curated chloroplast assemblies.

---

## 2. Ancestral Subpopulation Structure (K = 9)

The 3,024 accessions partition into the following primary and fine-scale nuclear sub-taxa based on high-density SNP analysis:

| Subpopulation Code | Varietal Group | Total Accessions | Percentage (%) | Evolutionary / Agronomic Characteristics |
| :--- | :--- | :---: | :---: | :--- |
| `Indica admixed` | Indica admixed | 615 | 20.34% | Lowland tropical indica gene pool, high tillering, submergence/heat tolerance |
| `Indica 3` | Indica 3 | 475 | 15.71% | Lowland tropical indica gene pool, high tillering, submergence/heat tolerance |
| `GJ-trp` | GJ-trp | 372 | 12.30% | Major varietal group |
| `Temperate japonica` | Temperate japonica | 287 | 9.49% | High-latitude cold-tolerant japonica adapted to East Asia and European margins |
| `Indica 2` | Indica 2 | 285 | 9.42% | Lowland tropical indica gene pool, high tillering, submergence/heat tolerance |
| `Indica 1A` | Indica 1A | 209 | 6.91% | Lowland tropical indica gene pool, high tillering, submergence/heat tolerance |
| `Indica 1B` | Indica 1B | 205 | 6.78% | Lowland tropical indica gene pool, high tillering, submergence/heat tolerance |
| `cA (Aus)` | cA (Aus) | 201 | 6.65% | Early-maturing, drought-tolerant, nitrogen-efficient South Asian lineage |
| `GJ-sbtrp` | GJ-sbtrp | 113 | 3.74% | Major varietal group |
| `admix` | admix | 103 | 3.41% | Major varietal group |
| `Japonica admixed` | Japonica admixed | 82 | 2.71% | Major varietal group |
| `cB (Bas)` | cB (Bas) | 76 | 2.51% | Major varietal group |
| `na` | na | 1 | 0.03% | Major varietal group |

---

## 3. Geographic Provenance Distribution

The germplasm encompasses 89 nations spanning all rice-growing agro-climatic zones:

| Geographic Region Code | Description | Accession Count | Percentage (%) |
| :--- | :--- | :---: | :---: |
| `SAC` | Other / Unspecified region | 486 | 16.07% |
| `EAS` | East Asia (China, Japan, South Korea) | 474 | 15.67% |
| `SEA` | Southeast Asia (Philippines, Thailand, Vietnam, Indonesia, Myanmar) | 468 | 15.48% |
| `SER` | Other / Unspecified region | 454 | 15.01% |
| `SAE` | Other / Unspecified region | 205 | 6.78% |
| `WAF` | Other / Unspecified region | 160 | 5.29% |
| `EAR` | East Asia Regional / Peripheral | 135 | 4.46% |
| `SER-IRRI` | Other / Unspecified region | 123 | 4.07% |
| `IOC` | Other / Unspecified region | 121 | 4.00% |
| `EUR` | Europe (Russia, Italy, Spain, Portugal, Greece) | 113 | 3.74% |
| `SAM` | Other / Unspecified region | 101 | 3.34% |
| `NAM` | Other / Unspecified region | 51 | 1.69% |
| `SAW` | Other / Unspecified region | 36 | 1.19% |
| `EAF` | Other / Unspecified region | 23 | 0.76% |
| `CAM` | Other / Unspecified region | 20 | 0.66% |
| `WAS` | West Asia (Iran, Iraq, Turkey) | 17 | 0.56% |
| `OCE` | Oceania (Australia, Papua New Guinea) | 16 | 0.53% |
| `na` | Other / Unspecified region | 12 | 0.40% |
| `NAF` | Other / Unspecified region | 9 | 0.30% |

---

## 4. Phenotypic Dataset Overview

- **Evaluated Accessions:** 2,266 varieties evaluated across uniform trial nurseries.
- **Key Phenotypic Dimensions for Hypotheses $H_2$ & $H_3$:**
  1. **Days to Heading / Flowering (`APANTH_REPRO`):** Critical for photoperiod and temperature sensitivity.
  2. **Seedling Height Class (`SDHT_CODE`):** Early seedling vigor, cold water tolerance at establishment.
  3. **Spikelet Fertility (`SPKF`):** Cold-induced floret sterility at heading; reproductive heat stress index.
  4. **Panicle Structure & Awn Presence (`PTY`, `APSH`, `SECOND_BR_REPRO`):** Domestication and adaptive harvest traits.
  5. **Seed Coat Pigmentation (`SCCO_REV`):** Flavonoid/proanthocyanidin accumulation associated with antioxidant defense.

---

## 5. Remote Download Manifest for Paired Nuclear SNP Datasets

Because nuclear SNP matrices encompass billions of genotypes across 3,024 samples (exceeding gigabytes in size), data streaming via AWS Open Data S3 (`s3://3kricegenome`) is pre-configured:

| Dataset Resource | Format | Target Scale | Access Path / URI |
| :--- | :--- | :--- | :--- |
| 3K RGP 4.8M Core Biallelic SNP Matrix (PLINK BED) | `PLINK Binary (.bed/.bim/.fam)` | ~4.8 Million SNPs | [`s3://3kricegenome/snpseek-dl/3krg-base-filt-core-v0.7/core_v0.7.bed.gz`](https://s3.amazonaws.com/3kricegenome/snpseek-dl/3krg-base-filt-core-v0.7/core_v0.7.bed.gz) |
| 3K RGP 4.8M Core SNP BIM Coordinate Table | `PLINK Variant Table (.bim)` | ~4.8 Million SNPs | [`s3://3kricegenome/snpseek-dl/3krg-base-filt-core-v0.7/core_v0.7.bim.gz`](https://s3.amazonaws.com/3kricegenome/snpseek-dl/3krg-base-filt-core-v0.7/core_v0.7.bim.gz) |
| 3K RGP 4.8M Core FAM Sample ID Table | `PLINK Sample Table (.fam)` | ~4.8 Million SNPs | [`s3://3kricegenome/snpseek-dl/3krg-base-filt-core-v0.7/core_v0.7.fam.gz`](https://s3.amazonaws.com/3kricegenome/snpseek-dl/3krg-base-filt-core-v0.7/core_v0.7.fam.gz) |
| 3K RGP 1M LD-Pruned SNP Matrix (pruned_v2.1) | `PLINK Binary (.bed/.bim/.fam)` | 1,000,000 SNPs | [`s3://3kricegenome/snpseek-dl/3k-pruned-v2.1/pruned_v2.1.bed`](https://3kricegenome.s3.dualstack.us-east-1.amazonaws.com/snpseek-dl/3k-pruned-v2.1/pruned_v2.1.bed) |
| 3K RGP Genomic Kinship Matrix (K) | `Bzip2 Matrix (.cXX.txt.bz2)` | 1M Pruned | [`s3://3kricegenome/snpseek-dl/3k-pruned-v2.1/pruned_v2.1.kinship/result.cXX.txt.bz2`](https://3kricegenome.s3.dualstack.us-east-1.amazonaws.com/snpseek-dl/3k-pruned-v2.1/pruned_v2.1.kinship/result.cXX.txt.bz2) |
| 3K RGP Ancestry Q-Matrix (K=9 Admixture) | `CSV (.csv)` | Pan-genome SNPs | [`s3://3kricegenome/snpseek-dl/Qmatrix-k9-3kRG.csv`](https://3kricegenome.s3.amazonaws.com/snpseek-dl/Qmatrix-k9-3kRG.csv) |
| 3K RGP Phenotypic Evaluation Matrix | `Excel Spreadsheet (.xlsx)` | 53 Evaluated Traits | [`s3://3kricegenome/snpseek-dl/phenotype/3kRG_PhenotypeData_v20170411.xlsx`](https://3kricegenome.s3.amazonaws.com/snpseek-dl/phenotype/3kRG_PhenotypeData_v20170411.xlsx) |
| 3K RGP Raw Sequence Read Archive (ENA) | `FASTQ (paired-end reads)` | Whole Genome Sequencing (~14x depth) | [`s3://3kricegenome/reads_AG43/`](https://www.ebi.ac.uk/ena/browser/view/PRJEB6180) |

---

## 6. Verification and Compliance

- **Data Immutability:** Raw download manifests and sample registries strictly housed in `data/raw/` with accompanying SHA-256 cryptographic manifests.
- **Deterministic Seeds:** Downstream population structure modeling and GWAS pipelines will reference sample indices locked in `results/tables/rice_3k_passport_metadata.tsv`.

