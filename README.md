# DNA methylation and the frailty index (FI): preprocessing, EWAS and downstream analyses

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24238651.svg)](https://doi.org/10.5281/zenodo.24238651)

R code for the analyses of blood DNA methylation (Illumina Infinium MethylationEPIC v1) and the
frailty index (FI) described in the accompanying manuscript.

**This repository contains code only. No participant data, sample identifiers, clinical information
or analysis outputs are included.** All file paths in the scripts are placeholders (see
[Configuration](#configuration)).

## Repository structure

Scripts are numbered in the order of the sections of the manuscript's Methods.

```
.
├── R/
│   ├── 01_preprocessing_EWAS.R      # Genome-wide DNA methylation profiling
│   ├── 02_sensitivity_analysis.R    # Sensitivity analyses
│   ├── 03_replication_Mak2024.R     # Overlap with previously reported frailty-associated CpGs
│   ├── 04_genomic_distribution.R    # Genomic distribution and enrichment analyses (distribution)
│   ├── 05_enrichment_missMethyl.R   # Genomic distribution and enrichment analyses (enrichment)
│   ├── 06_epigenetic_clocks.R       # Estimation of epigenetic age
│   └── 07_nested_cv_FI.R            # FI prediction (nested cross-validation)
├── nested_cv/                       # generic nested-CV implementation (Zenodo v1.0.0, unchanged)
│   ├── nested_cv.R
│   ├── run_analysis.R
│   ├── tests/test_synthetic.R
│   ├── README.md
│   └── VALIDATION.md
├── .gitignore
├── LICENSE
└── README.md
```

## Workflow

`01` must be run first. Scripts `02`–`07` only read the outputs of `01` and can be run in any order.

| Manuscript section | Script | Main inputs | Main outputs (under `Results_directory/`) |
|---|---|---|---|
| Genome-wide DNA methylation profiling | `01_preprocessing_EWAS.R` | IDATs, sample sheet, SNP / cross-reactive probe lists | QC tables, cell proportions, filtered beta / M-value matrices, primary and secondary EWAS results |
| Sensitivity analyses | `02_sensitivity_analysis.R` | M-values, cell proportions, primary EWAS, clinical information | `sensitivity/`: models with additional covariates, concordance tables |
| Overlap with previously reported frailty-associated CpGs | `03_replication_Mak2024.R` | primary and secondary EWAS, Mak et al. Supplement 2 | `replication/`: concordance tables, replication summary |
| Genomic distribution and enrichment analyses | `04_genomic_distribution.R` | primary EWAS | `genomic_features/`: category proportions, Hyper vs Hypo tests, bar plots |
| | `05_enrichment_missMethyl.R` | primary EWAS, Reactome / GO BP GMT files | `enrichment/`: pathway summaries, gene tables, dot plots |
| Estimation of epigenetic age | `06_epigenetic_clocks.R` | ssNoob + BMIQ beta, phenotypes, PC clock coefficients | `epigenetic_clocks/`: clock estimates, EAA, QC tables, FI correlations |
| FI prediction | `07_nested_cv_FI.R` | ssNoob + BMIQ beta, phenotypes | `nested_cv_FI/`: out-of-fold predictions, metrics, selected CpGs |

### Analysis summary

1. **Genome-wide DNA methylation profiling (01).** IDATs are read with minfi. Samples whose predicted
   sex does not match the recorded sex are excluded. Probes are excluded if detection P ≥ 0.05 in any
   sample, if they overlap known SNPs (published lists and `minfi::dropLociWithSnps`, MAF > 0) or if
   they are cross-reactive. Beta values (M / (M + U)) are normalized with BMIQ (wateRmelon) and
   converted to M-values. Cell proportions are estimated with `estimateCellCounts2` and the EPIC IDOL
   reference (FlowSorted.Blood.EPIC). EWAS with limma:
   `M ~ FI + age + sex + CD8T + CD4T + NK + Bcell + Mono` (neutrophils = reference), Benjamini–Hochberg
   FDR; secondary model with frailty status (Non-frail / Pre-frail / frail, reference = Non-frail).
   A single-sample Noob (ssNoob) + BMIQ dataset with the same samples and probes is saved for `06` and `07`.
2. **Sensitivity analyses (02).** Smoking, diabetes and polypharmacy are each added to the primary model;
   each model is compared with the base model fitted on the same samples.
3. **Overlap with previously reported CpGs (03).** CpGs reported by Mak et al. (2024) for continuous FI
   (Table S3) and frail vs non-frail (Table S4) are compared with the primary and secondary EWAS
   (replication: P < 0.05 and same direction).
4. **Genomic distribution and enrichment (04, 05).** FDR < 0.1 CpGs, split into hyper- and
   hypomethylated, are compared across CpG island, DNase hypersensitivity and gene-region categories
   (04). Gene set testing uses `missMethyl::gsameth` on promoter CpGs with `prior.prob = TRUE` (05).
5. **Estimation of epigenetic age (06).** methylclock (Horvath, Hannum, Skin & Blood, PhenoAge),
   meffonym (GrimAge, GrimAge v2), PC clocks (methylCIPHER) and DunedinPACE. EAA is the residual of
   each clock regressed on chronological age; DunedinPACE is used without residualization. FI is
   correlated with each measure using Pearson correlation.
6. **FI prediction (07).** Nested cross-validation with training-fold-only limma screening of the top
   800 CpGs and elastic net regression; see [`nested_cv/README.md`](nested_cv/README.md).

## Software

Analyses were run in R 4.6.1. Exact package versions are written to output files at run time
(`minfi` version message in `01`, `package_versions.txt` in `06`, `sessionInfo.txt` in `07`).

| Package | Used for | Script(s) | Source |
|---|---|---|---|
| minfi | IDAT loading, detection P, sex prediction, ssNoob, SNP filtering | 01 | [GitHub](https://github.com/hansenlab/minfi) · [Bioconductor](https://bioconductor.org/packages/minfi) |
| FlowSorted.Blood.EPIC | `estimateCellCounts2`, EPIC IDOL reference | 01 | [GitHub](https://github.com/immunomethylomics/FlowSorted.Blood.EPIC) · [Bioconductor](https://bioconductor.org/packages/FlowSorted.Blood.EPIC) |
| wateRmelon | BMIQ normalization | 01 | [GitHub](https://github.com/schalkwyk/wateRmelon) · [Bioconductor mirror](https://github.com/bioc/wateRmelon) · [Bioconductor](https://bioconductor.org/packages/wateRmelon) |
| limma | EWAS (moderated statistics), nested-CV screening | 01, 02, 07 | [Bioconductor](https://bioconductor.org/packages/limma) |
| IlluminaHumanMethylationEPICmanifest, IlluminaHumanMethylationEPICanno.ilm10b4.hg19 | EPIC v1 manifest and annotation (hg19) | 01, 04, 05 | [Bioconductor](https://bioconductor.org/packages/IlluminaHumanMethylationEPICanno.ilm10b4.hg19) |
| missMethyl | gene set testing (`gsameth`) | 05 | [GitHub](https://github.com/Oshlack/missMethyl) · [Bioconductor](https://bioconductor.org/packages/missMethyl) |
| org.Hs.eg.db | gene symbol to Entrez ID mapping | 05 | [Bioconductor](https://bioconductor.org/packages/org.Hs.eg.db) |
| methylclock | Horvath, Hannum, Skin & Blood, PhenoAge | 06 | [GitHub](https://github.com/isglobal-brge/methylclock) · [Bioconductor](https://bioconductor.org/packages/methylclock) |
| meffonym | GrimAge, GrimAge v2 (and cross-check of other clocks) | 06 | [GitHub](https://github.com/perishky/meffonym) |
| methylCIPHER | PC clocks (`calcPCClocks`) | 06 | [GitHub](https://github.com/MorganLevineLab/methylCIPHER) · [PC-Clocks](https://github.com/MorganLevineLab/PC-Clocks) |
| DunedinPACE | pace of aging | 06 | [GitHub](https://github.com/danbelsky/DunedinPACE) |
| glmnet | elastic net | 07 | [CRAN](https://cran.r-project.org/package=glmnet) · [website](https://glmnet.stanford.edu/) |
| readxl, openxlsx, ggplot2 | file input, plots | 01–05 | CRAN |

```r
## CRAN
install.packages(c("readxl", "openxlsx", "ggplot2", "glmnet", "devtools", "BiocManager"))

## Bioconductor
BiocManager::install(c(
  "minfi", "IlluminaHumanMethylationEPICmanifest",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "FlowSorted.Blood.EPIC",
  "wateRmelon", "limma", "missMethyl", "org.Hs.eg.db", "methylclock"
))

## GitHub
devtools::install_github("perishky/meffonym")
devtools::install_github("MorganLevineLab/methylCIPHER")
devtools::install_github("danbelsky/DunedinPACE")
```

## External resources (not redistributed here)

Download these separately and point the corresponding `cfg` entries to them:

| Resource | Script | `cfg` entry |
|---|---|---|
| SNP-overlapping probe lists (sheets `CpG_rs`, `SBE_rs`, `Probe_rs`, `stem_cg`) and cross-reactive probe list for EPIC, e.g. the annotated lists in the Additional files of [Pidsley et al. (2016) *Genome Biology* 17:208](https://doi.org/10.1186/s13059-016-1066-1) | 01 | `snp_file`, `cross_files` |
| Mak et al. (2024) *Aging Cell* 23:e14135, Supplement 2 (Tables S3 and S4), [doi:10.1111/acel.14135](https://doi.org/10.1111/acel.14135) | 03 | `rep_xlsx` |
| Enrichr gene set libraries `Reactome_Pathways_2024` and `GO_Biological_Process_2026` (GMT), [Enrichr libraries](https://maayanlab.cloud/Enrichr/#libraries) | 05 | `reactome_gmt`, `gobp_gmt` |
| PC clock coefficients (`CalcAllPCClocks.RData`), see [methylCIPHER](https://github.com/MorganLevineLab/methylCIPHER) | 06 | `pcclock_dir` |

## Configuration

Every script has a `cfg` list at the top. Paths are the only thing that needs to be changed.
Placeholders:

| Placeholder | Meaning |
|---|---|
| `idat_directory` | folder containing the IDAT files |
| `clinical_sheet_file` | sample sheet with phenotypes (xlsx), used in 01 |
| `clinical_information_file` | clinical information with additional covariates (xlsx), used in 02 |
| `filtering_option` | SNP probe list workbook |
| `cross_active_probe` | cross-reactive probe list workbook (column `probe`) |
| `Results_directory` | output folder of 01; every other script reads from and writes into it (`v2_dir`) |
| `paper_sensitivity_rds` | (optional) previously generated sensitivity results, read-only comparison in 02 |
| `replication_supplement_xlsx` | Mak et al. (2024) Supplement 2 |
| `reactome_gmt_file`, `gobp_gmt_file` | Enrichr GMT files |
| `pcclock_directory` | folder containing `CalcAllPCClocks.RData` |

### Expected input format

**Sample sheet (`clinical_sheet_file`, 01).** Illumina sample sheet layout: six lines are skipped and
the next row holds the column names. One row per array, with at least:

| Column | Content |
|---|---|
| `Sentrix_ID`, `Sentrix_Position` | used to build IDAT basenames (`<Sentrix_ID>_<Sentrix_Position>`) |
| `Sample_Name` | sample label |
| `Sex` | `1` = male, `0` = female |
| `Age` | years |
| `FI` | frailty index (continuous) |
| `Frailty_status` | `Non-frail`, `Pre-frail`, `frail` (labels must match `frail_levels` exactly) |

**Clinical information (`clinical_information_file`, 02).** One row per sample with `Sentrix_id`
(matching `<Sentrix_ID>_<Sentrix_Position>`), `FI`, `Age`, `Sex` and the additional covariates
`Smoking` (`9` or blank = missing), `Diabetes` and `Polypharmacy`.

Column names can be changed in `cfg` (`col_FI`, `col_age`, `col_sex`, `col_frail`, `id_col`, `sens`).

## Data availability and privacy

- Do **not** commit IDATs, clinical sheets or anything under `Results_directory/`. Outputs contain
  sample identifiers and sample-level values (e.g. `pheno_with_cellprop.rds`, `epigenetic_clocks.txt`,
  nested-CV predictions). The provided `.gitignore` excludes these file types.
- Participant-level data are not publicly available; see the Data Availability statement of the
  manuscript.

## Third-party code

`nested_cv/` is the generic nested cross-validation implementation archived on Zenodo
(record [24211256](https://zenodo.org/records/24211256), `generic_nested_cv_v1.0.0.zip`,
CC BY 4.0), included without modification. Its own README and validation notes are in that folder;
run `Rscript nested_cv/tests/test_synthetic.R` for the synthetic self-test.

## License

All code in this repository, including `nested_cv/`, is released under the
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](LICENSE).
Copyright (c) 2026 BIOJHK, except `nested_cv/` (see [Third-party code](#third-party-code)).

## Citation

If you use this code, please cite the accompanying manuscript (citation to be added) and this
code release:

- Code (v1.0.0): Zenodo, [doi:10.5281/zenodo.24238651](https://doi.org/10.5281/zenodo.24238651)
- Nested cross-validation: the Zenodo record listed under [Third-party code](#third-party-code)
