# FRED-AMF25 Results Summary (FIXED)

## Why prior OLS was problematic and what was fixed
- Prior OLS treated sample rows as independent despite clustering by `study_id` and repeated species-level traits, inflating nominal significance.
- This analysis uses mixed-effects models with random intercept structure (study + species variance component) where feasible, with cluster-robust fallback when mixed fitting fails.
- Species-level weighted models (weights=`n_samples`) provide a second inference layer aligned with species-level trait measurement.
- Richness was residualized against sequencing depth (`log_total_reads`) and source to reduce confounding.

## QC and coverage
- Rows: 5,226
- Unique samples: 5,226
- Unique species: 568
- Unique studies: 976
- Rows by source: {'EcoBank': 2705, 'GlobalAMFungi': 2521}

### Missingness (overall)
| variable | n_missing | frac_missing |
|---|---:|---:|
| Root diameter | 4476 | 0.856 |
| Root tissue density (RTD) | 4731 | 0.905 |
| Specific root area (SRA) | 4962 | 0.949 |
| Root N content | 4392 | 0.840 |
| Root P content | 4769 | 0.913 |
| Mycorrhiza_Fraction of root length or tips colonized | 5062 | 0.969 |
| MAT | 0 | 0.000 |
| MAP | 0 | 0.000 |
| log_total_reads | 1 | 0.000 |

### Coverage by source (fraction non-missing)
| source | variable | frac_non_missing |
|---|---|---:|
| EcoBank | Root diameter | 0.099 |
| EcoBank | Root tissue density (RTD) | 0.071 |
| EcoBank | Specific root area (SRA) | 0.032 |
| EcoBank | Root N content | 0.102 |
| EcoBank | Root P content | 0.042 |
| EcoBank | Mycorrhiza_Fraction of root length or tips colonized | 0.007 |
| EcoBank | MAT | 1.000 |
| EcoBank | MAP | 1.000 |
| EcoBank | log_total_reads | 1.000 |
| GlobalAMFungi | Root diameter | 0.191 |
| GlobalAMFungi | Root tissue density (RTD) | 0.121 |
| GlobalAMFungi | Specific root area (SRA) | 0.070 |
| GlobalAMFungi | Root N content | 0.222 |
| GlobalAMFungi | Root P content | 0.136 |
| GlobalAMFungi | Mycorrhiza_Fraction of root length or tips colonized | 0.057 |
| GlobalAMFungi | MAT | 1.000 |
| GlobalAMFungi | MAP | 1.000 |
| GlobalAMFungi | log_total_reads | 1.000 |

## Sequencing depth confounding
- Depth model (seq richness) R2: 0.970
- Depth model (genus richness) R2: 0.222
- Residual outcomes (`richness_resid_seq`, `richness_resid_genus`) were used in sensitivity models and figures to reduce read-depth bias.

## Main results (Root diameter and Root N)
- Sample-level mixed/robust (pooled), Root diameter -> seq richness: beta=-0.2873, p=0.1427, N=750
- Sample-level mixed/robust (pooled), Root N -> seq richness: beta=-0.006415, p=0.3078, N=834
- Species-level weighted (pooled), Root diameter -> seq richness: beta=0.1304, p=0.8904, N=34
- Species-level weighted (pooled), Root N -> seq richness: beta=-0.005243, p=0.8427, N=37
- Persistence across pooled mixed/robust models, species-level weighted models, and source-stratified models should be interpreted as stronger evidence than any single model alone.

## Woodiness interactions
| trait | beta_main_trait | p_value | N |
|---|---:|---:|---:|
| Root diameter | -0.1436 | 0.1429 | 750 |
| Root tissue density (RTD) | -0.01984 | 0.8631 | 495 |
| Specific root area (SRA) | -4.388e-05 | 0.5879 | 264 |
| Root N content | -1.006 | 1 | 834 |
| Root P content | -0.01939 | 0.3222 | 457 |

## PCA validity
- PCA-2trait completed on species level with N=15.
- PCA-RES4 skipped due to insufficient coverage (species N=3, genus N=2).

## Model result files
- `../Output/fred_amf25_model_results_FIXED.csv`
- `../Output/fred_amf25_model_results_specieslevel_FIXED.csv`

## Figure inventory (_FIXED)
- `../Output/Figures/Fig01_global_sample_map_FIXED.png`
- `../Output/Figures/Fig02_richness_histogram_FIXED.png`
- `../Output/Figures/Fig03_log_reads_histogram_by_source_FIXED.png`
- `../Output/Figures/Fig04_RD_vs_richness_sample_FIXED.png`
- `../Output/Figures/Fig05_RTD_vs_richness_sample_FIXED.png`
- `../Output/Figures/Fig06_SRA_vs_richness_sample_FIXED.png`
- `../Output/Figures/Fig07_RootN_vs_richness_sample_FIXED.png`
- `../Output/Figures/Fig08_RootP_vs_richness_sample_FIXED.png`
- `../Output/Figures/Fig09_RD_vs_richness_species_FIXED.png`
- `../Output/Figures/Fig10_RTD_vs_richness_species_FIXED.png`
- `../Output/Figures/Fig11_SRA_vs_richness_species_FIXED.png`
- `../Output/Figures/Fig12_RootN_vs_richness_species_FIXED.png`
- `../Output/Figures/Fig13_RootP_vs_richness_species_FIXED.png`
- `../Output/Figures/Fig15_RD_woody_nonwoody_species_FIXED.png`
- `../Output/Figures/Fig16_RootN_woody_nonwoody_species_FIXED.png`
- `../Output/Figures/Fig20_MAT_vs_richness_resid_FIXED.png`
- `../Output/Figures/Fig21_MAP_vs_richness_resid_FIXED.png`
- `../Output/Figures/Fig14_PCA_2trait_FIXED.png`

## Limitations
- Trait missingness remains high for several RES dimensions, limiting precision and some multivariate analyses.
- Species-level traits are reused across samples; mixed-effects and species-level analyses reduce, but do not eliminate, dependence concerns.
- Plant phylogeny is not explicitly modeled; trait effects may partially proxy phylogenetic structure.
- EcoBank coverage and sequencing characteristics differ from GlobalAMFungi; source-stratified results are therefore critical context.

## Notebook log
```text
Reference found: FRED.Rmd
Reference found: FRED-AMF25 (1).docx
Rows: 5,226
Unique samples: 5,226
Unique species: 568
Unique studies: 976
Rows by source: {'EcoBank': 2705, 'GlobalAMFungi': 2521}
Depth models fitted.
Seq depth model R2: 0.970
Genus depth model R2: 0.222
Wrote: ../Output/fred_amf25_model_results_FIXED.csv
Sample-level model rows: 35
Wrote: ../Output/fred_amf25_model_results_specieslevel_FIXED.csv
Species-level model rows: 30
Wrote figure: ../Output/Figures/Fig01_global_sample_map_FIXED.png
Wrote figure: ../Output/Figures/Fig02_richness_histogram_FIXED.png
Wrote figure: ../Output/Figures/Fig03_log_reads_histogram_by_source_FIXED.png
Wrote figure: ../Output/Figures/Fig04_RD_vs_richness_sample_FIXED.png
Wrote figure: ../Output/Figures/Fig05_RTD_vs_richness_sample_FIXED.png
Wrote figure: ../Output/Figures/Fig06_SRA_vs_richness_sample_FIXED.png
Wrote figure: ../Output/Figures/Fig07_RootN_vs_richness_sample_FIXED.png
Wrote figure: ../Output/Figures/Fig08_RootP_vs_richness_sample_FIXED.png
Wrote figure: ../Output/Figures/Fig09_RD_vs_richness_species_FIXED.png
Wrote figure: ../Output/Figures/Fig10_RTD_vs_richness_species_FIXED.png
Wrote figure: ../Output/Figures/Fig11_SRA_vs_richness_species_FIXED.png
Wrote figure: ../Output/Figures/Fig12_RootN_vs_richness_species_FIXED.png
Wrote figure: ../Output/Figures/Fig13_RootP_vs_richness_species_FIXED.png
Wrote figure: ../Output/Figures/Fig15_RD_woody_nonwoody_species_FIXED.png
Wrote figure: ../Output/Figures/Fig16_RootN_woody_nonwoody_species_FIXED.png
Wrote figure: ../Output/Figures/Fig20_MAT_vs_richness_resid_FIXED.png
Wrote figure: ../Output/Figures/Fig21_MAP_vs_richness_resid_FIXED.png
Generated fixed figures: 17
Wrote figure: ../Output/Figures/Fig14_PCA_2trait_FIXED.png
PCA-2trait completed on species level with N=15.
PCA-RES4 skipped due to insufficient coverage (species N=3, genus N=2).
Root diameter: nonwoody had <5 species; panel not plotted.
Root N content: nonwoody had <5 species; panel not plotted.
```