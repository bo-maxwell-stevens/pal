# FRED-AMF25 Trait-Richness Results Summary

## Dataset overview
- Total samples: 5,226
- Samples per source: {'EcoBank': 2705, 'GlobalAMFungi': 2521}
- Total species: 568
- Total genera: 387

## Trait coverage (fraction non-missing)

| trait | fraction_non_missing |
|---|---:|
| Root diameter | 0.144 |
| Root tissue density (RTD) | 0.095 |
| Specific root area (SRA) | 0.051 |
| Root N content | 0.160 |
| Root P content | 0.087 |
| Mycorrhiza_Fraction of root length or tips colonized | 0.031 |

## Model results

| model_type | trait | beta | se | p | R2 | N |
|---|---|---:|---:|---:|---:|---:|
| seq_richness | Root diameter | -0.2873 | 0.1007 | 0.004325 | 0.943 | 750 |
| genus_richness | Root diameter | -0.3043 | 0.111 | 0.0061 | 0.379 | 481 |
| seq_richness | Root tissue density (RTD) | -0.03968 | 0.1024 | 0.6983 | 0.9409 | 495 |
| genus_richness | Root tissue density (RTD) | 0.07475 | 0.1157 | 0.5183 | 0.4308 | 304 |
| seq_richness | Specific root area (SRA) | -8.776e-05 | 0.000134 | 0.5126 | 0.9321 | 264 |
| genus_richness | Specific root area (SRA) | 0.0001568 | 0.0002053 | 0.4451 | 0.4905 | 177 |
| seq_richness | Root N content | -0.007574 | 0.003259 | 0.02012 | 0.9586 | 834 |
| genus_richness | Root N content | -0.007053 | 0.003069 | 0.02154 | 0.2678 | 559 |
| seq_richness | Root P content | -0.03878 | 0.02975 | 0.1924 | 0.9565 | 457 |
| genus_richness | Root P content | 0.1151 | 0.04979 | 0.02076 | 0.2845 | 343 |
| interaction_woodiness | Root diameter | 0.4707 | 0.1717 | 0.00613 | 0.7774 | 750 |
| interaction_woodiness | Root tissue density (RTD) | -1.218 | 0.1925 | 2.535e-10 | 0.7663 | 495 |
| interaction_woodiness | Specific root area (SRA) | -0.0003047 | 0.000241 | 0.2061 | 0.8046 | 264 |
| interaction_woodiness | Root N content | 0.04077 | 0.01748 | 0.01965 | 0.8617 | 834 |
| interaction_woodiness | Root P content | 0.052 | 0.05465 | 0.3414 | 0.8775 | 457 |

## Interpretation relative to RES hypotheses

- Root diameter: negative association (beta=-0.287, p=0.00432)
- Root N content: negative association (beta=-0.00757, p=0.0201)
- These patterns are exploratory and not fully phylogenetically controlled; they should be interpreted as broad associations.
- Consistency checks across sample-level, species-level, and woody/nonwoody plots provide qualitative support where signs agree.

## Printed output log

```text
Input dataset: ../Output/amf_traits_merged_sample_level.csv
Figure directory: /home/stevens/projects/pal/Output/Figures
Summary markdown: /home/stevens/projects/pal/Output/fred_amf25_results_summary.md
Total rows: 5,226
Sources: {'EcoBank': 2705, 'GlobalAMFungi': 2521}
Unique species: 568
Unique genera: 387
Wrote figure: ../Output/Figures/Fig01_global_sample_map.png
Wrote figure: ../Output/Figures/Fig02_richness_histogram.png
Wrote figure: ../Output/Figures/Fig03_log_reads_histogram.png
Wrote figure: ../Output/Figures/Fig04_RD_vs_richness_sample.png
Wrote figure: ../Output/Figures/Fig05_RTD_vs_richness_sample.png
Wrote figure: ../Output/Figures/Fig06_SRA_vs_richness_sample.png
Wrote figure: ../Output/Figures/Fig07_RootN_vs_richness_sample.png
Wrote figure: ../Output/Figures/Fig08_RootP_vs_richness_sample.png
Wrote figure: ../Output/Figures/Fig09_RD_vs_richness_species.png
Wrote figure: ../Output/Figures/Fig10_RTD_vs_richness_species.png
Wrote figure: ../Output/Figures/Fig11_SRA_vs_richness_species.png
Wrote figure: ../Output/Figures/Fig12_RootN_vs_richness_species.png
Wrote figure: ../Output/Figures/Fig13_RootP_vs_richness_species.png
Wrote figure: ../Output/Figures/Fig14_RES_PCA.png
PCA genera used: 2
Wrote figure: ../Output/Figures/Fig15_RD_vs_richness_woody.png
Wrote figure: ../Output/Figures/Fig16_RD_vs_richness_nonwoody.png
Wrote figure: ../Output/Figures/Fig17_RTD_vs_richness_woody.png
Wrote figure: ../Output/Figures/Fig18_RTD_vs_richness_nonwoody.png
Wrote figure: ../Output/Figures/Fig19_SRA_vs_richness_woody.png
Wrote figure: ../Output/Figures/Fig19b_SRA_vs_richness_nonwoody.png
Wrote figure: ../Output/Figures/Fig19c_RootN_vs_richness_woody.png
Wrote figure: ../Output/Figures/Fig19d_RootN_vs_richness_nonwoody.png
Wrote figure: ../Output/Figures/Fig19e_RootP_vs_richness_woody.png
Wrote figure: ../Output/Figures/Fig19f_RootP_vs_richness_nonwoody.png
Wrote figure: ../Output/Figures/Fig20_MAT_vs_richness.png
Wrote figure: ../Output/Figures/Fig21_MAP_vs_richness.png
Model rows generated: 15
```