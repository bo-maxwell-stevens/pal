# GlobalAMFungi Data Prep Summary (Full)

## Input files
- GlobalAMFungi: `../Data/globalamf.csv`
- FRED: `../Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv`

## Dataset dimensions by stage

| stage | rows | unique_samples | unique_plant_species |
|---|---:|---:|---:|
| initial | 2171488 | 5524 | 441 |
| root_only | 1096982 | 2691 | 351 |
| root_singleplant_nonnull | 1079452 | 2521 | 285 |

## Richness statistics (sample-level)

| metric | count | mean | std | min | 25% | 50% | 75% | max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| amf_seq_richness | 2521 | 261.9 | 983.3 | 1 | 4 | 18 | 120 | 2.048e+04 |
| amf_genus_richness | 2521 | 2.011 | 1.319 | 1 | 1 | 2 | 3 | 10 |
| total_reads | 2521 | 4.178e+06 | 4.455e+07 | 1 | 19 | 761 | 6.076e+04 | 1.346e+09 |
| MAT | 2521 | 11.44 | 7.321 | -5.9 | 7.3 | 10.2 | 16 | 28.7 |
| MAP | 2521 | 898.1 | 630.5 | 78 | 517 | 657 | 1197 | 4062 |
| pH | 2112 | 6.038 | 1.27 | 2.79 | 5.18 | 5.98 | 7.05 | 8.76 |

## FRED matching statistics

- total_global_species: **285**
- total_fred_species: **896**
- exact_matches_raw: **46**
- matches_canonical: **46**
- genus_overlap: **96**
- match_percentage: **16.14%**

## Warnings / caveats

- None

## Exact printed notebook output

```text
Global input: ../Data/globalamf.csv
FRED input: ../Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv
Output directory: /home/stevens/projects/pal/Output
--- Initial GlobalAMFungi summary ---
rows: 2,171,488
unique samples: 5,524
unique plant species: 441
initial: rows=2,171,488, unique_samples=5,524, unique_plant_species=441
root_only: rows=1,096,982, unique_samples=2,691, unique_plant_species=351
root_singleplant_nonnull: rows=1,079,452, unique_samples=2,521, unique_plant_species=285
--- Sample-level summary ---
rows: 2,521
mean amf_seq_richness: 261.934
median amf_seq_richness: 18.000
mean total_reads: 4178282.573
Wrote: ../Output/globalamfungi_sample_level_full.csv
--- Plant-species-level summary ---
total plant species: 285
total samples: 2,521
Wrote: ../Output/globalamfungi_species_level_full.csv
FRED rows (data-only): 6,271
FRED unique species (canonical): 896
--- GlobalAMFungi ↔ FRED matching ---
global_unique_species_after_filters: 285
fred_unique_species: 896
exact_matches_raw: 46
canonical_matches: 46
genus_overlap: 96
match_percentage: 16.14%
Wrote figure: ../Output/globalamf_hist_amf_seq_richness_full.png
Wrote figure: ../Output/globalamf_scatter_richness_vs_MAT_full.png
Wrote figure: ../Output/globalamf_map_points_full.png
Wrote markdown summary: ../Output/globalamfungi_data_prep_summary_full.md
Final dataset summary:
Total samples: 2,521
Total plant species: 285
Total GlobalAMFungi–FRED matches: 46
Match percentage: 16.14%
```