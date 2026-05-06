# GlobalAMFungi Data Prep Summary

## Inputs
- GlobalAMFungi: `../Data/globalamf.csv`
- FRED: `../Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv`

## Dataset dimensions at each filter step

| stage | rows | unique_samples | unique_plant_species |
|---|---:|---:|---:|
| initial | 2171488 | 5524 | 441 |
| root_only | 1096982 | 2691 | 351 |
| root_singleplant | 1079452 | 2521 | 285 |
| root_singleplant_singleton_species | 658 | 50 | 50 |

## Descriptive statistics (sample-level)

| metric | count | mean | std | min | 25% | 50% | 75% | max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| amf_seq_richness | 50 | 13.16 | 26.85 | 1 | 2 | 4.5 | 9 | 156 |
| amf_genus_richness | 50 | 1.42 | 0.7309 | 1 | 1 | 1 | 2 | 4 |
| total_reads | 50 | 1757 | 5568 | 1 | 2 | 25 | 160.8 | 2.746e+04 |
| MAT | 50 | 16.29 | 9.888 | -5.9 | 10.35 | 19.05 | 24.58 | 27 |
| MAP | 50 | 1135 | 696.8 | 78 | 643 | 1022 | 1490 | 2921 |
| pH | 36 | 5.249 | 1.098 | 3.5 | 4.3 | 5.3 | 5.7 | 8 |

## GlobalAMFungi ↔ FRED matching

- global_unique_species_after_filters: **50**
- fred_unique_species_raw: **896**
- raw_exact_matches: **4**
- canonical_matches: **4**
- global_unique_genera_canonical: **48**
- fred_unique_genera_canonical: **506**
- genus_overlap_canonical: **12**

## Warnings / caveats

- None

## Exact text output printed in notebook

```text
Global input: ../Data/globalamf.csv
FRED input: ../Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv
Output dir: /home/stevens/projects/pal/Output
Global raw shape: (2171488, 27)
FRED raw (header=1, data-only) shape: (6271, 64)
initial: rows=2,171,488, unique_samples=5,524, unique_plant_species=441
root_only: rows=1,096,982, unique_samples=2,691, unique_plant_species=351
root_singleplant: rows=1,079,452, unique_samples=2,521, unique_plant_species=285
root_singleplant_singleton_species: rows=658, unique_samples=50, unique_plant_species=50
Sample-level table rows: 50
Wrote: ../Output/globalamfungi_sample_level.csv
Plant-species summary rows: 50
Wrote: ../Output/globalamfungi_species_level.csv
--- GlobalAMFungi ↔ FRED matching ---
global_unique_species_after_filters: 50
fred_unique_species_raw: 896
raw_exact_matches: 4
canonical_matches: 4
global_unique_genera_canonical: 48
fred_unique_genera_canonical: 506
genus_overlap_canonical: 12
Wrote figure: ../Output/globalamf_hist_amf_seq_richness.png
Wrote figure: ../Output/globalamf_scatter_richness_vs_MAT.png
Wrote figure: ../Output/globalamf_map_points.png
Wrote descriptive stats: ../Output/globalamf_descriptive_stats.csv
```