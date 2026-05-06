# Collaborator Guide

This guide is intended for collaborators reviewing this repository on GitHub.

## What this project does

We integrate:
- plant root trait data from FRED,
- AMF metrics from EcoBank and GlobalAMFungi,

to evaluate trait-richness and trait-guild relationships.

Primary FRED source is full FRED4:
- `Data/FRED_4_full_20260312_2.csv`

Sensitivity FRED sources:
- `Data/FRED3_fineRoots.csv`
- `Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv`
- `Data/FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv`

## Recommended reading order

1. `README.md`
2. `docs/RESULTS_SPECIES_MATCH.md`
3. `Output/Reports/fred_amf_next_analysis_report.md`
4. `Output/Reports/fred_amf_figure_trends_summary.md`
5. `Output/Figures/Fig1_dataset_overlap_summary.png` through `Fig8_genus_level_sensitivity.png`

## Reproduce results

```bash
PYTHONPATH=src /home/stevens/miniconda3/bin/python scripts/run_fred4_amf_analysis.py
```

## How to understand figure roles

- `Fig1`: source overlap and match performance.
- `Fig2`: geographic context.
- `Fig3`: trait coverage/distributions.
- `Fig4`: root economics PCA structure.
- `Fig5`: species-level richness-trait trends (exploratory).
- `Fig6`: species-level guild-trait trends (exploratory; rhizophilic often near upper bound).
- `Fig7`: woody vs non-woody sensitivity check.
- `Fig8`: genus-level robustness check.

## Known limitations (important)

- Trait and AMF values are not necessarily co-located in space/time.
- Complete-case overlap for species-level multivariate models is limited.
- Root N and Root P are environmentally plastic.
- Root traits are phylogenetically structured; current models are non-phylogenetic.
- Fuzzy matches are not auto-accepted.
