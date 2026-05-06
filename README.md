# PAL: FRED x AMF Analysis

This repository contains a reproducible analysis pipeline linking plant root traits (FRED) to AMF community metrics (EcoBank and GlobalAMFungi).

Primary objective:
- Use **full FRED4** (`Data/FRED_4_full_20260312_2.csv`) as the main trait source.
- Use FRED3 and filtered FRED4 files as sensitivity sources.

## Start Here

- Project overview for collaborators: `docs/COLLABORATOR_GUIDE.md`
- Reproducible workflow: `docs/WORKFLOW.md`
- Species matching summary: `docs/RESULTS_SPECIES_MATCH.md`
- Final analysis report: `Output/Reports/fred_amf_next_analysis_report.md`
- Figure trends interpretation: `Output/Reports/fred_amf_figure_trends_summary.md`

## Repository Layout

- `Data/`: input datasets (FRED, EcoBank, GlobalAMFungi, reference files)
- `src/pal/`: reusable pipeline and plotting modules
- `scripts/`: runnable analysis scripts
- `Output/Tables/`: final tables
- `Output/Figures/`: final figures (`Fig1`-`Fig8`, PNG)
- `Output/Reports/`: markdown reports
- `Code/`: original notebooks kept for provenance

## Main Commands

Run full analysis pipeline (recommended):

```bash
PYTHONPATH=src /home/stevens/miniconda3/bin/python scripts/run_fred4_amf_analysis.py
```

Run species-match comparison only:

```bash
PYTHONPATH=src /home/stevens/miniconda3/bin/python scripts/compare_fred_sources.py
```

## Key Outputs

Tables:
- `Output/Tables/fred_amf_species_overlap_summary.csv`
- `Output/Tables/fred4_amf_species_master_table.csv`
- `Output/Tables/fred4_amf_genus_master_table.csv`
- `Output/Tables/species_level_model_summaries.csv`
- `Output/Tables/genus_level_model_summaries.csv`

Figures:
- `Output/Figures/Fig1_dataset_overlap_summary.png`
- `Output/Figures/Fig2_geographic_coverage.png`
- `Output/Figures/Fig3_root_trait_distributions.png`
- `Output/Figures/Fig4_root_economics_pca.png`
- `Output/Figures/Fig5_species_amf_richness_vs_root_traits.png`
- `Output/Figures/Fig6_amf_guild_proportions_vs_root_traits.png`
- `Output/Figures/Fig7_woody_herbaceous_sensitivity.png`
- `Output/Figures/Fig8_genus_level_sensitivity.png`

Reports:
- `Output/Reports/fred_amf_next_analysis_report.md`
- `Output/Reports/fred_amf_figure_trends_summary.md`

## Current Interpretation Status

- Robust: overlap summaries and source-comparison conclusions (`Fig1`, overlap tables).
- Exploratory: species-level multivariate trend panels (`Fig5`-`Fig7`) due to limited complete-case n.
- Caution: rhizophilic proportion models can be numerically fragile because values are often near 100%.
