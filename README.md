# PAL Analysis Pipeline

This repository now includes script-based workflows so analyses can be rerun repeatedly without re-downloading source data.

## What is script-converted

- `scripts/01_prepare_global_full.py`
  - Converts the core logic from `Code/01_data_prep_global_full.ipynb` into a script.
- `scripts/compare_fred_sources.py`
  - Reusable species matching comparison across old and new FRED files.
- `scripts/run_match_analysis.sh`
  - Convenience wrapper for repeated runs.

The notebooks in `Code/` remain as reference/provenance.

## Run matching comparison (recommended first step)

```bash
PYTHON_BIN=/home/stevens/miniconda3/bin/python scripts/run_match_analysis.sh
```

Outputs:
- `Output/fred_species_match_summary.csv`
- `Output/fred_species_likely_matches.csv`
- `Output/fred_species_match_report.md`

## Run full FRED4-primary analysis

```bash
PYTHONPATH=src /home/stevens/miniconda3/bin/python scripts/run_fred4_amf_analysis.py
```

Primary outputs:
- `Output/Tables/fred_amf_species_overlap_summary.csv`
- `Output/Tables/fred4_amf_species_master_table.csv`
- `Output/Tables/fred4_amf_genus_master_table.csv`
- `Output/Tables/species_level_model_summaries.csv`
- `Output/Tables/genus_level_model_summaries.csv`
- `Output/Figures/Fig1_...` through `Fig8_...` in PNG and PDF
- `Output/Reports/fred_amf_next_analysis_report.md`

## Notes

- No download step is required for this comparison.
- The script reads local files in `Data/` and `Output/` only.
- FRED parsing is robust to encodings (`utf-8`, `latin-1`, `cp1252`) and uses the established header-row convention from your notebooks.
