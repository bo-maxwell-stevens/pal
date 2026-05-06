# Workflow

## Goal

Run species-match analysis repeatedly without re-downloading source data.

## One-time setup

- Ensure Python environment has `pandas` (current tested binary: `/home/stevens/miniconda3/bin/python`).

## Repeated run

```bash
PYTHONPATH=src /home/stevens/miniconda3/bin/python scripts/compare_fred_sources.py
```

Outputs are regenerated in `Output/`:
- `fred_species_match_summary.csv`
- `fred_species_likely_matches.csv`
- `fred_species_match_report.md`

## Full FRED4-primary analysis run

```bash
PYTHONPATH=src /home/stevens/miniconda3/bin/python scripts/run_fred4_amf_analysis.py
```

This generates:
- overlap summary table,
- species-level and genus-level master tables,
- species/genus model summary tables,
- multi-panel figures `Fig1` to `Fig8` in PNG and PDF,
- report: `Output/Reports/fred_amf_next_analysis_report.md`.

## Notebook conversion

To regenerate Python script versions of all notebooks in `Code/`:

```bash
/home/stevens/miniconda3/bin/python scripts/convert_notebooks_to_py.py
```

Generated files are written to `scripts/notebook_converted/`.
