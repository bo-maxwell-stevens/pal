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

## Notebook conversion

To regenerate Python script versions of all notebooks in `Code/`:

```bash
/home/stevens/miniconda3/bin/python scripts/convert_notebooks_to_py.py
```

Generated files are written to `scripts/notebook_converted/`.
