# Species Matching Results

Run date: 2026-05-06

This document compares species matching between plant datasets and available FRED sources, including the new full file `Data/FRED_4_full_20260312_2.csv`.

For broader context and downstream interpretation, see:
- `docs/COLLABORATOR_GUIDE.md`
- `Output/Reports/fred_amf_next_analysis_report.md`

## Exact canonical match counts

| dataset | fred_source | dataset_species | fred_species | exact_matches | match_percent |
|---|---|---:|---:|---:|---:|
| EcoBank | FRED4_full | 500 | 6307 | 160 | 32.00 |
| EcoBank | FRED3 | 500 | 176 | 89 | 17.80 |
| EcoBank | FRED4_filtered_fineroot_lt2mm | 500 | 896 | 41 | 8.20 |
| EcoBank | FRED4_filtered_1storder | 500 | 535 | 13 | 2.60 |
| GlobalAMFungi | FRED4_full | 275 | 6307 | 121 | 44.00 |
| GlobalAMFungi | FRED3 | 275 | 176 | 68 | 24.73 |
| GlobalAMFungi | FRED4_filtered_fineroot_lt2mm | 275 | 896 | 46 | 16.73 |
| GlobalAMFungi | FRED4_filtered_1storder | 275 | 535 | 19 | 6.91 |

## Interpretation

- The new full FRED has the largest species coverage and gives the highest exact match counts for both datasets.
- Relative gain is substantial vs the previously used filtered files.
- `FRED3_fineRoots.csv` still performs well despite being much smaller, likely because of overlap in sampled taxa.
- The current pipeline uses strict canonical binomials; counts may differ from earlier runs that retained non-binomial names.

## Likely (non-exact) matches

Likely name matches are produced with `difflib.get_close_matches`, genus-prioritized, cutoff `0.86`.

See:
- `Output/fred_species_likely_matches.csv`
- `Output/fred_species_match_report.md`

## Reproducibility

Command used:

```bash
PYTHONPATH=src /home/stevens/miniconda3/bin/python scripts/compare_fred_sources.py
```

Related interpretation report:
- `Output/Reports/fred_amf_figure_trends_summary.md`
