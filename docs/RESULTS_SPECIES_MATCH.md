# Species Matching Results

Run date: 2026-05-06

This document compares species matching between plant datasets and available FRED sources, including the new full file `Data/FRED_4_full_20260312_2.csv`.

## Exact canonical match counts

| dataset | fred_source | dataset_species | fred_species | exact_matches | match_percent |
|---|---|---:|---:|---:|---:|
| EcoBank | FRED_4_full_20260312_2.csv | 529 | 7005 | 153 | 28.92 |
| EcoBank | FRED3_fineRoots.csv | 529 | 204 | 83 | 15.69 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | 529 | 992 | 37 | 6.99 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv | 529 | 558 | 11 | 2.08 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | 285 | 7005 | 112 | 39.30 |
| GlobalAMFungi | FRED3_fineRoots.csv | 285 | 204 | 65 | 22.81 |
| GlobalAMFungi | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | 285 | 992 | 45 | 15.79 |
| GlobalAMFungi | FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv | 285 | 558 | 18 | 6.32 |

## Interpretation

- The new full FRED has the largest species coverage and gives the highest exact match counts for both datasets.
- Relative gain is substantial vs the previously used filtered files.
- `FRED3_fineRoots.csv` still performs well despite being much smaller, likely because of overlap in sampled taxa.

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
