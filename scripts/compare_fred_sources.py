#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pal.loaders import first_present, read_fred_table
from pal.species_match import compute_match_stats, likely_matches
from pal.text import canonicalize_species


def species_from_global(path: Path) -> set[str]:
    df = pd.read_csv(path)
    col = first_present(df.columns, ["plant_species", "species"])
    if not col:
        raise ValueError(f"No species column in {path}")
    canon = df[col].map(canonicalize_species)
    return {x for x in canon.dropna().unique().tolist() if x}


def species_from_ecobank(path: Path) -> set[str]:
    df = pd.read_csv(path, sep=";", low_memory=False)
    col = first_present(df.columns, ["plant_species", "species"])
    if not col:
        raise ValueError(f"No species column in {path}")
    canon = df[col].map(canonicalize_species)
    return {x for x in canon.dropna().unique().tolist() if x}


def species_from_fred(path: Path) -> set[str]:
    df = read_fred_table(path)
    species_col = first_present(
        df.columns,
        [
            "Plant taxonomy_Species_Data source",
            "Plant taxonomy_Species name unplaced",
            "Name",
        ],
    )
    genus_col = first_present(
        df.columns,
        [
            "Plant taxonomy_Accepted genus_WFO",
            "Plant taxonomy_Genus_Data Source",
            "Genus",
        ],
    )
    if not species_col:
        raise ValueError(f"No recognized FRED species column in {path}")

    species_raw = df[species_col].astype(str)
    if genus_col and genus_col in df.columns:
        genus_raw = df[genus_col].astype(str)
        full_name = []
        for g, s in zip(genus_raw, species_raw):
            gs = str(g).strip()
            ss = str(s).strip()
            if " " in ss:
                full_name.append(ss)
            elif gs and ss and ss.lower() not in {"nan", "none"}:
                full_name.append(f"{gs} {ss}")
            else:
                full_name.append(ss)
        canon = pd.Series(full_name).map(canonicalize_species)
    else:
        canon = species_raw.map(canonicalize_species)

    return {x for x in canon.dropna().unique().tolist() if x}


def markdown_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        vals = [str(r[c]) for c in cols]
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare species matching across FRED source files.")
    parser.add_argument("--global-sample", default="Output/globalamfungi_sample_level_full.csv")
    parser.add_argument("--ecobank", default="Data/ecobank_full.csv")
    parser.add_argument(
        "--fred",
        nargs="+",
        default=[
            "Data/FRED3_fineRoots.csv",
            "Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv",
            "Data/FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv",
            "Data/FRED_4_full_20260312_2.csv",
        ],
    )
    parser.add_argument("--outdir", default="Output")
    parser.add_argument("--cutoff", type=float, default=0.86)
    parser.add_argument("--top-n", type=int, default=30)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    global_species = species_from_global(Path(args.global_sample))
    ecobank_species = species_from_ecobank(Path(args.ecobank))

    match_rows = []
    likely_rows = []
    for fred in [Path(x) for x in args.fred]:
        fred_species = species_from_fred(fred)
        for dataset_name, species_set in [("GlobalAMFungi", global_species), ("EcoBank", ecobank_species)]:
            stat = compute_match_stats(dataset_name, fred.name, species_set, fred_species)
            row = {
                "dataset": stat.dataset,
                "fred_source": stat.fred_source,
                "n_dataset_species": stat.n_dataset_species,
                "n_fred_species": stat.n_fred_species,
                "exact_matches": stat.exact_matches,
                "match_percent": round(100.0 * stat.exact_matches / max(stat.n_dataset_species, 1), 2),
                "genus_overlap": stat.genus_overlap,
            }
            match_rows.append(row)

            unmatched = species_set - fred_species
            suggestions = likely_matches(unmatched, fred_species, cutoff=args.cutoff, top_n=args.top_n)
            for s in suggestions:
                likely_rows.append(
                    {
                        "dataset": dataset_name,
                        "fred_source": fred.name,
                        **s,
                    }
                )

    summary = pd.DataFrame(match_rows).sort_values(["dataset", "exact_matches"], ascending=[True, False])
    likely_df = pd.DataFrame(likely_rows)

    summary_csv = outdir / "fred_species_match_summary.csv"
    likely_csv = outdir / "fred_species_likely_matches.csv"
    summary.to_csv(summary_csv, index=False)
    likely_df.to_csv(likely_csv, index=False)

    lines = [
        "# FRED Species Match Comparison",
        "",
        "This report compares exact canonical species matching from each dataset to multiple FRED source files.",
        "",
        "## Match summary",
        "",
        markdown_table(summary),
        "",
        "## Likely name matches",
        "",
        "Likely matches are from `difflib.get_close_matches` with genus-prioritized matching and cutoff "
        f"`{args.cutoff}`.",
        "",
    ]

    for (dataset, fred), block in likely_df.groupby(["dataset", "fred_source"], sort=False):
        lines.append(f"### {dataset} vs {fred}")
        lines.append("")
        lines.append(markdown_table(block.head(args.top_n)))
        lines.append("")

    report_md = outdir / "fred_species_match_report.md"
    report_md.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote: {summary_csv}")
    print(f"Wrote: {likely_csv}")
    print(f"Wrote: {report_md}")


if __name__ == "__main__":
    main()
