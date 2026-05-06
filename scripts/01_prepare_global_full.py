#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import pandas as pd


def n_plants(value: object) -> int:
    if value is None:
        return 0
    s = str(value).strip()
    if not s or s.lower() in {"nan", "none"}:
        return 0
    return len([x for x in s.split(";") if x.strip()])


def main() -> None:
    data_fp = Path("Data/globalamf.csv")
    out_dir = Path("Output")
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(data_fp, low_memory=False)
    df["sample_type_norm"] = df["sample_type"].astype(str).str.strip().str.lower()
    root = df[df["sample_type_norm"] == "root"].copy()
    root["n_plants"] = root["plants_dominant"].map(n_plants)
    one = root[root["n_plants"] == 1].copy()
    one["plant_species"] = one["plants_dominant"].astype(str).str.strip()
    one = one[one["plant_species"].notna() & (one["plant_species"] != "")].copy()

    sample = (
        one.groupby(["id", "plant_species"], as_index=False)
        .agg(
            latitude=("latitude", "first"),
            longitude=("longitude", "first"),
            MAT=("MAT", "first"),
            MAP=("MAP", "first"),
            pH=("pH", "first"),
            amf_seq_richness=("sequence", "nunique"),
            amf_genus_richness=("Genus", lambda s: pd.Series(s).dropna().astype(str).nunique()),
            total_reads=("abundances", lambda s: pd.to_numeric(s, errors="coerce").fillna(0).sum()),
        )
        .rename(columns={"id": "sample_id"})
    )

    species = (
        sample.groupby("plant_species", as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            mean_amf_seq_richness=("amf_seq_richness", "mean"),
            mean_amf_genus_richness=("amf_genus_richness", "mean"),
            mean_total_reads=("total_reads", "mean"),
        )
        .sort_values("n_samples", ascending=False)
    )

    sample_fp = out_dir / "globalamfungi_sample_level_full.csv"
    species_fp = out_dir / "globalamfungi_species_level_full.csv"
    summary_fp = out_dir / "globalamfungi_data_prep_summary_full.md"

    sample.to_csv(sample_fp, index=False)
    species.to_csv(species_fp, index=False)

    lines = [
        "# GlobalAMFungi Data Prep Summary (Script)",
        "",
        f"- Input: `{data_fp}`",
        f"- Output sample-level: `{sample_fp}`",
        f"- Output species-level: `{species_fp}`",
        "",
        "## Counts",
        f"- initial rows: {len(df):,}",
        f"- root rows: {len(root):,}",
        f"- root + single-dominant-plant rows: {len(one):,}",
        f"- sample-level rows: {len(sample):,}",
        f"- species-level rows: {len(species):,}",
    ]
    summary_fp.write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote: {sample_fp}")
    print(f"Wrote: {species_fp}")
    print(f"Wrote: {summary_fp}")


if __name__ == "__main__":
    main()
