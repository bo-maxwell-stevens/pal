#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from pal.fred_amf_figures import (
    fig1_overlap,
    fig2_geo,
    fig3_trait_distributions,
    fig4_pca,
    fig5_richness_vs_traits,
    fig6_guilds_vs_traits,
    fig7_woody_sensitivity,
    fig8_genus_sensitivity,
)
from pal.fred_amf_pipeline import (
    FRED_SOURCES,
    MANUAL_CORRECTIONS,
    add_pca,
    build_fred4_species_master,
    build_genus_master,
    ensure_dirs,
    load_ecobank,
    load_fred_species,
    load_globalamf,
    overlap_summary,
    run_models_genus,
    run_models_species,
)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    paths = ensure_dirs(root)

    eco_sample, eco_species = load_ecobank(root)
    glo_sample, glo_species = load_globalamf(root)

    fred_map = {k: load_fred_species(root, root / v) for k, v in FRED_SOURCES.items()}

    overlap = overlap_summary(
        set(eco_species["canonical_species"].dropna()),
        set(glo_species["canonical_species"].dropna()),
        fred_map,
    )
    overlap_fp = paths["tables"] / "fred_amf_species_overlap_summary.csv"
    overlap.to_csv(overlap_fp, index=False)

    fred4 = fred_map["FRED4_full"]
    species_master = build_fred4_species_master(fred4, eco_species, glo_species)
    species_master = add_pca(species_master)
    species_fp = paths["tables"] / "fred4_amf_species_master_table.csv"
    species_master.to_csv(species_fp, index=False)

    genus_master = build_genus_master(species_master)
    genus_fp = paths["tables"] / "fred4_amf_genus_master_table.csv"
    genus_master.to_csv(genus_fp, index=False)

    species_models = run_models_species(species_master)
    species_models_fp = paths["tables"] / "species_level_model_summaries.csv"
    species_models.to_csv(species_models_fp, index=False)

    genus_models = run_models_genus(genus_master)
    genus_models_fp = paths["tables"] / "genus_level_model_summaries.csv"
    genus_models.to_csv(genus_models_fp, index=False)

    fig1_overlap(overlap, paths["figures"] / "Fig1_dataset_overlap_summary")
    fig2_geo(fred4, eco_sample, glo_sample, paths["figures"] / "Fig2_geographic_coverage")
    fig3_trait_distributions(species_master, paths["figures"] / "Fig3_root_trait_distributions")
    fig4_pca(species_master, paths["figures"] / "Fig4_root_economics_pca")
    fig5_richness_vs_traits(species_master, paths["figures"] / "Fig5_species_amf_richness_vs_root_traits")
    fig6_guilds_vs_traits(species_master, paths["figures"] / "Fig6_amf_guild_proportions_vs_root_traits")
    fig7_woody_sensitivity(species_master, paths["figures"] / "Fig7_woody_herbaceous_sensitivity")
    fig8_genus_sensitivity(genus_master, paths["figures"] / "Fig8_genus_level_sensitivity")

    report = paths["reports"] / "fred_amf_next_analysis_report.md"
    lines = [
        "# FRED-AMF Next Analysis Report",
        "",
        "## Dataset files used",
        f"- EcoBank: `Data/ecobank_full.csv` + `Data/ecobank_cultured.csv`",
        f"- GlobalAMFungi: `Output/globalamfungi_sample_level_full.csv`",
        f"- Primary FRED: `{FRED_SOURCES['FRED4_full']}`",
        "- Sensitivity FRED: FRED3 + filtered FRED4 files",
        "",
        "## Matching summary",
        f"- Overlap table: `{overlap_fp}`",
        f"- Best overlap source is FRED4_full in both datasets.",
        "",
        "## Accepted manual name corrections",
    ]
    for k, v in MANUAL_CORRECTIONS.items():
        lines.append(f"- `{k}` -> `{v}`")

    lines += [
        "",
        "## Number of species retained",
        f"- Species master rows: {len(species_master):,}",
        f"- With EcoBank overlap: {int((species_master['in_EcoBank'] == 1).sum()):,}",
        f"- With GlobalAMFungi overlap: {int((species_master['in_GlobalAMFungi'] == 1).sum()):,}",
        f"- Genus master rows: {len(genus_master):,}",
        "",
        "## Figure list",
        "- Fig1_dataset_overlap_summary.png",
        "- Fig2_geographic_coverage.png",
        "- Fig3_root_trait_distributions.png",
        "- Fig4_root_economics_pca.png",
        "- Fig5_species_amf_richness_vs_root_traits.png",
        "- Fig6_amf_guild_proportions_vs_root_traits.png",
        "- Fig7_woody_herbaceous_sensitivity.png",
        "- Fig8_genus_level_sensitivity.png",
        "",
        "## Main preliminary statistical results",
        f"- Species-level models: `{species_models_fp}`",
        f"- Genus-level models: `{genus_models_fp}`",
        "",
        "## Caveats",
        "- FRED trait values and AMF observations are not necessarily co-located.",
        "- Root N and Root P can be strongly environmentally variable.",
        "- Root traits are partly phylogenetically conserved.",
        "- Genus-level aggregation increases overlap but adds trait heterogeneity.",
        "- Fuzzy matches are not automatically accepted.",
    ]
    report.write_text("\n".join(lines), encoding="utf-8")

    print("Generated outputs:")
    for fp in [
        overlap_fp,
        species_fp,
        genus_fp,
        species_models_fp,
        genus_models_fp,
        report,
    ]:
        print(f"- {fp}")


if __name__ == "__main__":
    main()
