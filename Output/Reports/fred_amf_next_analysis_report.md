# FRED-AMF Next Analysis Report

## Dataset files used
- EcoBank: `Data/ecobank_full.csv` + `Data/ecobank_cultured.csv`
- GlobalAMFungi: `Output/globalamfungi_sample_level_full.csv`
- Primary FRED: `Data/FRED_4_full_20260312_2.csv`
- Sensitivity FRED: FRED3 + filtered FRED4 files

## Matching summary
- Overlap table: `/home/stevens/projects/pal/Output/Tables/fred_amf_species_overlap_summary.csv`
- Best overlap source is FRED4_full in both datasets.

## Accepted manual name corrections
- `deschampsia caespitosa` -> `deschampsia cespitosa`
- `hypochoeris radicata` -> `hypochaeris radicata`

## Number of species retained
- Species master rows: 6,695
- With EcoBank overlap: 500
- With GlobalAMFungi overlap: 275
- Genus master rows: 2,204

## Figure list
- Fig1_dataset_overlap_summary.(png,pdf)
- Fig2_geographic_coverage.(png,pdf)
- Fig3_root_trait_distributions.(png,pdf)
- Fig4_root_economics_pca.(png,pdf)
- Fig5_species_amf_richness_vs_root_traits.(png,pdf)
- Fig6_amf_guild_proportions_vs_root_traits.(png,pdf)
- Fig7_woody_herbaceous_sensitivity.(png,pdf)
- Fig8_genus_level_sensitivity.(png,pdf)

## Main preliminary statistical results
- Species-level models: `/home/stevens/projects/pal/Output/Tables/species_level_model_summaries.csv`
- Genus-level models: `/home/stevens/projects/pal/Output/Tables/genus_level_model_summaries.csv`

## Caveats
- FRED trait values and AMF observations are not necessarily co-located.
- Root N and Root P can be strongly environmentally variable.
- Root traits are partly phylogenetically conserved.
- Genus-level aggregation increases overlap but adds trait heterogeneity.
- Fuzzy matches are not automatically accepted.