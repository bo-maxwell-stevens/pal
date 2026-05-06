# AMF Traits Merged Summary

## Inputs used
- GlobalAMFungi sample-level: `../Output/globalamfungi_sample_level_full.csv`
- EcoBank source: derived from `../Data/ecobank_full.csv` + `../Data/ecobank_cultured.csv`
- FRED: `../Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv`

## Source counts before/after FRED merge

- GlobalAMFungi before merge rows: 2,521
- EcoBank before merge rows: 2,705
- GlobalAMFungi after merge rows: 2,521
- EcoBank after merge rows: 2,705

## Unique species per source

- GlobalAMFungi species: 285
- EcoBank species: 461

## Species matched to FRED

- GlobalAMFungi exact raw matches: 46
- GlobalAMFungi canonical matches: 46 / 285 (16.14%)
- EcoBank exact raw matches: 38
- EcoBank canonical matches: 38 / 461 (8.24%)

## Final merged dataset dimensions

- Rows: 5,226
- Columns: 70
- Unique samples (source+sample_id): 5,226
- Unique species (canonical): 568
- Fraction rows with >=1 non-null core trait: 0.221

## Missingness summary for key trait columns

| column | missing_count | missing_frac |
|---|---:|---:|
| Root N content | 4392 | 0.840 |
| Root P content | 4769 | 0.913 |
| Root diameter | 4476 | 0.856 |
| Root tissue density (RTD) | 4731 | 0.905 |
| Specific root area (SRA) | 4962 | 0.949 |
| Mycorrhiza_Fraction of root length or tips colonized | 5062 | 0.969 |
| Soil pH_Water | 4962 | 0.949 |
| Soil C content | 5202 | 0.995 |
| Soil N content | 5115 | 0.979 |
| Mean annual precipitation (MAP) | 4525 | 0.866 |
| Mean annual air temperature (MAT) | 4620 | 0.884 |

## Top 30 unmatched species

### GlobalAMFungi unmatched
- Zea mays: 126
- Daucus carota: 116
- Botrychium lunaria: 71
- Stachys sylvatica: 68
- Asclepias speciosa: 59
- Tanacetum vulgare: 48
- NA_: 48
- Solanum tuberosum: 42
- Centaurea stoebe: 36
- Solidago virgaurea: 34
- Alnus glutinosa: 32
- Eragostis tef: 27
- Voyria tenella: 27
- Nicotiana glauca: 24
- Artemisia vulgaris: 22
- Leymus mollis: 22
- Tamarix aphylla: 21
- Ammophila arenaria: 21
- Knautia arvensis: 20
- Populus deltoides: 20
- Themeda triandra: 18
- Digitaria macroblephara: 17
- Sorghum bicolor: 16
- Calamagrostis epigejos: 16
- Cenchrus setaceus: 15
- Talipariti tiliaceum: 13
- Hyparrhenia hirta: 13
- Pancratium maritimum: 12
- Symphiotrichum novi-belgii: 12
- Echinops sphaerocephalus: 12

### EcoBank unmatched
- Themeda triandra: 47
- Araucaria araucana: 40
- Festuca pallescens: 40
- Digitaria macroblephara: 40
- Hepatica nobilis: 40
- Potentilla chiloensis: 39
- Manihot esculenta: 36
- Galium verum: 32
- Lobelia dortmanna: 25
- Microsorum scolopendria: 24
- Paris quadrifolia: 22
- Filipendula vulgaris: 21
- Anemone multifida: 20
- Ranunculus polyanthemos: 20
- Helictotrichon pratense: 20
- Ipomoea nil: 19
- Nephrolepis hirsutula: 19
- Poa attenuata: 18
- Cestrum parqui: 18
- Sesleria caerulea: 18
- Hypochoeris radicata: 17
- Calamagrostis purpurea: 17
- Hibiscus tiliaceus: 17
- Solidago virgaurea: 16
- Leptorhynchos squamatus: 15
- Convallaria majalis: 14
- Viola mirabilis: 14
- Centaurea jacea: 14
- Veronica chamaedrys: 13
- Deschampsia flexuosa: 12

## Close-name diagnostics for unmatched species

- Method: `difflib.get_close_matches` with `SequenceMatcher` ratio cutoff `0.86` (genus-prioritized, then global fallback).
- Why close but not exact: most likely synonym updates, accepted-name changes between taxonomic references/versions, spelling variants, or rank/qualifier differences (e.g., `sp.`, `cf.`).

- GlobalAMFungi unmatched canonical species: 239
- EcoBank unmatched canonical species: 423
- GlobalAMFungi close-name suggestions: 1
- EcoBank close-name suggestions: 5

### GlobalAMFungi unmatched species with close FRED names

| unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---:|---:|
| filipendula vulgaris | filipendula ulmaria | 0.872 | 1 |

### EcoBank unmatched species with close FRED names

| unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---:|---:|
| deschampsia caespitosa | deschampsia cespitosa | 0.977 | 1 |
| hypochoeris radicata | hypochaeris radicata | 0.950 | 0 |
| cirsium arvense | cerastium arvense | 0.875 | 0 |
| artemisia rupestris | artemisia campestris | 0.872 | 1 |
| filipendula vulgaris | filipendula ulmaria | 0.872 | 1 |

## Warnings

- EcoBank total_reads missing; using amf_seq_richness as proxy total_reads.

## Printed output

```text
Using GlobalAMFungi file: ../Output/globalamfungi_sample_level_full.csv
Global rows: 2,521
Global unique samples: 2,521
Global unique species: 285
No EcoBank sample-level output file found; deriving from raw ecobank files as fallback.
WARNING: EcoBank total_reads missing; using amf_seq_richness as proxy total_reads.
EcoBank rows: 2,705
EcoBank unique samples: 2,705
EcoBank unique species: 461
GlobalAMFungi amf_seq_richness mean/median: 261.934/18.000
GlobalAMFungi total_reads mean/median: 4178282.573/761.000
EcoBank amf_seq_richness mean/median: 18.015/16.000
EcoBank total_reads mean/median: 18.021/16.000
FRED rows before aggregation: 5,720
FRED species before aggregation: 896
FRED species after aggregation: 896
Number of FRED trait columns carried: 56
GlobalAMFungi exact raw species matches: 46
GlobalAMFungi canonical species matches: 46 / 285 (16.14%)
EcoBank exact raw species matches: 38
EcoBank canonical species matches: 38 / 461 (8.24%)
Global rows after FRED merge: 2,521
EcoBank rows after FRED merge: 2,705
--- Final merged dataset QA ---
total rows: 5,226
unique samples overall: 5,226
unique species overall: 568
EcoBank rows: 2,705, unique_samples: 2,705, unique_species: 461
GlobalAMFungi rows: 2,521, unique_samples: 2,521, unique_species: 285
fraction rows with >=1 non-null core trait: 0.221
Wrote merged dataset: ../Output/amf_traits_merged_sample_level.csv
GlobalAMFungi unmatched canonical species: 239
EcoBank unmatched canonical species: 423
GlobalAMFungi close-name suggestions (cutoff=0.86): 1
EcoBank close-name suggestions (cutoff=0.86): 5
Wrote close-name diagnostics: ../Output/amf_traits_close_name_diagnostics.md
```