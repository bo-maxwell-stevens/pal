# FRED Species Match Comparison

This report compares exact canonical species matching from each dataset to multiple FRED source files.

## Match summary

| dataset | fred_source | n_dataset_species | n_fred_species | exact_matches | match_percent | genus_overlap |
|---|---|---|---|---|---|---|
| EcoBank | FRED_4_full_20260312_2.csv | 529 | 7005 | 153 | 28.92 | 232 |
| EcoBank | FRED3_fineRoots.csv | 529 | 204 | 83 | 15.69 | 88 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | 529 | 992 | 37 | 6.99 | 115 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv | 529 | 558 | 11 | 2.08 | 50 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | 285 | 7005 | 112 | 39.3 | 171 |
| GlobalAMFungi | FRED3_fineRoots.csv | 285 | 204 | 65 | 22.81 | 76 |
| GlobalAMFungi | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | 285 | 992 | 45 | 15.79 | 97 |
| GlobalAMFungi | FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv | 285 | 558 | 18 | 6.32 | 46 |

## Likely name matches

Likely matches are from `difflib.get_close_matches` with genus-prioritized matching and cutoff `0.86`.

### GlobalAMFungi vs FRED3_fineRoots.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| GlobalAMFungi | FRED3_fineRoots.csv | dichrostachys cinerea | dichronostachys cinerea | 0.955 | 0 |

### EcoBank vs FRED3_fineRoots.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| EcoBank | FRED3_fineRoots.csv | dichrostachys cinerea | dichronostachys cinerea | 0.955 | 0 |
| EcoBank | FRED3_fineRoots.csv | epilobium latifolium | epilobium angustifolium | 0.884 | 1 |

### GlobalAMFungi vs FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| GlobalAMFungi | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | filipendula vulgaris | filipendula ulmaria | 0.872 | 1 |

### EcoBank vs FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | deschampsia caespitosa | deschampsia cespitosa | 0.977 | 1 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | hypochoeris radicata | hypochaeris radicata | 0.95 | 0 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | epilobium latifolium | epilobium angustifolium | 0.884 | 1 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | nepeta sp | neea sp | 0.875 | 0 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | artemisia rupestris | artemisia campestris | 0.872 | 1 |
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv | filipendula vulgaris | filipendula ulmaria | 0.872 | 1 |

### GlobalAMFungi vs FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| GlobalAMFungi | FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv | poncirus trifoliata | citrus trifoliata | 0.889 | 0 |
| GlobalAMFungi | FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv | dryopteris sp | dryopteris sparsa | 0.867 | 1 |

### EcoBank vs FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| EcoBank | FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv | dryopteris sp | dryopteris sparsa | 0.867 | 1 |

### GlobalAMFungi vs FRED_4_full_20260312_2.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| GlobalAMFungi | FRED_4_full_20260312_2.csv | calamagrostis epigejos | calamagrostis epigeios | 0.955 | 1 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | crossopterix febrifuga | crossopteryx febrifuga | 0.955 | 0 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | loropetalum chinense | loropetalum chenense | 0.95 | 1 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | megathyrsus maximus | megathyrsus maximum | 0.947 | 1 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | cenchrus setaceus | cenchrus setaceum | 0.941 | 1 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | helictochloa pratensis | helictochloa pratense | 0.93 | 1 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | plantago afra | plantago varia | 0.889 | 1 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | poncirus trifoliata | citrus trifoliata | 0.889 | 0 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | chimonanthus salicifolius | helianthus salicifolius | 0.875 | 0 |
| GlobalAMFungi | FRED_4_full_20260312_2.csv | dryopteris sp | dryopteris sparsa | 0.867 | 1 |

### EcoBank vs FRED_4_full_20260312_2.csv

| dataset | fred_source | unmatched_name | candidate_fred_name | similarity_ratio | same_genus |
|---|---|---|---|---|---|
| EcoBank | FRED_4_full_20260312_2.csv | hypochoeris radicata | hypochaeris radicata | 0.95 | 0 |
| EcoBank | FRED_4_full_20260312_2.csv | digitaria sanguinalis | digitaria sanguinale | 0.927 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | plantago afra | plantago varia | 0.889 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | epilobium latifolium | epilobium angustifolium | 0.884 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | artemisia persica | artemisia sericea | 0.882 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | artemisia macrocephala | artemisia sphaerocephala | 0.87 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | poa tibetica | poa iberica | 0.87 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | dryopteris sp | dryopteris sparsa | 0.867 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | potentilla sp | potentilla supina | 0.867 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | allium oreoprasum | allium scorodoprasum | 0.865 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | potentilla curviseta | potentilla grisea | 0.865 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | chaerophyllum villosum | chaerophyllum bulbosum | 0.864 | 1 |
| EcoBank | FRED_4_full_20260312_2.csv | pedicularis cheilanthifolia | pedicularis anthemifolia | 0.863 | 1 |
