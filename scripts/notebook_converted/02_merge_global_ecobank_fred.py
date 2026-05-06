# Auto-converted from 02_merge_global_ecobank_fred.ipynb


# %% [cell 1] type=markdown
# # 02_merge_global_ecobank_fred
# 
# Reproducible merge pipeline for modeling-ready sample-level data combining:
# 1) GlobalAMFungi
# 2) EcoBank
# 3) FRED 4.0 traits


# %% [cell 2] type=markdown
# ## Block 1: Imports, paths, logging, and helpers


# %% [cell 3] type=code
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path('.')
DATA_DIR = ROOT / '../Data'
OUT_DIR = ROOT / '../Output'
OUT_DIR.mkdir(parents=True, exist_ok=True)

FRED_FP = DATA_DIR / 'FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv'
assert FRED_FP.exists(), f'Missing {FRED_FP}'

MERGED_OUT = OUT_DIR / 'amf_traits_merged_sample_level.csv'
SUMMARY_OUT = OUT_DIR / 'amf_traits_merged_summary.md'

LOG = []
WARN = []

def emit(msg):
    s = str(msg)
    print(s)
    LOG.append(s)

def warn(msg):
    s = f'WARNING: {msg}'
    print(s)
    LOG.append(s)
    WARN.append(str(msg))

def canonicalize(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    s = re.sub(r'\s+', ' ', s)
    s = re.sub(r'[\.;:,]+$', '', s)
    return s

def clean_spaces(x):
    if pd.isna(x):
        return np.nan
    return re.sub(r'\s+', ' ', str(x).strip())

def first_present(df, names):
    for n in names:
        if n in df.columns:
            return n
    return None


# %% [cell 4] type=markdown
# ## Block 2: Load and standardize GlobalAMFungi sample-level


# %% [cell 5] type=code
global_candidates = [
    OUT_DIR / 'globalamfungi_sample_level_full.csv',
    OUT_DIR / 'globalamfungi_sample_level.csv',
]
global_fp = next((p for p in global_candidates if p.exists()), None)
if global_fp is None:
    raise FileNotFoundError('No GlobalAMFungi sample-level file found in expected locations.')

g = pd.read_csv(global_fp)
emit(f'Using GlobalAMFungi file: {global_fp}')

sample_col = first_present(g, ['sample_id', 'id', 'sample'])
species_col = first_present(g, ['plant_species', 'species'])
if sample_col is None or species_col is None:
    raise ValueError('GlobalAMFungi sample-level must contain sample id and plant species columns.')

g = g.rename(columns={sample_col: 'sample_id', species_col: 'plant_species'})

# study_id preference, with fallback to raw global file paper_id by id
study_col = first_present(g, ['study_id', 'paper_id', 'dataset_id'])
if study_col is not None:
    g['study_id'] = g[study_col].astype(str)
else:
    raw_global_fp = DATA_DIR / 'globalamf.csv'
    if raw_global_fp.exists():
        raw_g = pd.read_csv(raw_global_fp, usecols=['id', 'paper_id'])
        raw_g = raw_g.drop_duplicates(subset=['id'])
        g = g.merge(raw_g, left_on='sample_id', right_on='id', how='left')
        g['study_id'] = g['paper_id'].fillna('unknown').map(lambda x: f'GlobalAMFungi_{x}')
        g = g.drop(columns=['id'])
    else:
        g['study_id'] = 'GlobalAMFungi_unknown'

for col in ['MAT', 'MAP', 'pH', 'latitude', 'longitude', 'amf_seq_richness', 'amf_genus_richness', 'total_reads']:
    if col not in g.columns:
        g[col] = np.nan

g['source'] = 'GlobalAMFungi'
g['plant_species'] = g['plant_species'].map(clean_spaces)
g['plant_species_canon'] = g['plant_species'].map(canonicalize)

for col in ['MAT', 'MAP', 'pH', 'latitude', 'longitude', 'amf_seq_richness', 'amf_genus_richness', 'total_reads']:
    g[col] = pd.to_numeric(g[col], errors='coerce')

g['total_reads'] = g['total_reads'].where(g['total_reads'] > 0, np.nan)
g['log_total_reads'] = np.log(g['total_reads'])

emit(f'Global rows: {len(g):,}')
emit(f"Global unique samples: {g['sample_id'].nunique():,}")
emit(f"Global unique species: {g['plant_species_canon'].dropna().nunique():,}")


# %% [cell 6] type=markdown
# ## Block 3: Load/derive and standardize EcoBank sample-level


# %% [cell 7] type=code
eco_candidates = [
    OUT_DIR / 'ecobank_sample_level.csv',
    OUT_DIR / 'ecobank_sample_level_full.csv',
    OUT_DIR / 'ecobank_samples.csv',
    OUT_DIR / 'sample_level_ecobank.csv',
]
eco_fp = next((p for p in eco_candidates if p.exists()), None)

if eco_fp is not None:
    e = pd.read_csv(eco_fp)
    emit(f'Using EcoBank sample-level file: {eco_fp}')
else:
    emit('No EcoBank sample-level output file found; deriving from raw ecobank files as fallback.')
    eco_raw_fp = DATA_DIR / 'ecobank_full.csv'
    vt_fp = DATA_DIR / 'ecobank_cultured.csv'
    assert eco_raw_fp.exists(), eco_raw_fp
    assert vt_fp.exists(), vt_fp

    eco_raw = pd.read_csv(eco_raw_fp, sep=';')
    vt = pd.read_csv(vt_fp, sep=';')
    eco_raw = eco_raw.sort_values('sample')
    vt = vt.sort_values('sample')
    roots = eco_raw.merge(vt, on='sample', how='inner')

    if roots['isRoot'].dtype == object:
        is_root = roots['isRoot'].astype(str).str.upper().isin(['TRUE', 'T', '1', 'YES'])
    else:
        is_root = roots['isRoot'].astype(bool)
    roots = roots.loc[is_root].copy()

    # Sample-level already one row per sample in this merged table
    e = roots.copy()

sample_col = first_present(e, ['sample_id', 'id', 'sample'])
species_col = first_present(e, ['plant_species', 'species'])
if sample_col is None or species_col is None:
    raise ValueError('EcoBank data needs sample and species columns.')

e = e.rename(columns={sample_col: 'sample_id', species_col: 'plant_species'})

# richness/reads fallback logic
if 'amf_seq_richness' not in e.columns:
    if 'richness' in e.columns:
        e['amf_seq_richness'] = pd.to_numeric(e['richness'], errors='coerce')
    else:
        e['amf_seq_richness'] = np.nan
        warn('EcoBank has no richness column; amf_seq_richness set to NaN.')

if 'amf_genus_richness' not in e.columns:
    e['amf_genus_richness'] = np.nan

if 'total_reads' not in e.columns:
    # proxy because EcoBank table lacks read counts in this fallback
    e['total_reads'] = pd.to_numeric(e['amf_seq_richness'], errors='coerce')
    warn('EcoBank total_reads missing; using amf_seq_richness as proxy total_reads.')

if 'study_id' not in e.columns:
    study_src = first_present(e, ['sample_code', 'dataset_id', 'paper_id', 'site'])
    if study_src is not None:
        e['study_id'] = e[study_src].astype(str).map(lambda x: f'EcoBank_{x}')
    else:
        e['study_id'] = 'EcoBank_unknown'

    if 'latitude' not in e.columns and 'lat' in e.columns:
        e['latitude'] = e['lat']
    if 'longitude' not in e.columns and 'lon' in e.columns:
        e['longitude'] = e['lon']

for col in ['MAT', 'MAP', 'pH', 'latitude', 'longitude', 'amf_seq_richness', 'amf_genus_richness', 'total_reads']:
    if col not in e.columns:
        e[col] = np.nan

e['source'] = 'EcoBank'
e['plant_species'] = e['plant_species'].map(clean_spaces)
e['plant_species_canon'] = e['plant_species'].map(canonicalize)

for col in ['MAT', 'MAP', 'pH', 'latitude', 'longitude', 'amf_seq_richness', 'amf_genus_richness', 'total_reads']:
    e[col] = pd.to_numeric(e[col], errors='coerce')

e['total_reads'] = e['total_reads'].where(e['total_reads'] > 0, np.nan)
e['log_total_reads'] = np.log(e['total_reads'])

emit(f'EcoBank rows: {len(e):,}')
emit(f"EcoBank unique samples: {e['sample_id'].nunique():,}")
emit(f"EcoBank unique species: {e['plant_species_canon'].dropna().nunique():,}")


# %% [cell 8] type=markdown
# ## Block 4: Unified schema and source diagnostics


# %% [cell 9] type=code
required_schema = [
    'sample_id', 'plant_species', 'plant_species_canon', 'source', 'study_id',
    'latitude', 'longitude', 'MAT', 'MAP', 'pH',
    'amf_seq_richness', 'amf_genus_richness', 'total_reads', 'log_total_reads',
]

def to_schema(df, name):
    out = df.copy()
    for col in required_schema:
        if col not in out.columns:
            out[col] = np.nan
    out = out[required_schema].copy()
    out['sample_id'] = out['sample_id'].astype(str)
    out['study_id'] = out['study_id'].astype(str)
    dups = out['sample_id'].duplicated().sum()
    if dups > 0:
        warn(f'{name} has {dups:,} duplicated sample_id rows (kept as requested).')
    return out

g_std = to_schema(g, 'GlobalAMFungi')
e_std = to_schema(e, 'EcoBank')

def dist_line(df, label):
    seq = pd.to_numeric(df['amf_seq_richness'], errors='coerce')
    reads = pd.to_numeric(df['total_reads'], errors='coerce')
    emit(f'{label} amf_seq_richness mean/median: {seq.mean():.3f}/{seq.median():.3f}')
    emit(f'{label} total_reads mean/median: {reads.mean():.3f}/{reads.median():.3f}')

dist_line(g_std, 'GlobalAMFungi')
dist_line(e_std, 'EcoBank')


# %% [cell 10] type=markdown
# ## Block 5: Load FRED 4.0 and aggregate to one row per species


# %% [cell 11] type=code
# Match 01_data_prep_4 header-row handling
fred_raw = pd.read_csv(FRED_FP, header=None)
fred_raw.columns = fred_raw.iloc[1].astype(str)
fred = fred_raw.iloc[2:].copy().reset_index(drop=True)

# Keep only true data rows (exclude text metadata rows)
if 'Notes_Row ID' in fred.columns:
    fred = fred[pd.to_numeric(fred['Notes_Row ID'], errors='coerce').notna()].copy()

if 'Plant taxonomy_Accepted genus_WFO' not in fred.columns or 'Plant Taxonomy_Accepted species_WFO' not in fred.columns:
    raise ValueError('FRED is missing accepted genus/species columns used in 01_data_prep_4.ipynb.')

fred['fred_species_raw'] = (
    fred['Plant taxonomy_Accepted genus_WFO'].astype(str).map(clean_spaces) + ' ' +
    fred['Plant Taxonomy_Accepted species_WFO'].astype(str).map(clean_spaces)
)
fred['fred_species_raw'] = fred['fred_species_raw'].str.replace(r'\s+', ' ', regex=True).str.strip()
fred['fred_species_raw'] = fred['fred_species_raw'].replace({'nan nan': np.nan})
fred['plant_species_canon'] = fred['fred_species_raw'].map(canonicalize)
fred = fred[fred['plant_species_canon'].notna()].copy()

drop_id_like = ['Notes_Row ID', 'Abbreviated article citation', 'Data source_Citation', 'Data source_DOI', 'Data set_Citation', 'Data set_DOI', 'Notes_Site ID']
trait_cols = [c for c in fred.columns if c not in ['fred_species_raw', 'plant_species_canon'] + drop_id_like]

def agg_one(s):
    num = pd.to_numeric(s, errors='coerce')
    if num.notna().sum() > 0:
        return num.mean()
    s2 = s.dropna().astype(str).str.strip()
    s2 = s2[s2.ne('')]
    if len(s2) == 0:
        return np.nan
    mode = s2.mode()
    return mode.iloc[0] if len(mode) else s2.iloc[0]

fred_agg = fred.groupby('plant_species_canon', as_index=False)[trait_cols].agg(agg_one)
emit(f'FRED rows before aggregation: {len(fred):,}')
emit(f"FRED species before aggregation: {fred['plant_species_canon'].nunique():,}")
emit(f'FRED species after aggregation: {len(fred_agg):,}')
emit(f'Number of FRED trait columns carried: {len(trait_cols):,}')


# %% [cell 12] type=markdown
# ## Block 6: Match diagnostics and merge each source with FRED


# %% [cell 13] type=code
fred_species = set(fred_agg['plant_species_canon'].dropna())

def match_stats(df, name):
    raw_species = set(df['plant_species'].dropna().astype(str).map(clean_spaces))
    can_species = set(df['plant_species_canon'].dropna())
    fred_raw = set(fred['fred_species_raw'].dropna().astype(str).map(clean_spaces))
    exact_raw = len(raw_species & fred_raw)
    can_match = len(can_species & fred_species)
    rate = (can_match / len(can_species) * 100) if len(can_species) else 0.0
    emit(f'{name} exact raw species matches: {exact_raw:,}')
    emit(f'{name} canonical species matches: {can_match:,} / {len(can_species):,} ({rate:.2f}%)')
    return {
        'source': name,
        'raw_exact': exact_raw,
        'canonical_match': can_match,
        'species_total': len(can_species),
        'match_rate_pct': rate,
    }

global_match = match_stats(g_std, 'GlobalAMFungi')
eco_match = match_stats(e_std, 'EcoBank')

g_enriched = g_std.merge(fred_agg, on='plant_species_canon', how='left')
e_enriched = e_std.merge(fred_agg, on='plant_species_canon', how='left')

emit(f'Global rows after FRED merge: {len(g_enriched):,}')
emit(f'EcoBank rows after FRED merge: {len(e_enriched):,}')


# %% [cell 14] type=markdown
# ## Block 7: Final concatenation and QA summaries


# %% [cell 15] type=code
merged = pd.concat([g_enriched, e_enriched], ignore_index=True)

required_final = ['sample_id', 'plant_species', 'source', 'study_id', 'MAT', 'MAP', 'pH', 'amf_seq_richness', 'amf_genus_richness', 'total_reads', 'log_total_reads']
missing_final = [c for c in required_final if c not in merged.columns]
if missing_final:
    raise ValueError(f'Final merged dataset missing required columns: {missing_final}')

trait_core_candidates = [
    'Root N content', 'Root P content', 'Root diameter', 'Root tissue density (RTD)',
    'Specific root area (SRA)', 'Mycorrhiza_Fraction of root length or tips colonized',
    'Soil pH_Water', 'Soil C content', 'Soil N content',
    'Mean annual precipitation (MAP)', 'Mean annual air temperature (MAT)'
]
core_traits = [c for c in trait_core_candidates if c in merged.columns]

if core_traits:
    nonnull_any_core = merged[core_traits].notna().any(axis=1).mean()
else:
    nonnull_any_core = np.nan
    warn('No core trait columns found in merged dataset for non-null trait fraction check.')

emit('--- Final merged dataset QA ---')
emit(f'total rows: {len(merged):,}')
emit(f"unique samples overall: {merged[['source', 'sample_id']].drop_duplicates().shape[0]:,}")
emit(f"unique species overall: {merged['plant_species_canon'].dropna().nunique():,}")
for src, sub in merged.groupby('source'):
    emit(f"{src} rows: {len(sub):,}, unique_samples: {sub['sample_id'].nunique():,}, unique_species: {sub['plant_species_canon'].dropna().nunique():,}")
if pd.notna(nonnull_any_core):
    emit(f'fraction rows with >=1 non-null core trait: {nonnull_any_core:.3f}')

merged.to_csv(MERGED_OUT, index=False)
emit(f'Wrote merged dataset: {MERGED_OUT}')


# %% [cell 16] type=markdown
# ## Block 8: Build markdown summary report


# %% [cell 17] type=code
def top_unmatched(df, n=30):
    x = df[df['plant_species_canon'].notna()].copy()
    x = x[~x['plant_species_canon'].isin(fred_species)]
    out = (x['plant_species']
           .fillna('')
           .astype(str)
           .str.strip()
           .value_counts()
           .head(n))
    return out

miss_rows = []
for c in core_traits:
    miss_rows.append({
        'column': c,
        'missing_count': int(merged[c].isna().sum()),
        'missing_frac': float(merged[c].isna().mean()),
    })
missing_df = pd.DataFrame(miss_rows)

global_unmatched = top_unmatched(g_std, n=30)
eco_unmatched = top_unmatched(e_std, n=30)

lines = []
lines.append('# AMF Traits Merged Summary')
lines.append('')
lines.append('## Inputs used')
lines.append(f'- GlobalAMFungi sample-level: `{global_fp}`')
lines.append(f'- EcoBank source: `{eco_fp}`' if eco_fp is not None else '- EcoBank source: derived from `../Data/ecobank_full.csv` + `../Data/ecobank_cultured.csv`')
lines.append(f'- FRED: `{FRED_FP}`')
lines.append('')

lines.append('## Source counts before/after FRED merge')
lines.append('')
lines.append(f'- GlobalAMFungi before merge rows: {len(g_std):,}')
lines.append(f'- EcoBank before merge rows: {len(e_std):,}')
lines.append(f'- GlobalAMFungi after merge rows: {len(g_enriched):,}')
lines.append(f'- EcoBank after merge rows: {len(e_enriched):,}')
lines.append('')

lines.append('## Unique species per source')
lines.append('')
lines.append(f"- GlobalAMFungi species: {g_std['plant_species_canon'].dropna().nunique():,}")
lines.append(f"- EcoBank species: {e_std['plant_species_canon'].dropna().nunique():,}")
lines.append('')

lines.append('## Species matched to FRED')
lines.append('')
lines.append(f"- GlobalAMFungi exact raw matches: {global_match['raw_exact']:,}")
lines.append(f"- GlobalAMFungi canonical matches: {global_match['canonical_match']:,} / {global_match['species_total']:,} ({global_match['match_rate_pct']:.2f}%)")
lines.append(f"- EcoBank exact raw matches: {eco_match['raw_exact']:,}")
lines.append(f"- EcoBank canonical matches: {eco_match['canonical_match']:,} / {eco_match['species_total']:,} ({eco_match['match_rate_pct']:.2f}%)")
lines.append('')

lines.append('## Final merged dataset dimensions')
lines.append('')
lines.append(f'- Rows: {len(merged):,}')
lines.append(f'- Columns: {merged.shape[1]:,}')
lines.append(f"- Unique samples (source+sample_id): {merged[['source', 'sample_id']].drop_duplicates().shape[0]:,}")
lines.append(f"- Unique species (canonical): {merged['plant_species_canon'].dropna().nunique():,}")
if pd.notna(nonnull_any_core):
    lines.append(f'- Fraction rows with >=1 non-null core trait: {nonnull_any_core:.3f}')
lines.append('')

lines.append('## Missingness summary for key trait columns')
lines.append('')
if len(missing_df):
    lines.append('| column | missing_count | missing_frac |')
    lines.append('|---|---:|---:|')
    for _, r in missing_df.iterrows():
        lines.append(f"| {r['column']} | {int(r['missing_count'])} | {r['missing_frac']:.3f} |")
else:
    lines.append('- No core trait columns available.')
lines.append('')

lines.append('## Top 30 unmatched species')
lines.append('')
lines.append('### GlobalAMFungi unmatched')
if len(global_unmatched):
    for sp, c in global_unmatched.items():
        lines.append(f'- {sp}: {c}')
else:
    lines.append('- None')
lines.append('')
lines.append('### EcoBank unmatched')
if len(eco_unmatched):
    for sp, c in eco_unmatched.items():
        lines.append(f'- {sp}: {c}')
else:
    lines.append('- None')
lines.append('')

# Close-match diagnostics for currently unmatched names
import difflib

def unmatched_set(df):
    s = set(df['plant_species_canon'].dropna().astype(str))
    return sorted(s - fred_species)

def best_close_matches(unmatched_names, fred_names, n=5, cutoff=0.86):
    out = []
    fred_names = sorted(set([x for x in fred_names if isinstance(x, str) and x]))
    fred_genus_map = {}
    for nm in fred_names:
        g = nm.split(' ', 1)[0] if ' ' in nm else nm
        fred_genus_map.setdefault(g, []).append(nm)

    for nm in unmatched_names:
        genus = nm.split(' ', 1)[0] if ' ' in nm else nm
        pool = fred_genus_map.get(genus, fred_names)
        cand = difflib.get_close_matches(nm, pool, n=n, cutoff=cutoff)
        if not cand and pool is not fred_names:
            # fallback to global pool in case genus changed between versions
            cand = difflib.get_close_matches(nm, fred_names, n=n, cutoff=cutoff)
        if cand:
            for c in cand:
                out.append({
                    'unmatched_name': nm,
                    'candidate_fred_name': c,
                    'similarity_ratio': round(difflib.SequenceMatcher(None, nm, c).ratio(), 3),
                    'same_genus': int((nm.split(' ',1)[0] if ' ' in nm else nm) == (c.split(' ',1)[0] if ' ' in c else c))
                })
    out_df = pd.DataFrame(out)
    if len(out_df):
        out_df = out_df.sort_values(['similarity_ratio', 'same_genus', 'unmatched_name'], ascending=[False, False, True]).reset_index(drop=True)
    return out_df

global_unmatched_set = unmatched_set(g_std)
eco_unmatched_set = unmatched_set(e_std)
global_close = best_close_matches(global_unmatched_set, fred_species, n=3, cutoff=0.86)
eco_close = best_close_matches(eco_unmatched_set, fred_species, n=3, cutoff=0.86)

emit(f'GlobalAMFungi unmatched canonical species: {len(global_unmatched_set):,}')
emit(f'EcoBank unmatched canonical species: {len(eco_unmatched_set):,}')
emit(f'GlobalAMFungi close-name suggestions (cutoff=0.86): {len(global_close):,}')
emit(f'EcoBank close-name suggestions (cutoff=0.86): {len(eco_close):,}')

lines.append('## Close-name diagnostics for unmatched species')
lines.append('')
lines.append('- Method: `difflib.get_close_matches` with `SequenceMatcher` ratio cutoff `0.86` (genus-prioritized, then global fallback).')
lines.append('- Why close but not exact: most likely synonym updates, accepted-name changes between taxonomic references/versions, spelling variants, or rank/qualifier differences (e.g., `sp.`, `cf.`).')
lines.append('')
lines.append(f'- GlobalAMFungi unmatched canonical species: {len(global_unmatched_set):,}')
lines.append(f'- EcoBank unmatched canonical species: {len(eco_unmatched_set):,}')
lines.append(f'- GlobalAMFungi close-name suggestions: {len(global_close):,}')
lines.append(f'- EcoBank close-name suggestions: {len(eco_close):,}')
lines.append('')

def append_close_table(lines_obj, title, df_close, top_n=200):
    lines_obj.append(f'### {title}')
    lines_obj.append('')
    if len(df_close) == 0:
        lines_obj.append('- None at current cutoff.')
        lines_obj.append('')
        return
    lines_obj.append('| unmatched_name | candidate_fred_name | similarity_ratio | same_genus |')
    lines_obj.append('|---|---|---:|---:|')
    for _, r in df_close.head(top_n).iterrows():
        lines_obj.append(f"| {r['unmatched_name']} | {r['candidate_fred_name']} | {r['similarity_ratio']:.3f} | {int(r['same_genus'])} |")
    lines_obj.append('')

append_close_table(lines, 'GlobalAMFungi unmatched species with close FRED names', global_close, top_n=200)
append_close_table(lines, 'EcoBank unmatched species with close FRED names', eco_close, top_n=200)

CLOSE_MD_OUT = OUT_DIR / 'amf_traits_close_name_diagnostics.md'
close_lines = []
close_lines.append('# Close-name diagnostics for AMF trait merge')
close_lines.append('')
close_lines.append('This report lists unmatched canonical plant species from GlobalAMFungi and EcoBank that are very similar to a FRED 4.0 canonical species name.')
close_lines.append('')
close_lines.append(f'- GlobalAMFungi unmatched canonical species: {len(global_unmatched_set):,}')
close_lines.append(f'- EcoBank unmatched canonical species: {len(eco_unmatched_set):,}')
close_lines.append(f'- GlobalAMFungi close suggestions (ratio >= 0.86): {len(global_close):,}')
close_lines.append(f'- EcoBank close suggestions (ratio >= 0.86): {len(eco_close):,}')
close_lines.append('')
close_lines.append('Interpretation: close names are typically unresolved synonymy/taxonomic updates across datasets or minor naming/orthography differences. FRED 4.0 uses accepted WFO names, so previously-used names from older references (including many FRED 3-era conventions) may not match exactly.')
close_lines.append('')
append_close_table(close_lines, 'GlobalAMFungi unmatched species with close FRED names', global_close, top_n=1000)
append_close_table(close_lines, 'EcoBank unmatched species with close FRED names', eco_close, top_n=1000)
CLOSE_MD_OUT.write_text('\n'.join(close_lines), encoding='utf-8')
emit(f'Wrote close-name diagnostics: {CLOSE_MD_OUT}')

lines.append('## Warnings')
lines.append('')
if WARN:
    for w in WARN:
        lines.append(f'- {w}')
else:
    lines.append('- None')
lines.append('')

lines.append('## Printed output')
lines.append('')
lines.append('```text')
lines.extend(LOG)
lines.append('```')

SUMMARY_OUT.write_text('\n'.join(lines), encoding='utf-8')
emit(f'Wrote summary report: {SUMMARY_OUT}')
