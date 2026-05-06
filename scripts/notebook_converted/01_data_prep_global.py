# Auto-converted from 01_data_prep_global.ipynb


# %% [cell 1] type=markdown
# # GlobalAMFungi Data Prep
# 
# This notebook reproduces the analysis style from `01_data_prep_4.ipynb`, adapted to GlobalAMFungi input (`../Data/globalamf.csv`) and FRED 4.0 (`../Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv`).


# %% [cell 2] type=markdown
# ## Block 1: Imports + paths (run once)
# 
# Sets paths, output directory, and helper logging utilities.


# %% [cell 3] type=code
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA_DIR = Path('../Data')
OUT_DIR = Path('../Output')
OUT_DIR.mkdir(parents=True, exist_ok=True)

GLOBAL_FP = DATA_DIR / 'globalamf.csv'
FRED_FP = DATA_DIR / 'FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv'

assert GLOBAL_FP.exists(), f'Missing input: {GLOBAL_FP}'
assert FRED_FP.exists(), f'Missing input: {FRED_FP}'

OUTPUT_LOG = []
WARNINGS = []

def log(msg: str):
    print(msg)
    OUTPUT_LOG.append(str(msg))

def warn(msg: str):
    line = f'WARNING: {msg}'
    print(line)
    OUTPUT_LOG.append(line)
    WARNINGS.append(msg)

def require_columns(df: pd.DataFrame, cols, df_name='dataframe'):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f'{df_name} is missing required columns: {missing}')

log(f'Global input: {GLOBAL_FP}')
log(f'FRED input: {FRED_FP}')
log(f'Output dir: {OUT_DIR.resolve()}')


# %% [cell 4] type=markdown
# ## Block 2: Load GlobalAMFungi + FRED
# 
# FRED is loaded using the same special handling pattern as the original notebook: `header=1` because the first row is field-code metadata.


# %% [cell 5] type=code
global_df = pd.read_csv(GLOBAL_FP)

# Drop unnamed helper/index columns if present
for c in ['Unnamed: 0', 'Unnamed: 0.1', '']:
    if c in global_df.columns:
        global_df = global_df.drop(columns=[c])

fred = pd.read_csv(FRED_FP, header=1)

# Keep only true data rows in FRED (exclude units/metadata rows).
fred['Notes_Row ID_num'] = pd.to_numeric(fred.get('Notes_Row ID'), errors='coerce')
fred = fred[fred['Notes_Row ID_num'].notna()].copy()

log(f'Global raw shape: {global_df.shape}')
log(f'FRED raw (header=1, data-only) shape: {fred.shape}')

require_columns(global_df, ['id', 'sample_type', 'plants_dominant', 'sequence', 'Genus', 'abundances'], 'global_df')
require_columns(
    fred,
    ['Plant taxonomy_Accepted genus_WFO', 'Plant Taxonomy_Accepted species_WFO'],
    'fred'
)


# %% [cell 6] type=markdown
# ## Block 3: Root + single-plant filtering
# 
# Implements the required filtering steps and logs dimensions at each checkpoint.


# %% [cell 7] type=code
def clean_spaces(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    s = re.sub(r'\s+', ' ', s)
    return s

stage_stats = []

def snapshot(name, df):
    n_samples = df['id'].nunique() if 'id' in df.columns else np.nan
    n_species = df['plant_species'].nunique(dropna=True) if 'plant_species' in df.columns else np.nan
    rec = {'stage': name, 'rows': int(len(df)), 'unique_samples': int(n_samples) if pd.notna(n_samples) else np.nan, 'unique_plant_species': int(n_species) if pd.notna(n_species) else np.nan}
    stage_stats.append(rec)
    log(f"{name}: rows={rec['rows']:,}, unique_samples={rec['unique_samples']:,}, unique_plant_species={rec['unique_plant_species']:,}")

g0 = global_df.copy()
g0['plant_species'] = g0['plants_dominant'].apply(clean_spaces)
snapshot('initial', g0)

g1 = g0[g0['sample_type'].astype(str).str.strip().str.lower().eq('root')].copy()
snapshot('root_only', g1)

# n_plants from semicolon-separated plants_dominant
def count_plants(x):
    if pd.isna(x):
        return 0
    tokens = [t.strip() for t in str(x).split(';') if t.strip()]
    return len(tokens)

g1['n_plants'] = g1['plants_dominant'].apply(count_plants)
g1['plant_species'] = g1['plants_dominant'].apply(clean_spaces)

g2 = g1[g1['plant_species'].notna() & g1['plant_species'].ne('') & (g1['n_plants'] == 1)].copy()
snapshot('root_singleplant', g2)

# Keep only plant species that appear in exactly one unique sample id
species_sample_counts = g2.groupby('plant_species')['id'].nunique()
singleton_species = species_sample_counts[species_sample_counts == 1].index
g3 = g2[g2['plant_species'].isin(singleton_species)].copy()
snapshot('root_singleplant_singleton_species', g3)

if g3.empty:
    raise ValueError('Filtered dataset is empty after singleton-species step.')

stage_df = pd.DataFrame(stage_stats)
stage_df


# %% [cell 8] type=markdown
# ## Block 4: Sample-level AMF summaries


# %% [cell 9] type=code
g3['abundances_num'] = pd.to_numeric(g3['abundances'], errors='coerce').fillna(0)
g3['Genus_filled'] = g3['Genus'].fillna('Unknown').astype(str)

meta_cols = ['paper_id', 'latitude', 'longitude', 'continent', 'Biome', 'MAT', 'MAP', 'pH', 'year_of_sampling', 'plant_species']
existing_meta_cols = [c for c in meta_cols if c in g3.columns]

sample_meta = g3.groupby('id', as_index=False)[existing_meta_cols].first()
sample_metrics = g3.groupby('id', as_index=False).agg(
    amf_seq_richness=('sequence', 'nunique'),
    amf_genus_richness=('Genus_filled', 'nunique'),
    total_reads=('abundances_num', 'sum')
)

sample_level = sample_meta.merge(sample_metrics, on='id', how='left')
sample_level = sample_level.sort_values('id').reset_index(drop=True)

log(f'Sample-level table rows: {len(sample_level):,}')

sample_out = OUT_DIR / 'globalamfungi_sample_level.csv'
sample_level.to_csv(sample_out, index=False)
log(f'Wrote: {sample_out}')

sample_level.head()


# %% [cell 10] type=markdown
# ## Block 5: Plant-species-level summary


# %% [cell 11] type=code
species_level = sample_level.groupby('plant_species', as_index=False).agg(
    n_samples=('id', 'nunique'),
    mean_amf_seq_richness=('amf_seq_richness', 'mean'),
    mean_amf_genus_richness=('amf_genus_richness', 'mean'),
    mean_total_reads=('total_reads', 'mean'),
    mean_MAT=('MAT', 'mean'),
    mean_MAP=('MAP', 'mean'),
    mean_pH=('pH', 'mean')
)

species_level = species_level.sort_values('plant_species').reset_index(drop=True)
species_out = OUT_DIR / 'globalamfungi_species_level.csv'
species_level.to_csv(species_out, index=False)
log(f'Plant-species summary rows: {len(species_level):,}')
log(f'Wrote: {species_out}')

species_level.head()


# %% [cell 12] type=markdown
# ## Block 6: FRED matching (raw and canonical)


# %% [cell 13] type=code
def build_binomial(genus, species):
    g = clean_spaces(genus)
    s = clean_spaces(species)
    if pd.isna(g) or pd.isna(s) or g == '' or s == '':
        return np.nan
    return f'{g} {s}'

fred['fred_species_raw'] = fred.apply(
    lambda r: build_binomial(r.get('Plant taxonomy_Accepted genus_WFO'), r.get('Plant Taxonomy_Accepted species_WFO')),
    axis=1
)

def canonical_binomial(s):
    if pd.isna(s):
        return np.nan
    s = str(s).strip().lower()
    s = re.sub(r'\s+', ' ', s)
    s = re.sub(r'[\.;:,]+$', '', s)
    return s

global_species_raw = set(g3['plant_species'].dropna().astype(str))
fred_species_raw = set(fred['fred_species_raw'].dropna().astype(str))

raw_match = global_species_raw & fred_species_raw

global_species_can = set(g3['plant_species'].dropna().map(canonical_binomial))
fred_species_can = set(fred['fred_species_raw'].dropna().map(canonical_binomial))

global_species_can = {x for x in global_species_can if x}
fred_species_can = {x for x in fred_species_can if x}
can_match = global_species_can & fred_species_can

global_genera = {x.split(' ', 1)[0] for x in global_species_can if ' ' in x}
fred_genera = {x.split(' ', 1)[0] for x in fred_species_can if ' ' in x}
genus_overlap = global_genera & fred_genera

match_counts = {
    'global_unique_species_after_filters': len(global_species_raw),
    'fred_unique_species_raw': len(fred_species_raw),
    'raw_exact_matches': len(raw_match),
    'canonical_matches': len(can_match),
    'global_unique_genera_canonical': len(global_genera),
    'fred_unique_genera_canonical': len(fred_genera),
    'genus_overlap_canonical': len(genus_overlap),
}

log('--- GlobalAMFungi ↔ FRED matching ---')
for k, v in match_counts.items():
    log(f'{k}: {v:,}')

pd.DataFrame([match_counts]).T.rename(columns={0: 'value'})


# %% [cell 14] type=markdown
# ## Block 7: Figures


# %% [cell 15] type=code
# Histogram: AMF sequence richness
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(sample_level['amf_seq_richness'].dropna(), bins=20, color='#2b8cbe', edgecolor='black')
ax.set_title('GlobalAMFungi root samples: AMF sequence richness')
ax.set_xlabel('amf_seq_richness')
ax.set_ylabel('Count of samples')
hist_fp = OUT_DIR / 'globalamf_hist_amf_seq_richness.png'
fig.tight_layout()
fig.savefig(hist_fp, dpi=150)
plt.close(fig)
log(f'Wrote figure: {hist_fp}')

# Scatter: richness vs MAT with linear fit
scatter_df = sample_level[['MAT', 'amf_seq_richness']].copy()
scatter_df['MAT'] = pd.to_numeric(scatter_df['MAT'], errors='coerce')
scatter_df['amf_seq_richness'] = pd.to_numeric(scatter_df['amf_seq_richness'], errors='coerce')
scatter_df = scatter_df.dropna()

fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(scatter_df['MAT'], scatter_df['amf_seq_richness'], alpha=0.7, s=20, color='#238b45')

if len(scatter_df) >= 2:
    x = scatter_df['MAT'].to_numpy()
    y = scatter_df['amf_seq_richness'].to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ys = slope * xs + intercept
    ax.plot(xs, ys, color='black', linewidth=2, label=f'Linear fit: y={slope:.3f}x+{intercept:.3f}')
    ax.legend(loc='best')
else:
    warn('Too few non-null points for MAT vs richness regression line.')

ax.set_title('AMF sequence richness vs MAT')
ax.set_xlabel('MAT')
ax.set_ylabel('amf_seq_richness')
scatter_fp = OUT_DIR / 'globalamf_scatter_richness_vs_MAT.png'
fig.tight_layout()
fig.savefig(scatter_fp, dpi=150)
plt.close(fig)
log(f'Wrote figure: {scatter_fp}')

# Simple global point map (lat/lon scatter)
map_df = sample_level[['longitude', 'latitude']].copy()
map_df['longitude'] = pd.to_numeric(map_df['longitude'], errors='coerce')
map_df['latitude'] = pd.to_numeric(map_df['latitude'], errors='coerce')
map_df = map_df.dropna()

if len(map_df) > 0:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.scatter(map_df['longitude'], map_df['latitude'], s=12, alpha=0.7, color='#88419d')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('GlobalAMFungi filtered root samples (simple global scatter)')
    ax.grid(alpha=0.3)
    map_fp = OUT_DIR / 'globalamf_map_points.png'
    fig.tight_layout()
    fig.savefig(map_fp, dpi=150)
    plt.close(fig)
    log(f'Wrote figure: {map_fp}')
else:
    warn('No finite latitude/longitude available for map figure.')


# %% [cell 16] type=markdown
# ## Block 8: Descriptive statistics


# %% [cell 17] type=code
desc_cols = ['amf_seq_richness', 'amf_genus_richness', 'total_reads', 'MAT', 'MAP', 'pH']
desc = sample_level[desc_cols].apply(pd.to_numeric, errors='coerce').describe().T
desc_out = OUT_DIR / 'globalamf_descriptive_stats.csv'
desc.to_csv(desc_out)
log(f'Wrote descriptive stats: {desc_out}')
desc


# %% [cell 18] type=markdown
# ## Block 9: Write markdown summary report


# %% [cell 19] type=code
summary_fp = OUT_DIR / 'globalamfungi_data_prep_summary.md'

lines = []
lines.append('# GlobalAMFungi Data Prep Summary')
lines.append('')
lines.append('## Inputs')
lines.append(f'- GlobalAMFungi: `{GLOBAL_FP}`')
lines.append(f'- FRED: `{FRED_FP}`')
lines.append('')

lines.append('## Dataset dimensions at each filter step')
lines.append('')
lines.append('| stage | rows | unique_samples | unique_plant_species |')
lines.append('|---|---:|---:|---:|')
for _, r in stage_df.iterrows():
    us = 0 if pd.isna(r['unique_samples']) else int(r['unique_samples'])
    ups = 0 if pd.isna(r['unique_plant_species']) else int(r['unique_plant_species'])
    lines.append(f"| {r['stage']} | {int(r['rows'])} | {us} | {ups} |")
lines.append('')

lines.append('## Descriptive statistics (sample-level)')
lines.append('')
lines.append('| metric | count | mean | std | min | 25% | 50% | 75% | max |')
lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
for metric, row in desc.iterrows():
    vals = [row.get(k, np.nan) for k in ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']]
    fmt = [f"{v:.4g}" if pd.notna(v) else '' for v in vals]
    lines.append(f'| {metric} | ' + ' | '.join(fmt) + ' |')
lines.append('')

lines.append('## GlobalAMFungi ↔ FRED matching')
lines.append('')
for k, v in match_counts.items():
    lines.append(f'- {k}: **{v:,}**')
lines.append('')

lines.append('## Warnings / caveats')
lines.append('')
if WARNINGS:
    for w in WARNINGS:
        lines.append(f'- {w}')
else:
    lines.append('- None')
lines.append('')

lines.append('## Exact text output printed in notebook')
lines.append('')
lines.append('```text')
lines.extend(OUTPUT_LOG)
lines.append('```')

summary_fp.write_text('\n'.join(lines), encoding='utf-8')
log(f'Wrote markdown summary: {summary_fp}')

print('--- FINAL GlobalAMFungi ↔ FRED match counts ---')
for k in ['global_unique_species_after_filters', 'raw_exact_matches', 'canonical_matches', 'genus_overlap_canonical']:
    print(f'{k}: {match_counts[k]:,}')
