# Auto-converted from 01_data_prep_global_full.ipynb


# %% [cell 1] type=markdown
# # 01_data_prep_global_full
# 
# Updated GlobalAMFungi data preparation pipeline using **all valid root samples with exactly one dominant plant species**.
# 
# Key change from prior notebook: **no singleton species filtering**.


# %% [cell 2] type=markdown
# ## Block 1: Imports + paths
# 
# Initialize dependencies, file paths, output directory, and logging helpers.


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

assert GLOBAL_FP.exists(), f'Missing file: {GLOBAL_FP}'
assert FRED_FP.exists(), f'Missing file: {FRED_FP}'

LOG_LINES = []
WARNINGS = []

def emit(msg):
    text = str(msg)
    print(text)
    LOG_LINES.append(text)

def warn(msg):
    text = f'WARNING: {msg}'
    print(text)
    LOG_LINES.append(text)
    WARNINGS.append(str(msg))

def require_columns(df: pd.DataFrame, required, df_name):
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f'{df_name} missing required columns: {missing}')

def clean_spaces(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip()
    s = re.sub(r'\s+', ' ', s)
    return s

def canonicalize(s):
    if pd.isna(s):
        return np.nan
    s = str(s).strip().lower()
    s = re.sub(r'\s+', ' ', s)
    s = re.sub(r'[\.;:,]+$', '', s)
    return s

emit(f'Global input: {GLOBAL_FP}')
emit(f'FRED input: {FRED_FP}')
emit(f'Output directory: {OUT_DIR.resolve()}')


# %% [cell 4] type=markdown
# ## Block 2: Load data + initial validation
# 
# Load GlobalAMFungi and FRED, validate required columns, and print initial dimensions.


# %% [cell 5] type=code
global_df = pd.read_csv(GLOBAL_FP)

# Drop unnamed index-like columns if present
unnamed_cols = [c for c in global_df.columns if str(c).startswith('Unnamed')]
if unnamed_cols:
    global_df = global_df.drop(columns=unnamed_cols)

required_global = [
    'id', 'sample_type', 'plants_dominant', 'sequence', 'Genus', 'abundances',
    'latitude', 'longitude', 'MAT', 'MAP', 'pH'
]
require_columns(global_df, required_global, 'GlobalAMFungi')

global_df['plant_species'] = global_df['plants_dominant'].apply(clean_spaces)

emit('--- Initial GlobalAMFungi summary ---')
emit(f'rows: {len(global_df):,}')
emit(f"unique samples: {global_df['id'].nunique():,}")
emit(f"unique plant species: {global_df['plant_species'].dropna().nunique():,}")


# %% [cell 6] type=markdown
# ## Block 3: Root + one-dominant-plant filtering
# 
# Filter to root samples with exactly one dominant plant species and non-null `plants_dominant`.


# %% [cell 7] type=code
stage_stats = []

def count_plants(x):
    if pd.isna(x):
        return 0
    parts = [p.strip() for p in str(x).split(';') if p.strip()]
    return len(parts)

def add_stage(name, df):
    n_species = df['plant_species'].dropna().nunique() if 'plant_species' in df.columns else 0
    stage_stats.append({
        'stage': name,
        'rows': int(len(df)),
        'unique_samples': int(df['id'].nunique()) if 'id' in df.columns else 0,
        'unique_plant_species': int(n_species),
    })
    emit(f"{name}: rows={len(df):,}, unique_samples={df['id'].nunique():,}, unique_plant_species={n_species:,}")

g0 = global_df.copy()
add_stage('initial', g0)

g1 = g0[g0['sample_type'].astype(str).str.strip().str.lower().eq('root')].copy()
g1['n_plants'] = g1['plants_dominant'].apply(count_plants)
g1['plant_species'] = g1['plants_dominant'].apply(clean_spaces)
add_stage('root_only', g1)

non_null_plant = g1['plant_species'].notna() & g1['plant_species'].astype(str).str.strip().ne('')
g2 = g1[non_null_plant & (g1['n_plants'] == 1)].copy()

# IMPORTANT: No singleton species filtering here (removed by design).
g_filtered = g2.copy()
add_stage('root_singleplant_nonnull', g_filtered)

g_filtered['plant_species_canonical'] = g_filtered['plant_species'].map(canonicalize)

stage_df = pd.DataFrame(stage_stats)
stage_df


# %% [cell 8] type=markdown
# ## Block 4: Sample-level AMF richness metrics
# 
# Compute sequence/genus richness and total reads per sample.


# %% [cell 9] type=code
g_filtered['abundances_num'] = pd.to_numeric(g_filtered['abundances'], errors='coerce').fillna(0)
g_filtered['Genus_filled'] = g_filtered['Genus'].fillna('Unknown').astype(str)

for c in ['latitude', 'longitude', 'MAT', 'MAP', 'pH']:
    g_filtered[c] = pd.to_numeric(g_filtered[c], errors='coerce')

meta_cols = ['plant_species', 'latitude', 'longitude', 'MAT', 'MAP', 'pH']
sample_meta = g_filtered.groupby('id', as_index=False)[meta_cols].first()

sample_metrics = g_filtered.groupby('id', as_index=False).agg(
    amf_seq_richness=('sequence', 'nunique'),
    amf_genus_richness=('Genus_filled', 'nunique'),
    total_reads=('abundances_num', 'sum')
)

sample_level = sample_meta.merge(sample_metrics, on='id', how='left')
sample_level = sample_level.sort_values('id').reset_index(drop=True)

sample_out = OUT_DIR / 'globalamfungi_sample_level_full.csv'
sample_level.to_csv(sample_out, index=False)

emit('--- Sample-level summary ---')
emit(f'rows: {len(sample_level):,}')
emit(f"mean amf_seq_richness: {sample_level['amf_seq_richness'].mean():.3f}")
emit(f"median amf_seq_richness: {sample_level['amf_seq_richness'].median():.3f}")
emit(f"mean total_reads: {sample_level['total_reads'].mean():.3f}")
emit(f'Wrote: {sample_out}')
sample_level.head()


# %% [cell 10] type=markdown
# ## Block 5: Plant-species-level summary
# 
# Aggregate sample-level richness by plant species.


# %% [cell 11] type=code
species_level = sample_level.groupby('plant_species', as_index=False).agg(
    number_of_samples=('id', 'nunique'),
    mean_amf_seq_richness=('amf_seq_richness', 'mean'),
    std_amf_seq_richness=('amf_seq_richness', 'std')
)
species_level = species_level.sort_values('plant_species').reset_index(drop=True)

species_out = OUT_DIR / 'globalamfungi_species_level_full.csv'
species_level.to_csv(species_out, index=False)

emit('--- Plant-species-level summary ---')
emit(f'total plant species: {len(species_level):,}')
emit(f"total samples: {sample_level['id'].nunique():,}")
emit(f'Wrote: {species_out}')
species_level.head()


# %% [cell 12] type=markdown
# ## Block 6: Load and prepare FRED 4.0
# 
# Use the same header-row approach as the reference notebook (`header=1`), then keep data rows only.


# %% [cell 13] type=code
fred = pd.read_csv(FRED_FP, header=1)

# Remove descriptive/unit rows by requiring numeric Row ID
fred['Notes_Row ID_num'] = pd.to_numeric(fred.get('Notes_Row ID'), errors='coerce')
fred = fred[fred['Notes_Row ID_num'].notna()].copy()

require_columns(
    fred,
    ['Plant taxonomy_Accepted genus_WFO', 'Plant Taxonomy_Accepted species_WFO'],
    'FRED'
)

def fred_binomial(row):
    g = clean_spaces(row.get('Plant taxonomy_Accepted genus_WFO'))
    s = clean_spaces(row.get('Plant Taxonomy_Accepted species_WFO'))
    if pd.isna(g) or pd.isna(s) or g == '' or s == '':
        return np.nan
    return f'{g} {s}'

fred['fred_species_raw'] = fred.apply(fred_binomial, axis=1)
fred['fred_species_canonical'] = fred['fred_species_raw'].map(canonicalize)

emit(f'FRED rows (data-only): {len(fred):,}')
emit(f"FRED unique species (canonical): {fred['fred_species_canonical'].dropna().nunique():,}")


# %% [cell 14] type=markdown
# ## Block 7: GlobalAMFungi ↔ FRED matching
# 
# Compute raw and canonical species matches, plus genus overlap.


# %% [cell 15] type=code
global_species_raw = set(sample_level['plant_species'].dropna().astype(str))
global_species_can = set(sample_level['plant_species'].dropna().map(canonicalize))
global_species_can = {x for x in global_species_can if x}

fred_species_raw = set(fred['fred_species_raw'].dropna().astype(str))
fred_species_can = set(fred['fred_species_canonical'].dropna().astype(str))

exact_matches_raw = global_species_raw & fred_species_raw
matches_canonical = global_species_can & fred_species_can

global_genera = {s.split(' ', 1)[0] for s in global_species_can if ' ' in s}
fred_genera = {s.split(' ', 1)[0] for s in fred_species_can if ' ' in s}
genus_overlap = global_genera & fred_genera

total_global_species = len(global_species_can)
total_fred_species = len(fred_species_can)
match_pct = (len(matches_canonical) / total_global_species * 100) if total_global_species else 0.0

match_counts = {
    'total_global_species': total_global_species,
    'total_fred_species': total_fred_species,
    'exact_matches_raw': len(exact_matches_raw),
    'matches_canonical': len(matches_canonical),
    'genus_overlap': len(genus_overlap),
    'match_percentage': match_pct,
}

emit('--- GlobalAMFungi ↔ FRED matching ---')
emit(f'global_unique_species_after_filters: {total_global_species:,}')
emit(f'fred_unique_species: {total_fred_species:,}')
emit(f'exact_matches_raw: {len(exact_matches_raw):,}')
emit(f'canonical_matches: {len(matches_canonical):,}')
emit(f'genus_overlap: {len(genus_overlap):,}')
emit(f'match_percentage: {match_pct:.2f}%')

pd.DataFrame([match_counts])


# %% [cell 16] type=markdown
# ## Block 8: Figures
# 
# Generate required histogram, scatter, and global point map.


# %% [cell 17] type=code
# Figure 1: Histogram of AMF sequence richness
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(sample_level['amf_seq_richness'].dropna(), bins=30, color='#2b8cbe', edgecolor='black')
ax.set_title('AMF sequence richness distribution (filtered root samples)')
ax.set_xlabel('amf_seq_richness')
ax.set_ylabel('Number of samples')
hist_fp = OUT_DIR / 'globalamf_hist_amf_seq_richness_full.png'
fig.tight_layout()
fig.savefig(hist_fp, dpi=150)
plt.close(fig)
emit(f'Wrote figure: {hist_fp}')

# Figure 2: Scatter richness vs MAT with linear trend
scatter_df = sample_level[['MAT', 'amf_seq_richness']].dropna().copy()
fig, ax = plt.subplots(figsize=(8, 5))
ax.scatter(scatter_df['MAT'], scatter_df['amf_seq_richness'], alpha=0.7, s=14, color='#238b45')
if len(scatter_df) >= 2:
    x = scatter_df['MAT'].to_numpy()
    y = scatter_df['amf_seq_richness'].to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    ax.plot(xs, slope*xs + intercept, color='black', linewidth=2)
else:
    warn('Too few points to fit MAT regression line.')
ax.set_title('AMF sequence richness vs MAT')
ax.set_xlabel('MAT')
ax.set_ylabel('amf_seq_richness')
scatter_fp = OUT_DIR / 'globalamf_scatter_richness_vs_MAT_full.png'
fig.tight_layout()
fig.savefig(scatter_fp, dpi=150)
plt.close(fig)
emit(f'Wrote figure: {scatter_fp}')

# Figure 3: Simple global map scatter (lon vs lat)
map_df = sample_level[['longitude', 'latitude']].dropna().copy()
fig, ax = plt.subplots(figsize=(10, 5))
ax.scatter(map_df['longitude'], map_df['latitude'], s=8, alpha=0.6, color='#88419d')
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
ax.set_title('GlobalAMFungi filtered sample coordinates')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.grid(alpha=0.25)
map_fp = OUT_DIR / 'globalamf_map_points_full.png'
fig.tight_layout()
fig.savefig(map_fp, dpi=150)
plt.close(fig)
emit(f'Wrote figure: {map_fp}')


# %% [cell 18] type=markdown
# ## Block 9: Write markdown summary report
# 
# Save stage dimensions, statistics, matching metrics, warnings, and captured notebook output.


# %% [cell 19] type=code
summary_fp = OUT_DIR / 'globalamfungi_data_prep_summary_full.md'
desc = sample_level[['amf_seq_richness', 'amf_genus_richness', 'total_reads', 'MAT', 'MAP', 'pH']].describe().T

lines = []
lines.append('# GlobalAMFungi Data Prep Summary (Full)')
lines.append('')
lines.append('## Input files')
lines.append(f'- GlobalAMFungi: `{GLOBAL_FP}`')
lines.append(f'- FRED: `{FRED_FP}`')
lines.append('')

lines.append('## Dataset dimensions by stage')
lines.append('')
lines.append('| stage | rows | unique_samples | unique_plant_species |')
lines.append('|---|---:|---:|---:|')
for _, r in stage_df.iterrows():
    lines.append(f"| {r['stage']} | {int(r['rows'])} | {int(r['unique_samples'])} | {int(r['unique_plant_species'])} |")
lines.append('')

lines.append('## Richness statistics (sample-level)')
lines.append('')
lines.append('| metric | count | mean | std | min | 25% | 50% | 75% | max |')
lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
for metric, row in desc.iterrows():
    vals = [row.get(k, np.nan) for k in ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']]
    vals_fmt = [f"{v:.4g}" if pd.notna(v) else '' for v in vals]
    lines.append(f'| {metric} | ' + ' | '.join(vals_fmt) + ' |')
lines.append('')

lines.append('## FRED matching statistics')
lines.append('')
for k in ['total_global_species', 'total_fred_species', 'exact_matches_raw', 'matches_canonical', 'genus_overlap']:
    lines.append(f'- {k}: **{match_counts[k]:,}**')
lines.append(f"- match_percentage: **{match_counts['match_percentage']:.2f}%**")
lines.append('')

lines.append('## Warnings / caveats')
lines.append('')
if WARNINGS:
    for w in WARNINGS:
        lines.append(f'- {w}')
else:
    lines.append('- None')
lines.append('')

lines.append('## Exact printed notebook output')
lines.append('')
lines.append('```text')
lines.extend(LOG_LINES)
lines.append('```')

summary_fp.write_text('\n'.join(lines), encoding='utf-8')
emit(f'Wrote markdown summary: {summary_fp}')


# %% [cell 20] type=markdown
# ## Block 10: Final required summary output


# %% [cell 21] type=code
final_samples = sample_level['id'].nunique()
final_species = sample_level['plant_species'].nunique()
final_matches = match_counts['matches_canonical']
final_pct = match_counts['match_percentage']

emit('Final dataset summary:')
emit(f'Total samples: {final_samples:,}')
emit(f'Total plant species: {final_species:,}')
emit(f'Total GlobalAMFungi–FRED matches: {final_matches:,}')
emit(f'Match percentage: {final_pct:.2f}%')

# Refresh summary so it includes final printed output as well
lines = []
lines.append('# GlobalAMFungi Data Prep Summary (Full)')
lines.append('')
lines.append('## Input files')
lines.append(f'- GlobalAMFungi: `{GLOBAL_FP}`')
lines.append(f'- FRED: `{FRED_FP}`')
lines.append('')
lines.append('## Dataset dimensions by stage')
lines.append('')
lines.append('| stage | rows | unique_samples | unique_plant_species |')
lines.append('|---|---:|---:|---:|')
for _, r in stage_df.iterrows():
    lines.append(f"| {r['stage']} | {int(r['rows'])} | {int(r['unique_samples'])} | {int(r['unique_plant_species'])} |")
lines.append('')
lines.append('## Richness statistics (sample-level)')
lines.append('')
lines.append('| metric | count | mean | std | min | 25% | 50% | 75% | max |')
lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
for metric, row in desc.iterrows():
    vals = [row.get(k, np.nan) for k in ['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max']]
    vals_fmt = [f"{v:.4g}" if pd.notna(v) else '' for v in vals]
    lines.append(f'| {metric} | ' + ' | '.join(vals_fmt) + ' |')
lines.append('')
lines.append('## FRED matching statistics')
lines.append('')
for k in ['total_global_species', 'total_fred_species', 'exact_matches_raw', 'matches_canonical', 'genus_overlap']:
    lines.append(f'- {k}: **{match_counts[k]:,}**')
lines.append(f"- match_percentage: **{match_counts['match_percentage']:.2f}%**")
lines.append('')
lines.append('## Warnings / caveats')
lines.append('')
if WARNINGS:
    for w in WARNINGS:
        lines.append(f'- {w}')
else:
    lines.append('- None')
lines.append('')
lines.append('## Exact printed notebook output')
lines.append('')
lines.append('```text')
lines.extend(LOG_LINES)
lines.append('```')
summary_fp.write_text('\n'.join(lines), encoding='utf-8')
emit(f'Refreshed markdown summary: {summary_fp}')
