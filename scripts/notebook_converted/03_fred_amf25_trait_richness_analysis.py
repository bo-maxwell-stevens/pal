# Auto-converted from 03_fred_amf25_trait_richness_analysis.ipynb


# %% [cell 1] type=markdown
# # 03_fred_amf25_trait_richness_analysis
# 
# Exploratory analyses and publication-quality figures for merged AMF-FRED sample-level dataset.
# 
# Reference intent adapted from `FRED.Rmd` and `FRED-AMF25 (1).docx`:
# - Root economics trait relationships with AMF richness
# - Woody vs nonwoody contrasts
# - Climate covariates and sequencing-depth-aware models


# %% [cell 2] type=markdown
# ## Block 1: Imports, paths, settings


# %% [cell 3] type=code
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from statsmodels.formula.api import ols
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

DATA_FP = Path('../Output/amf_traits_merged_sample_level.csv')
FIG_DIR = Path('../Output/Figures')
FIG_DIR.mkdir(parents=True, exist_ok=True)
SUMMARY_FP = Path('../Output/fred_amf25_results_summary.md')

assert DATA_FP.exists(), f'Missing input dataset: {DATA_FP}'

plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

COLORS = {'GlobalAMFungi': '#1f77b4', 'EcoBank': '#d62728'}

LOG = []

def emit(msg):
    s = str(msg)
    print(s)
    LOG.append(s)

def canon(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    s = ' '.join(s.split())
    return s

trait_cols = [
    'Root diameter',
    'Root tissue density (RTD)',
    'Specific root area (SRA)',
    'Root N content',
    'Root P content',
]

core_required = [
    'amf_seq_richness', 'amf_genus_richness', 'log_total_reads', 'total_reads',
    'MAT', 'MAP', 'sample_id', 'plant_species', 'plant_species_canon',
    'source', 'study_id', 'latitude', 'longitude', 'Plant woodiness_TRY'
] + trait_cols

emit(f'Input dataset: {DATA_FP}')
emit(f'Figure directory: {FIG_DIR.resolve()}')
emit(f'Summary markdown: {SUMMARY_FP.resolve()}')


# %% [cell 4] type=markdown
# ## Block 2: Load and preprocess data


# %% [cell 5] type=code
df = pd.read_csv(DATA_FP)

missing_required = [c for c in core_required if c not in df.columns]
if missing_required:
    raise ValueError(f'Missing required columns: {missing_required}')

# Coerce numerics used in analyses
for c in ['amf_seq_richness', 'amf_genus_richness', 'log_total_reads', 'total_reads', 'MAT', 'MAP', 'latitude', 'longitude'] + trait_cols + ['Mycorrhiza_Fraction of root length or tips colonized']:
    if c in df.columns:
        df[c] = pd.to_numeric(df[c], errors='coerce')

df['log_richness'] = np.log1p(df['amf_seq_richness'])
df['log_genus_richness'] = np.log1p(df['amf_genus_richness'])

df['plant_species_canon'] = df['plant_species_canon'].map(canon)
df['genus'] = df['plant_species_canon'].fillna('').astype(str).str.split().str[0].replace({'': np.nan})

wood_raw = df['Plant woodiness_TRY'].fillna('').astype(str).str.lower()
is_woody = wood_raw.str.contains('woody|tree|shrub', regex=True)
df['woodiness_group'] = np.where(is_woody, 'woody', 'nonwoody')

emit(f'Total rows: {len(df):,}')
emit(f"Sources: {df['source'].value_counts(dropna=False).to_dict()}")
emit(f"Unique species: {df['plant_species_canon'].nunique(dropna=True):,}")
emit(f"Unique genera: {df['genus'].nunique(dropna=True):,}")


# %% [cell 6] type=markdown
# ## Block 3: Plot helpers


# %% [cell 7] type=code
def savefig(fp):
    plt.tight_layout()
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close()
    emit(f'Wrote figure: {fp}')

def sample_trait_plot(data, xcol, ycol, fp, xlabel=None, ylabel='log(1 + AMF sequence richness)'):
    use = data[[xcol, ycol, 'source']].dropna()
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    for src, sub in use.groupby('source'):
        c = COLORS.get(src, '#555555')
        ax.scatter(sub[xcol], sub[ycol], s=14, alpha=0.35, color=c, label=src)
        if len(sub) >= 15:
            sm = lowess(sub[ycol].to_numpy(), sub[xcol].to_numpy(), frac=0.35, return_sorted=True)
            ax.plot(sm[:, 0], sm[:, 1], color=c, linewidth=2.2)
    ax.set_xlabel(xlabel if xlabel else xcol)
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False)
    savefig(fp)

def species_trait_plot(data, xcol, fp, title=None):
    use = data[['plant_species_canon', xcol, 'log_richness']].dropna()
    agg = use.groupby('plant_species_canon', as_index=False).agg(
        mean_trait=(xcol, 'mean'),
        mean_log_richness=('log_richness', 'mean'),
        sample_count=('log_richness', 'size')
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    sizes = 12 + 3.0 * np.sqrt(agg['sample_count'].to_numpy())
    ax.scatter(agg['mean_trait'], agg['mean_log_richness'], s=sizes, alpha=0.45, color='#2a9d8f', edgecolor='none')
    if len(agg) >= 15:
        sm = lowess(agg['mean_log_richness'].to_numpy(), agg['mean_trait'].to_numpy(), frac=0.45, return_sorted=True)
        ax.plot(sm[:, 0], sm[:, 1], color='black', linewidth=2)
    ax.set_xlabel(xcol)
    ax.set_ylabel('Species mean log(1 + AMF sequence richness)')
    if title:
        ax.set_title(title)
    savefig(fp)

def species_trait_by_wood_plot(data, xcol, wood_group, fp):
    use = data[data['woodiness_group'].eq(wood_group)][['plant_species_canon', xcol, 'log_richness']].dropna()
    agg = use.groupby('plant_species_canon', as_index=False).agg(
        mean_trait=(xcol, 'mean'),
        mean_log_richness=('log_richness', 'mean'),
        sample_count=('log_richness', 'size')
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    if len(agg):
        sizes = 14 + 3.2 * np.sqrt(agg['sample_count'].to_numpy())
        color = '#8c564b' if wood_group == 'woody' else '#17becf'
        ax.scatter(agg['mean_trait'], agg['mean_log_richness'], s=sizes, alpha=0.5, color=color, edgecolor='none')
        if len(agg) >= 10:
            sm = lowess(agg['mean_log_richness'].to_numpy(), agg['mean_trait'].to_numpy(), frac=0.5, return_sorted=True)
            ax.plot(sm[:, 0], sm[:, 1], color='black', linewidth=2)
    ax.set_xlabel(xcol)
    ax.set_ylabel(f'{wood_group.capitalize()} species mean log(1 + richness)')
    savefig(fp)


# %% [cell 8] type=markdown
# ## Block 4: Figures 1-3 (map and distributions)


# %% [cell 9] type=code
# Figure 1: global sample map
map_df = df[['longitude', 'latitude', 'source']].dropna()
fig, ax = plt.subplots(figsize=(10, 5.2))
for src, sub in map_df.groupby('source'):
    ax.scatter(sub['longitude'], sub['latitude'], s=10, alpha=0.45, color=COLORS.get(src, '#666666'), label=src)
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.legend(frameon=False, loc='lower left')
ax.grid(alpha=0.2)
savefig(FIG_DIR / 'Fig01_global_sample_map.png')

# Figure 2: richness distribution (x log scale)
r = df['amf_seq_richness'].dropna()
r = r[r > 0]
fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.hist(r, bins=50, color='#4c78a8', alpha=0.85, edgecolor='white')
ax.set_xscale('log')
ax.set_xlabel('AMF sequence richness (log scale)')
ax.set_ylabel('Sample count')
savefig(FIG_DIR / 'Fig02_richness_histogram.png')

# Figure 3: log reads distribution by source
fig, ax = plt.subplots(figsize=(7.2, 5.2))
for src, sub in df[['log_total_reads', 'source']].dropna().groupby('source'):
    ax.hist(sub['log_total_reads'], bins=45, alpha=0.5, color=COLORS.get(src, '#666666'), label=src)
ax.set_xlabel('log_total_reads')
ax.set_ylabel('Sample count')
ax.legend(frameon=False)
savefig(FIG_DIR / 'Fig03_log_reads_histogram.png')


# %% [cell 10] type=markdown
# ## Block 5: Figures 4-8 sample-level trait relationships


# %% [cell 11] type=code
sample_figs = [
    ('Root diameter', 'Fig04_RD_vs_richness_sample.png'),
    ('Root tissue density (RTD)', 'Fig05_RTD_vs_richness_sample.png'),
    ('Specific root area (SRA)', 'Fig06_SRA_vs_richness_sample.png'),
    ('Root N content', 'Fig07_RootN_vs_richness_sample.png'),
    ('Root P content', 'Fig08_RootP_vs_richness_sample.png'),
]
for trait, fname in sample_figs:
    sample_trait_plot(df, trait, 'log_richness', FIG_DIR / fname, xlabel=trait)


# %% [cell 12] type=markdown
# ## Block 6: Figures 9-13 species-level trait relationships


# %% [cell 13] type=code
species_figs = [
    ('Root diameter', 'Fig09_RD_vs_richness_species.png'),
    ('Root tissue density (RTD)', 'Fig10_RTD_vs_richness_species.png'),
    ('Specific root area (SRA)', 'Fig11_SRA_vs_richness_species.png'),
    ('Root N content', 'Fig12_RootN_vs_richness_species.png'),
    ('Root P content', 'Fig13_RootP_vs_richness_species.png'),
]
for trait, fname in species_figs:
    species_trait_plot(df, trait, FIG_DIR / fname, title=f'{trait} vs species mean richness')


# %% [cell 14] type=markdown
# ## Block 7: Figure 14 PCA of root economics space (genus level)


# %% [cell 15] type=code
pca_vars = ['Root diameter', 'Root tissue density (RTD)', 'Specific root area (SRA)', 'Root N content']
pca_data = df[['genus', 'woodiness_group'] + pca_vars].copy()
genus_agg = pca_data.groupby('genus', as_index=False).agg({
    'Root diameter': 'mean',
    'Root tissue density (RTD)': 'mean',
    'Specific root area (SRA)': 'mean',
    'Root N content': 'mean',
    'woodiness_group': lambda s: s.mode().iloc[0] if len(s.mode()) else s.iloc[0],
})
genus_agg = genus_agg.dropna(subset=pca_vars).copy().reset_index(drop=True)

X = genus_agg[pca_vars].to_numpy()
Xz = StandardScaler().fit_transform(X)
pca = PCA(n_components=2)
scores = pca.fit_transform(Xz)
load = pca.components_.T

fig, ax = plt.subplots(figsize=(8, 6))
for grp, sub_idx in genus_agg.groupby('woodiness_group').groups.items():
    idx = np.array(list(sub_idx))
    color = '#8c564b' if grp == 'woody' else '#17becf'
    ax.scatter(scores[idx, 0], scores[idx, 1], s=28, alpha=0.7, color=color, label=grp)

arrow_scale = 2.3
for i, v in enumerate(pca_vars):
    ax.arrow(0, 0, load[i, 0] * arrow_scale, load[i, 1] * arrow_scale, color='black', width=0.01, alpha=0.8)
    ax.text(load[i, 0] * arrow_scale * 1.08, load[i, 1] * arrow_scale * 1.08, v, fontsize=9)

ax.axhline(0, color='grey', linewidth=0.7)
ax.axvline(0, color='grey', linewidth=0.7)
ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)')
ax.legend(frameon=False)
savefig(FIG_DIR / 'Fig14_RES_PCA.png')
emit(f'PCA genera used: {len(genus_agg):,}')


# %% [cell 16] type=markdown
# ## Block 8: Figures 15+ woody vs nonwoody species-level relationships


# %% [cell 17] type=code
wood_fig_map = [
    ('Root diameter', 'woody', 'Fig15_RD_vs_richness_woody.png'),
    ('Root diameter', 'nonwoody', 'Fig16_RD_vs_richness_nonwoody.png'),
    ('Root tissue density (RTD)', 'woody', 'Fig17_RTD_vs_richness_woody.png'),
    ('Root tissue density (RTD)', 'nonwoody', 'Fig18_RTD_vs_richness_nonwoody.png'),
    ('Specific root area (SRA)', 'woody', 'Fig19_SRA_vs_richness_woody.png'),
    ('Specific root area (SRA)', 'nonwoody', 'Fig19b_SRA_vs_richness_nonwoody.png'),
    ('Root N content', 'woody', 'Fig19c_RootN_vs_richness_woody.png'),
    ('Root N content', 'nonwoody', 'Fig19d_RootN_vs_richness_nonwoody.png'),
    ('Root P content', 'woody', 'Fig19e_RootP_vs_richness_woody.png'),
    ('Root P content', 'nonwoody', 'Fig19f_RootP_vs_richness_nonwoody.png'),
]
for trait, grp, fname in wood_fig_map:
    species_trait_by_wood_plot(df, trait, grp, FIG_DIR / fname)


# %% [cell 18] type=markdown
# ## Block 9: Figures 20-21 climate vs richness


# %% [cell 19] type=code
sample_trait_plot(df, 'MAT', 'log_richness', FIG_DIR / 'Fig20_MAT_vs_richness.png', xlabel='MAT')
sample_trait_plot(df, 'MAP', 'log_richness', FIG_DIR / 'Fig21_MAP_vs_richness.png', xlabel='MAP')


# %% [cell 20] type=markdown
# ## Block 10: Statistical models (OLS HC3)


# %% [cell 21] type=code
def fit_trait_models(data, response, trait, interaction=False):
    if interaction:
        formula = f"{response} ~ Q('{trait}') * C(woodiness_group) + log_total_reads + MAT + MAP"
        needed = [response, trait, 'woodiness_group', 'log_total_reads', 'MAT', 'MAP']
    else:
        formula = f"{response} ~ Q('{trait}') + log_total_reads + MAT + MAP + C(source)"
        needed = [response, trait, 'log_total_reads', 'MAT', 'MAP', 'source']
    d = data[needed].dropna().copy()
    if len(d) < 30:
        return None, formula, len(d)
    model = ols(formula, data=d).fit(cov_type='HC3')
    return model, formula, len(d)

rows = []
for trait in trait_cols:
    for response, label in [('log_richness', 'seq_richness'), ('log_genus_richness', 'genus_richness')]:
        m, fml, n = fit_trait_models(df, response, trait, interaction=False)
        if m is None:
            rows.append({'model_type': label, 'trait': trait, 'formula': fml, 'N': n, 'beta': np.nan, 'se': np.nan, 'p': np.nan, 'R2': np.nan})
        else:
            coef_name = f"Q('{trait}')"
            rows.append({
                'model_type': label,
                'trait': trait,
                'formula': fml,
                'N': int(m.nobs),
                'beta': float(m.params.get(coef_name, np.nan)),
                'se': float(m.bse.get(coef_name, np.nan)),
                'p': float(m.pvalues.get(coef_name, np.nan)),
                'R2': float(m.rsquared),
            })

# Interaction models
for trait in trait_cols:
    m, fml, n = fit_trait_models(df, 'log_richness', trait, interaction=True)
    if m is None:
        rows.append({'model_type': 'interaction_woodiness', 'trait': trait, 'formula': fml, 'N': n, 'beta': np.nan, 'se': np.nan, 'p': np.nan, 'R2': np.nan})
    else:
        coef_name = f"Q('{trait}')"
        rows.append({
            'model_type': 'interaction_woodiness',
            'trait': trait,
            'formula': fml,
            'N': int(m.nobs),
            'beta': float(m.params.get(coef_name, np.nan)),
            'se': float(m.bse.get(coef_name, np.nan)),
            'p': float(m.pvalues.get(coef_name, np.nan)),
            'R2': float(m.rsquared),
        })

model_results = pd.DataFrame(rows)
emit(f'Model rows generated: {len(model_results):,}')
model_results


# %% [cell 22] type=markdown
# ## Block 11: Write markdown results summary


# %% [cell 23] type=code
trait_coverage = []
for t in trait_cols + ['Mycorrhiza_Fraction of root length or tips colonized']:
    trait_coverage.append({'trait': t, 'fraction_non_missing': float(df[t].notna().mean()) if t in df.columns else np.nan})
trait_cov_df = pd.DataFrame(trait_coverage)

# interpretation helpers
seq_models = model_results[model_results['model_type'].eq('seq_richness')].copy()
sig = seq_models[seq_models['p'] < 0.05].copy()
if len(sig):
    dir_lines = []
    for _, r in sig.iterrows():
        direction = 'positive' if r['beta'] > 0 else 'negative'
        dir_lines.append(f"- {r['trait']}: {direction} association (beta={r['beta']:.3g}, p={r['p']:.3g})")
else:
    dir_lines = ['- No trait effects reached p < 0.05 in primary richness models.']

lines = []
lines.append('# FRED-AMF25 Trait-Richness Results Summary')
lines.append('')
lines.append('## Dataset overview')
lines.append(f'- Total samples: {len(df):,}')
lines.append(f"- Samples per source: {df['source'].value_counts().to_dict()}")
lines.append(f"- Total species: {df['plant_species_canon'].nunique(dropna=True):,}")
lines.append(f"- Total genera: {df['genus'].nunique(dropna=True):,}")
lines.append('')

lines.append('## Trait coverage (fraction non-missing)')
lines.append('')
lines.append('| trait | fraction_non_missing |')
lines.append('|---|---:|')
for _, r in trait_cov_df.iterrows():
    lines.append(f"| {r['trait']} | {r['fraction_non_missing']:.3f} |")
lines.append('')

lines.append('## Model results')
lines.append('')
lines.append('| model_type | trait | beta | se | p | R2 | N |')
lines.append('|---|---|---:|---:|---:|---:|---:|')
for _, r in model_results.iterrows():
    b = '' if pd.isna(r['beta']) else f"{r['beta']:.4g}"
    se = '' if pd.isna(r['se']) else f"{r['se']:.4g}"
    p = '' if pd.isna(r['p']) else f"{r['p']:.4g}"
    r2 = '' if pd.isna(r['R2']) else f"{r['R2']:.4g}"
    n = '' if pd.isna(r['N']) else str(int(r['N']))
    lines.append(f"| {r['model_type']} | {r['trait']} | {b} | {se} | {p} | {r2} | {n} |")
lines.append('')

lines.append('## Interpretation relative to RES hypotheses')
lines.append('')
lines.extend(dir_lines)
lines.append('- These patterns are exploratory and not fully phylogenetically controlled; they should be interpreted as broad associations.')
lines.append('- Consistency checks across sample-level, species-level, and woody/nonwoody plots provide qualitative support where signs agree.')
lines.append('')

lines.append('## Printed output log')
lines.append('')
lines.append('```text')
lines.extend(LOG)
lines.append('```')

SUMMARY_FP.write_text('\n'.join(lines), encoding='utf-8')
emit(f'Wrote summary markdown: {SUMMARY_FP}')

print('OUTPUT_DATASET:', DATA_FP)
print('OUTPUT_NOTEBOOK:', Path('03_fred_amf25_trait_richness_analysis.ipynb').resolve())
print('OUTPUT_FIG_DIR:', FIG_DIR.resolve())
print('OUTPUT_SUMMARY:', SUMMARY_FP.resolve())
