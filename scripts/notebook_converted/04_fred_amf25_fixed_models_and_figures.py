# Auto-converted from 04_fred_amf25_fixed_models_and_figures.ipynb


# %% [cell 1] type=markdown
# # 04_fred_amf25_fixed_models_and_figures
# 
# Methodologically improved AMF-FRED analysis addressing pseudoreplication, sequencing-depth confounding,
# trait replication structure, PCA validity, and source imbalance.


# %% [cell 2] type=markdown
# ## Block 1: Imports, paths, global settings


# %% [cell 3] type=code

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

DATA_FP = Path('../Output/amf_traits_merged_sample_level.csv')
FIG_DIR = Path('../Output/Figures')
FIG_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_FP = Path('../Output/fred_amf25_results_summary_FIXED.md')
MODEL_FP = Path('../Output/fred_amf25_model_results_FIXED.csv')
MODEL_SPECIES_FP = Path('../Output/fred_amf25_model_results_specieslevel_FIXED.csv')

NOTEBOOK_FP = Path('04_fred_amf25_fixed_models_and_figures.ipynb')

assert DATA_FP.exists(), f'Missing dataset: {DATA_FP}'

plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

COLORS = {'GlobalAMFungi': '#1f77b4', 'EcoBank': '#d62728'}

LOG = []
FIG_FILES = []
WRITTEN_FILES = []

def emit(msg):
    s = str(msg)
    print(s)
    LOG.append(s)

def wrote(path):
    p = Path(path)
    WRITTEN_FILES.append(str(p))
    emit(f'Wrote: {p}')

def canon(x):
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    s = ' '.join(s.split())
    return s

def to_numeric_safe(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors='coerce')
    return df


# %% [cell 4] type=markdown
# ## Block 2: Load reference files and dataset QC


# %% [cell 5] type=code

# Reference files required by task
ref_rmd = Path('FRED.Rmd')
ref_docx = Path('FRED-AMF25 (1).docx')
assert ref_rmd.exists(), ref_rmd
assert ref_docx.exists(), ref_docx
emit(f'Reference found: {ref_rmd}')
emit(f'Reference found: {ref_docx}')

df = pd.read_csv(DATA_FP)

required_cols = [
    'amf_seq_richness', 'amf_genus_richness', 'total_reads', 'log_total_reads',
    'Root diameter', 'Root tissue density (RTD)', 'Specific root area (SRA)',
    'Root N content', 'Root P content',
    'Mycorrhiza_Fraction of root length or tips colonized',
    'Plant woodiness_TRY',
    'MAT', 'MAP',
    'sample_id', 'plant_species', 'plant_species_canon',
    'source', 'study_id', 'latitude', 'longitude'
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f'Missing required columns: {missing}')

trait_cols = [
    'Root diameter',
    'Root tissue density (RTD)',
    'Specific root area (SRA)',
    'Root N content',
    'Root P content',
]

numeric_cols = [
    'amf_seq_richness', 'amf_genus_richness', 'total_reads', 'log_total_reads',
    'MAT', 'MAP', 'latitude', 'longitude',
    'Mycorrhiza_Fraction of root length or tips colonized'
] + trait_cols

df = to_numeric_safe(df, numeric_cols)

# Impossible values handling
if (df['amf_seq_richness'].dropna() < 0).any() or (df['amf_genus_richness'].dropna() < 0).any():
    raise ValueError('Negative richness values detected; aborting as requested.')

df['total_reads'] = df['total_reads'].where(df['total_reads'] > 0, np.nan)
df['log_total_reads'] = np.log(df['total_reads'])

df['plant_species_canon'] = df['plant_species_canon'].map(canon)
df['genus'] = df['plant_species_canon'].fillna('').astype(str).str.split().str[0].replace({'': np.nan})

wood_raw = df['Plant woodiness_TRY'].fillna('').astype(str).str.lower()
df['woody_binary'] = np.where(wood_raw.str.contains('woody|tree|shrub', regex=True), 1, 0)
df['woodiness_group'] = np.where(df['woody_binary'] == 1, 'woody', 'nonwoody')

df['log_richness'] = np.log1p(df['amf_seq_richness'])
df['log_genus_richness'] = np.log1p(df['amf_genus_richness'])

qc_rows = len(df)
qc_unique_samples = df['sample_id'].nunique(dropna=True)
qc_unique_species = df['plant_species_canon'].nunique(dropna=True)
qc_unique_studies = df['study_id'].nunique(dropna=True)

emit(f'Rows: {qc_rows:,}')
emit(f'Unique samples: {qc_unique_samples:,}')
emit(f'Unique species: {qc_unique_species:,}')
emit(f'Unique studies: {qc_unique_studies:,}')
emit(f'Rows by source: {df["source"].value_counts(dropna=False).to_dict()}')

miss_cols = trait_cols + ['Mycorrhiza_Fraction of root length or tips colonized', 'MAT', 'MAP', 'log_total_reads']
missingness = pd.DataFrame({
    'column': miss_cols,
    'n_missing': [int(df[c].isna().sum()) for c in miss_cols],
    'frac_missing': [float(df[c].isna().mean()) for c in miss_cols],
})

reads_by_source = df.groupby('source')['log_total_reads'].describe()
missingness


# %% [cell 6] type=markdown
# ## Block 3: Depth confounding and residual outcomes


# %% [cell 7] type=code

depth_df = df[['log_richness', 'log_genus_richness', 'log_total_reads', 'source']].dropna().copy()

model_depth_seq = smf.ols('log_richness ~ log_total_reads + C(source)', data=depth_df).fit(cov_type='HC3')
model_depth_gen = smf.ols('log_genus_richness ~ log_total_reads + C(source)', data=depth_df).fit(cov_type='HC3')

model_depth_seq_nosrc = smf.ols('log_richness ~ log_total_reads', data=depth_df).fit(cov_type='HC3')
model_depth_gen_nosrc = smf.ols('log_genus_richness ~ log_total_reads', data=depth_df).fit(cov_type='HC3')

df['richness_resid_seq'] = np.nan
df['richness_resid_genus'] = np.nan
df['richness_resid_seq_depthonly'] = np.nan
df['richness_resid_genus_depthonly'] = np.nan

idx = depth_df.index
df.loc[idx, 'richness_resid_seq'] = model_depth_seq.resid
df.loc[idx, 'richness_resid_genus'] = model_depth_gen.resid
df.loc[idx, 'richness_resid_seq_depthonly'] = model_depth_seq_nosrc.resid
df.loc[idx, 'richness_resid_genus_depthonly'] = model_depth_gen_nosrc.resid

emit('Depth models fitted.')
emit(f'Seq depth model R2: {model_depth_seq.rsquared:.3f}')
emit(f'Genus depth model R2: {model_depth_gen.rsquared:.3f}')


# %% [cell 8] type=markdown
# ## Block 4: Model fitting helpers (MixedLM + robust fallback)


# %% [cell 9] type=code

def fit_mixed_or_fallback(formula, data):
    # primary attempt: random intercept for study_id + species variance component
    try:
        md = smf.mixedlm(
            formula=formula,
            data=data,
            groups=data['study_id'].astype(str),
            re_formula='1',
            vc_formula={'species': '0 + C(plant_species_canon)'}
        )
        res = md.fit(method='lbfgs', reml=False, maxiter=400, disp=False)
        return 'MixedLM', res, None
    except Exception as e1:
        # fallback: OLS with cluster robust SE by study_id
        try:
            ols_res = smf.ols(formula, data=data).fit(
                cov_type='cluster',
                cov_kwds={'groups': data['study_id'].astype(str)}
            )
            return 'OLS_cluster_study', ols_res, str(e1)
        except Exception as e2:
            return 'FAILED', None, f'{e1} | {e2}'


def extract_trait_result(res_type, res, trait_term, formula, outcome, trait, model_family, model_scope, n_used, fail_reason=None):
    base = {
        'model_family': model_family,
        'model_scope': model_scope,
        'outcome': outcome,
        'trait': trait,
        'formula': formula,
        'estimator': res_type,
        'N': int(n_used),
        'beta': np.nan,
        'se': np.nan,
        'p_value': np.nan,
        'ci_low': np.nan,
        'ci_high': np.nan,
        'aic_or_na': np.nan,
        'r2_or_na': np.nan,
        'var_study_or_na': np.nan,
        'var_species_or_na': np.nan,
        'note': fail_reason if fail_reason else ''
    }
    if res is None:
        return base

    term = trait_term
    if term in res.params.index:
        base['beta'] = float(res.params[term])
    if term in getattr(res, 'bse', pd.Series(dtype=float)).index:
        base['se'] = float(res.bse[term])
    if term in getattr(res, 'pvalues', pd.Series(dtype=float)).index:
        base['p_value'] = float(res.pvalues[term])

    try:
        ci = res.conf_int()
        if term in ci.index:
            base['ci_low'] = float(ci.loc[term, 0])
            base['ci_high'] = float(ci.loc[term, 1])
    except Exception:
        pass

    if hasattr(res, 'aic'):
        try:
            base['aic_or_na'] = float(res.aic)
        except Exception:
            pass
    if hasattr(res, 'rsquared'):
        try:
            base['r2_or_na'] = float(res.rsquared)
        except Exception:
            pass

    if res_type == 'MixedLM':
        try:
            base['var_study_or_na'] = float(np.diag(res.cov_re)[0])
        except Exception:
            pass
        try:
            if getattr(res, 'vcomp', None) is not None and len(res.vcomp) > 0:
                base['var_species_or_na'] = float(res.vcomp[0])
        except Exception:
            pass

    return base


# %% [cell 10] type=markdown
# ## Block 5: Primary sample-level models (pooled + per-source)


# %% [cell 11] type=code

sample_model_rows = []

primary_specs = [
    ('log_richness', 'primary_seq'),
    ('log_genus_richness', 'primary_genus'),
]

resid_specs = [
    ('richness_resid_seq', 'residual_seq'),
    ('richness_resid_genus', 'residual_genus'),
]

for trait in trait_cols:
    trait_term = f"Q('{trait}')"

    # pooled primary models
    for outcome, fam in primary_specs:
        formula = f"{outcome} ~ {trait_term} + log_total_reads + MAT + MAP + C(source) + woody_binary"
        need = [outcome, trait, 'log_total_reads', 'MAT', 'MAP', 'source', 'woody_binary', 'study_id', 'plant_species_canon']
        d = df[need].dropna().copy()
        est, res, reason = fit_mixed_or_fallback(formula, d)
        row = extract_trait_result(est, res, trait_term, formula, outcome, trait, fam, 'pooled', len(d), reason)
        sample_model_rows.append(row)

    # pooled interaction model
    formula_i = f"log_richness ~ {trait_term} * woody_binary + log_total_reads + MAT + MAP + C(source)"
    need_i = ['log_richness', trait, 'woody_binary', 'log_total_reads', 'MAT', 'MAP', 'source', 'study_id', 'plant_species_canon']
    d_i = df[need_i].dropna().copy()
    est_i, res_i, reason_i = fit_mixed_or_fallback(formula_i, d_i)
    row_i = extract_trait_result(est_i, res_i, trait_term, formula_i, 'log_richness', trait, 'interaction_seq', 'pooled', len(d_i), reason_i)
    sample_model_rows.append(row_i)

    # pooled residual sensitivity
    for outcome, fam in resid_specs:
        formula_r = f"{outcome} ~ {trait_term} + MAT + MAP + C(source) + woody_binary"
        need_r = [outcome, trait, 'MAT', 'MAP', 'source', 'woody_binary', 'study_id', 'plant_species_canon']
        d_r = df[need_r].dropna().copy()
        est_r, res_r, reason_r = fit_mixed_or_fallback(formula_r, d_r)
        row_r = extract_trait_result(est_r, res_r, trait_term, formula_r, outcome, trait, fam, 'pooled', len(d_r), reason_r)
        sample_model_rows.append(row_r)

    # stratified by source (primary seq only for concise robustness)
    for src in sorted(df['source'].dropna().unique()):
        ds = df[df['source'] == src].copy()
        formula_s = f"log_richness ~ {trait_term} + log_total_reads + MAT + MAP + woody_binary"
        need_s = ['log_richness', trait, 'log_total_reads', 'MAT', 'MAP', 'woody_binary', 'study_id', 'plant_species_canon']
        ds = ds[need_s].dropna().copy()
        est_s, res_s, reason_s = fit_mixed_or_fallback(formula_s, ds)
        row_s = extract_trait_result(est_s, res_s, trait_term, formula_s, 'log_richness', trait, 'primary_seq', f'source={src}', len(ds), reason_s)
        sample_model_rows.append(row_s)

sample_model_results = pd.DataFrame(sample_model_rows)
sample_model_results.to_csv(MODEL_FP, index=False)
wrote(MODEL_FP)
emit(f'Sample-level model rows: {len(sample_model_results):,}')
sample_model_results.head(10)


# %% [cell 12] type=markdown
# ## Block 6: Species-level weighted models


# %% [cell 13] type=code

species_agg = (
    df.groupby('plant_species_canon', dropna=True)
      .agg(
          species_mean_log_richness=('log_richness', 'mean'),
          species_mean_log_genus_richness=('log_genus_richness', 'mean'),
          species_mean_resid_seq=('richness_resid_seq', 'mean'),
          species_mean_resid_genus=('richness_resid_genus', 'mean'),
          mean_log_total_reads=('log_total_reads', 'mean'),
          mean_MAT=('MAT', 'mean'),
          mean_MAP=('MAP', 'mean'),
          woody_binary=('woody_binary', lambda s: int(np.nanmean(s) >= 0.5) if s.notna().any() else np.nan),
          n_samples=('sample_id', 'nunique'),
          **{f'mean_{t}': (t, 'mean') for t in trait_cols}
      )
      .reset_index()
)

species_src_agg = (
    df.groupby(['plant_species_canon', 'source'], dropna=True)
      .agg(
          species_mean_log_richness=('log_richness', 'mean'),
          species_mean_log_genus_richness=('log_genus_richness', 'mean'),
          mean_log_total_reads=('log_total_reads', 'mean'),
          mean_MAT=('MAT', 'mean'),
          mean_MAP=('MAP', 'mean'),
          woody_binary=('woody_binary', lambda s: int(np.nanmean(s) >= 0.5) if s.notna().any() else np.nan),
          n_samples=('sample_id', 'nunique'),
          **{f'mean_{t}': (t, 'mean') for t in trait_cols}
      )
      .reset_index()
)

species_rows = []

def fit_wls(formula, data, weights):
    m = smf.wls(formula, data=data, weights=weights).fit(cov_type='HC3')
    return m

for trait in trait_cols:
    tcol = f'mean_{trait}'
    term = f"Q('{tcol}')"

    for out, fam in [('species_mean_log_richness', 'species_seq'), ('species_mean_log_genus_richness', 'species_genus')]:
        form = f"{out} ~ {term} + mean_log_total_reads + mean_MAT + mean_MAP + woody_binary"
        use_cols = [out, tcol, 'mean_log_total_reads', 'mean_MAT', 'mean_MAP', 'woody_binary', 'n_samples']
        d = species_agg[use_cols].dropna().copy()
        if len(d) < 20:
            species_rows.append({'model_family': fam, 'scope': 'pooled', 'trait': trait, 'outcome': out, 'N': len(d), 'beta': np.nan, 'se': np.nan, 'p_value': np.nan, 'ci_low': np.nan, 'ci_high': np.nan, 'r2': np.nan})
        else:
            res = fit_wls(form, d, d['n_samples'])
            ci = res.conf_int().loc[term] if term in res.params.index else [np.nan, np.nan]
            species_rows.append({'model_family': fam, 'scope': 'pooled', 'trait': trait, 'outcome': out, 'N': int(res.nobs), 'beta': float(res.params.get(term, np.nan)), 'se': float(res.bse.get(term, np.nan)), 'p_value': float(res.pvalues.get(term, np.nan)), 'ci_low': float(ci[0]), 'ci_high': float(ci[1]), 'r2': float(res.rsquared)})

        # per-source species-level
        for src in sorted(species_src_agg['source'].dropna().unique()):
            ds = species_src_agg[species_src_agg['source'] == src].copy()
            use_cols_s = [out, tcol, 'mean_log_total_reads', 'mean_MAT', 'mean_MAP', 'woody_binary', 'n_samples']
            ds = ds[use_cols_s].dropna().copy()
            if len(ds) < 15:
                species_rows.append({'model_family': fam, 'scope': f'source={src}', 'trait': trait, 'outcome': out, 'N': len(ds), 'beta': np.nan, 'se': np.nan, 'p_value': np.nan, 'ci_low': np.nan, 'ci_high': np.nan, 'r2': np.nan})
            else:
                res_s = fit_wls(form, ds, ds['n_samples'])
                ci_s = res_s.conf_int().loc[term] if term in res_s.params.index else [np.nan, np.nan]
                species_rows.append({'model_family': fam, 'scope': f'source={src}', 'trait': trait, 'outcome': out, 'N': int(res_s.nobs), 'beta': float(res_s.params.get(term, np.nan)), 'se': float(res_s.bse.get(term, np.nan)), 'p_value': float(res_s.pvalues.get(term, np.nan)), 'ci_low': float(ci_s[0]), 'ci_high': float(ci_s[1]), 'r2': float(res_s.rsquared)})

species_model_results = pd.DataFrame(species_rows)
species_model_results.to_csv(MODEL_SPECIES_FP, index=False)
wrote(MODEL_SPECIES_FP)
emit(f'Species-level model rows: {len(species_model_results):,}')
species_model_results.head(10)


# %% [cell 14] type=markdown
# ## Block 7: Figure helpers


# %% [cell 15] type=code

def savefig(fp):
    plt.tight_layout()
    plt.savefig(fp, dpi=300, bbox_inches='tight')
    plt.close()
    FIG_FILES.append(str(fp))
    emit(f'Wrote figure: {fp}')

def scatter_lowess_by_source(data, x, y, fp, xlabel=None, ylabel=None):
    d = data[[x, y, 'source']].dropna().copy()
    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    for src, sub in d.groupby('source'):
        c = COLORS.get(src, '#666666')
        ax.scatter(sub[x], sub[y], s=14, alpha=0.35, color=c, label=src)
        if len(sub) >= 25:
            sm = lowess(sub[y].to_numpy(), sub[x].to_numpy(), frac=0.35, return_sorted=True)
            ax.plot(sm[:, 0], sm[:, 1], color=c, linewidth=2)
    ax.set_xlabel(xlabel if xlabel else x)
    ax.set_ylabel(ylabel if ylabel else y)
    ax.legend(frameon=False)
    savefig(fp)

def species_scatter(data, trait, ycol, fp):
    d = data[['plant_species_canon', trait, ycol, 'woody_binary', 'sample_id']].dropna().copy()
    agg = d.groupby('plant_species_canon', as_index=False).agg(
        trait_mean=(trait, 'mean'),
        y_mean=(ycol, 'mean'),
        n_samples=('sample_id', 'nunique'),
        woody_binary=('woody_binary', lambda s: int(np.nanmean(s) >= 0.5) if s.notna().any() else 0)
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    sizes = 14 + 3 * np.sqrt(agg['n_samples'].to_numpy())
    colors = np.where(agg['woody_binary'].to_numpy() == 1, '#8c564b', '#17becf')
    ax.scatter(agg['trait_mean'], agg['y_mean'], s=sizes, c=colors, alpha=0.55, edgecolor='none')
    if len(agg) >= 20:
        sm = lowess(agg['y_mean'].to_numpy(), agg['trait_mean'].to_numpy(), frac=0.45, return_sorted=True)
        ax.plot(sm[:, 0], sm[:, 1], color='black', linewidth=2)
    ax.set_xlabel(trait)
    ax.set_ylabel('Species mean residual richness (seq)')
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', label='nonwoody', markerfacecolor='#17becf', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='woody', markerfacecolor='#8c564b', markersize=8)
    ]
    ax.legend(handles=legend_elements, frameon=False)
    savefig(fp)

def woody_nonwoody_panel(data, trait, fp):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    notes = []
    for i, grp in enumerate(['woody', 'nonwoody']):
        d = data[data['woodiness_group'] == grp][['plant_species_canon', trait, 'richness_resid_seq', 'sample_id']].dropna().copy()
        agg = d.groupby('plant_species_canon', as_index=False).agg(
            trait_mean=(trait, 'mean'),
            y_mean=('richness_resid_seq', 'mean'),
            n_samples=('sample_id', 'nunique')
        )
        ax = axes[i]
        color = '#8c564b' if grp == 'woody' else '#17becf'
        if len(agg) < 5:
            ax.text(0.5, 0.5, f'{grp}: <5 species', ha='center', va='center', transform=ax.transAxes)
            notes.append(f'{trait}: {grp} had <5 species; panel not plotted.')
        else:
            sz = 16 + 3 * np.sqrt(agg['n_samples'].to_numpy())
            ax.scatter(agg['trait_mean'], agg['y_mean'], s=sz, color=color, alpha=0.55, edgecolor='none')
            if len(agg) >= 10:
                sm = lowess(agg['y_mean'].to_numpy(), agg['trait_mean'].to_numpy(), frac=0.5, return_sorted=True)
                ax.plot(sm[:, 0], sm[:, 1], color='black', linewidth=2)
            ax.set_xlabel(trait)
        ax.set_title(grp)
    axes[0].set_ylabel('Species mean residual richness (seq)')
    savefig(fp)
    return notes


# %% [cell 16] type=markdown
# ## Block 8: Generate _FIXED figures


# %% [cell 17] type=code

# Fig01 global sample map
map_df = df[['longitude', 'latitude', 'source']].dropna().copy()
fig, ax = plt.subplots(figsize=(10, 5.2))
for src, sub in map_df.groupby('source'):
    ax.scatter(sub['longitude'], sub['latitude'], s=10, alpha=0.45, color=COLORS.get(src, '#666666'), label=src)
ax.set_xlim(-180, 180)
ax.set_ylim(-90, 90)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.legend(frameon=False)
ax.grid(alpha=0.2)
savefig(FIG_DIR / 'Fig01_global_sample_map_FIXED.png')

# Fig02 richness distribution
rr = df['amf_seq_richness'].dropna()
rr = rr[rr > 0]
fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.hist(rr, bins=55, color='#4c78a8', alpha=0.9, edgecolor='white')
ax.set_xscale('log')
ax.set_xlabel('AMF sequence richness (log x-scale)')
ax.set_ylabel('Samples')
savefig(FIG_DIR / 'Fig02_richness_histogram_FIXED.png')

# Fig03 log reads by source
fig, ax = plt.subplots(figsize=(7.2, 5.2))
for src, sub in df[['log_total_reads', 'source']].dropna().groupby('source'):
    ax.hist(sub['log_total_reads'], bins=45, alpha=0.5, color=COLORS.get(src, '#666666'), label=src)
ax.set_xlabel('log_total_reads')
ax.set_ylabel('Samples')
ax.legend(frameon=False)
savefig(FIG_DIR / 'Fig03_log_reads_histogram_by_source_FIXED.png')

# Trait -> residual richness, sample-level
sample_figs = [
    ('Root diameter', 'Fig04_RD_vs_richness_sample_FIXED.png'),
    ('Root tissue density (RTD)', 'Fig05_RTD_vs_richness_sample_FIXED.png'),
    ('Specific root area (SRA)', 'Fig06_SRA_vs_richness_sample_FIXED.png'),
    ('Root N content', 'Fig07_RootN_vs_richness_sample_FIXED.png'),
    ('Root P content', 'Fig08_RootP_vs_richness_sample_FIXED.png'),
]
for trait, fname in sample_figs:
    scatter_lowess_by_source(df, trait, 'richness_resid_seq', FIG_DIR / fname, xlabel=trait, ylabel='Residual richness (seq)')

# Trait -> residual richness, species-level
species_figs = [
    ('Root diameter', 'Fig09_RD_vs_richness_species_FIXED.png'),
    ('Root tissue density (RTD)', 'Fig10_RTD_vs_richness_species_FIXED.png'),
    ('Specific root area (SRA)', 'Fig11_SRA_vs_richness_species_FIXED.png'),
    ('Root N content', 'Fig12_RootN_vs_richness_species_FIXED.png'),
    ('Root P content', 'Fig13_RootP_vs_richness_species_FIXED.png'),
]
for trait, fname in species_figs:
    species_scatter(df, trait, 'richness_resid_seq', FIG_DIR / fname)

# Woody/nonwoody panels (minimum requested)
panel_notes = []
panel_notes.extend(woody_nonwoody_panel(df, 'Root diameter', FIG_DIR / 'Fig15_RD_woody_nonwoody_species_FIXED.png'))
panel_notes.extend(woody_nonwoody_panel(df, 'Root N content', FIG_DIR / 'Fig16_RootN_woody_nonwoody_species_FIXED.png'))

# Climate vs residual richness
scatter_lowess_by_source(df, 'MAT', 'richness_resid_seq', FIG_DIR / 'Fig20_MAT_vs_richness_resid_FIXED.png', ylabel='Residual richness (seq)')
scatter_lowess_by_source(df, 'MAP', 'richness_resid_seq', FIG_DIR / 'Fig21_MAP_vs_richness_resid_FIXED.png', ylabel='Residual richness (seq)')

emit(f'Generated fixed figures: {len(FIG_FILES):,}')


# %% [cell 18] type=markdown
# ## Block 9: PCA redesign (2-trait + conditional 4-trait)


# %% [cell 19] type=code

pca_notes = []

def run_and_plot_pca(table, vars_, fp, title, min_n=20):
    t = table.dropna(subset=vars_).copy()
    n_taxa = len(t)
    if n_taxa < min_n:
        return None, n_taxa
    X = StandardScaler().fit_transform(t[vars_].to_numpy())
    pca = PCA(n_components=2)
    scores = pca.fit_transform(X)
    load = pca.components_.T

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = np.where(t['woody_binary'].to_numpy() == 1, '#8c564b', '#17becf')
    ax.scatter(scores[:, 0], scores[:, 1], s=26, alpha=0.65, c=colors)

    scale = 2.2
    for i, v in enumerate(vars_):
        ax.arrow(0, 0, load[i, 0] * scale, load[i, 1] * scale, color='black', width=0.01, alpha=0.8)
        ax.text(load[i, 0] * scale * 1.08, load[i, 1] * scale * 1.08, v, fontsize=9)

    ax.axhline(0, color='grey', linewidth=0.7)
    ax.axvline(0, color='grey', linewidth=0.7)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)')
    ax.set_title(f'{title} (N taxa={n_taxa})')
    savefig(fp)
    return pca, n_taxa

# Taxa-level tables for PCA coverage comparison
species_pca_tbl = (
    df.groupby('plant_species_canon', dropna=True)
      .agg(
          woody_binary=('woody_binary', lambda s: int(np.nanmean(s) >= 0.5) if s.notna().any() else 0),
          **{t: (t, 'mean') for t in ['Root diameter', 'Root N content', 'Root tissue density (RTD)', 'Specific root area (SRA)']}
      )
      .reset_index()
)

genus_pca_tbl = (
    df.groupby('genus', dropna=True)
      .agg(
          woody_binary=('woody_binary', lambda s: int(np.nanmean(s) >= 0.5) if s.notna().any() else 0),
          **{t: (t, 'mean') for t in ['Root diameter', 'Root N content', 'Root tissue density (RTD)', 'Specific root area (SRA)']}
      )
      .reset_index()
)

# PCA 2-trait chooses taxa level with better complete-case coverage
p2_vars = ['Root diameter', 'Root N content']
n_species_2 = species_pca_tbl.dropna(subset=p2_vars).shape[0]
n_genus_2 = genus_pca_tbl.dropna(subset=p2_vars).shape[0]

if n_species_2 >= n_genus_2:
    p2_table = species_pca_tbl
    p2_level = 'species'
else:
    p2_table = genus_pca_tbl
    p2_level = 'genus'

p2, n2 = run_and_plot_pca(p2_table, p2_vars, FIG_DIR / 'Fig14_PCA_2trait_FIXED.png', f'PCA 2-trait ({p2_level})', min_n=5)
if p2 is None:
    pca_notes.append(f'PCA-2trait skipped: insufficient taxa (N={n2}, minimum=5).')
else:
    pca_notes.append(f'PCA-2trait completed on {p2_level} level with N={n2}.')

# PCA RES 4-trait (conditional)
p4_vars = ['Root diameter', 'Root tissue density (RTD)', 'Specific root area (SRA)', 'Root N content']
n_species_4 = species_pca_tbl.dropna(subset=p4_vars).shape[0]
n_genus_4 = genus_pca_tbl.dropna(subset=p4_vars).shape[0]

if max(n_species_4, n_genus_4) < 20:
    pca_notes.append(f'PCA-RES4 skipped due to insufficient coverage (species N={n_species_4}, genus N={n_genus_4}).')
else:
    if n_species_4 >= n_genus_4:
        p4_table = species_pca_tbl
        p4_level = 'species'
    else:
        p4_table = genus_pca_tbl
        p4_level = 'genus'
    p4, n4 = run_and_plot_pca(p4_table, p4_vars, FIG_DIR / 'Fig14b_PCA_RES4trait_FIXED.png', f'PCA RES 4-trait ({p4_level})')
    if p4 is None:
        pca_notes.append(f'PCA-RES4 skipped at execution stage (N={n4}).')
    else:
        pca_notes.append(f'PCA-RES4 completed on {p4_level} level with N={n4}.')

for note in pca_notes:
    emit(note)

for note in panel_notes:
    emit(note)


# %% [cell 20] type=markdown
# ## Block 10: Build improved markdown summary


# %% [cell 21] type=code

def md_table_from_df(df_in, cols=None, float_cols=None, max_rows=None):
    d = df_in.copy()
    if cols is not None:
        d = d[cols]
    if max_rows is not None:
        d = d.head(max_rows)
    if float_cols:
        for c in float_cols:
            if c in d.columns:
                d[c] = d[c].map(lambda x: '' if pd.isna(x) else f'{x:.4g}')
    return d

coverage_by_source = (
    df.groupby('source')[trait_cols + ['Mycorrhiza_Fraction of root length or tips colonized', 'MAT', 'MAP', 'log_total_reads']]
      .apply(lambda x: x.notna().mean())
      .reset_index()
)

# concise finding extraction
primary_pooled = sample_model_results[(sample_model_results['model_family'] == 'primary_seq') & (sample_model_results['model_scope'] == 'pooled')]
species_pooled = species_model_results[(species_model_results['model_family'] == 'species_seq') & (species_model_results['scope'] == 'pooled')]

def find_trait_row(df_in, trait):
    d = df_in[df_in['trait'] == trait]
    return d.iloc[0] if len(d) else None

rd_sample = find_trait_row(primary_pooled, 'Root diameter')
rn_sample = find_trait_row(primary_pooled, 'Root N content')
rd_species = find_trait_row(species_pooled, 'Root diameter')
rn_species = find_trait_row(species_pooled, 'Root N content')

lines = []
lines.append('# FRED-AMF25 Results Summary (FIXED)')
lines.append('')
lines.append('## Why prior OLS was problematic and what was fixed')
lines.append('- Prior OLS treated sample rows as independent despite clustering by `study_id` and repeated species-level traits, inflating nominal significance.')
lines.append('- This analysis uses mixed-effects models with random intercept structure (study + species variance component) where feasible, with cluster-robust fallback when mixed fitting fails.')
lines.append('- Species-level weighted models (weights=`n_samples`) provide a second inference layer aligned with species-level trait measurement.')
lines.append('- Richness was residualized against sequencing depth (`log_total_reads`) and source to reduce confounding.')
lines.append('')

lines.append('## QC and coverage')
lines.append(f'- Rows: {qc_rows:,}')
lines.append(f'- Unique samples: {qc_unique_samples:,}')
lines.append(f'- Unique species: {qc_unique_species:,}')
lines.append(f'- Unique studies: {qc_unique_studies:,}')
lines.append(f"- Rows by source: {df['source'].value_counts(dropna=False).to_dict()}")
lines.append('')
lines.append('### Missingness (overall)')
lines.append('| variable | n_missing | frac_missing |')
lines.append('|---|---:|---:|')
for _, r in missingness.iterrows():
    lines.append(f"| {r['column']} | {int(r['n_missing'])} | {r['frac_missing']:.3f} |")
lines.append('')

lines.append('### Coverage by source (fraction non-missing)')
lines.append('| source | variable | frac_non_missing |')
lines.append('|---|---|---:|')
for _, r in coverage_by_source.iterrows():
    src = r['source']
    for c in [x for x in coverage_by_source.columns if x != 'source']:
        lines.append(f'| {src} | {c} | {r[c]:.3f} |')
lines.append('')

lines.append('## Sequencing depth confounding')
lines.append(f'- Depth model (seq richness) R2: {model_depth_seq.rsquared:.3f}')
lines.append(f'- Depth model (genus richness) R2: {model_depth_gen.rsquared:.3f}')
lines.append('- Residual outcomes (`richness_resid_seq`, `richness_resid_genus`) were used in sensitivity models and figures to reduce read-depth bias.')
lines.append('')

lines.append('## Main results (Root diameter and Root N)')
if rd_sample is not None:
    lines.append(f"- Sample-level mixed/robust (pooled), Root diameter -> seq richness: beta={rd_sample['beta']:.4g}, p={rd_sample['p_value']:.4g}, N={int(rd_sample['N'])}")
if rn_sample is not None:
    lines.append(f"- Sample-level mixed/robust (pooled), Root N -> seq richness: beta={rn_sample['beta']:.4g}, p={rn_sample['p_value']:.4g}, N={int(rn_sample['N'])}")
if rd_species is not None:
    lines.append(f"- Species-level weighted (pooled), Root diameter -> seq richness: beta={rd_species['beta']:.4g}, p={rd_species['p_value']:.4g}, N={int(rd_species['N'])}")
if rn_species is not None:
    lines.append(f"- Species-level weighted (pooled), Root N -> seq richness: beta={rn_species['beta']:.4g}, p={rn_species['p_value']:.4g}, N={int(rn_species['N'])}")
lines.append('- Persistence across pooled mixed/robust models, species-level weighted models, and source-stratified models should be interpreted as stronger evidence than any single model alone.')
lines.append('')

lines.append('## Woodiness interactions')
inter = sample_model_results[(sample_model_results['model_family'] == 'interaction_seq') & (sample_model_results['model_scope'] == 'pooled')]
if len(inter):
    lines.append('| trait | beta_main_trait | p_value | N |')
    lines.append('|---|---:|---:|---:|')
    for _, r in inter.iterrows():
        b = '' if pd.isna(r['beta']) else f"{r['beta']:.4g}"
        p = '' if pd.isna(r['p_value']) else f"{r['p_value']:.4g}"
        lines.append(f"| {r['trait']} | {b} | {p} | {int(r['N'])} |")
else:
    lines.append('- No interaction model estimates available.')
lines.append('')

lines.append('## PCA validity')
for note in pca_notes:
    lines.append(f'- {note}')
lines.append('')

lines.append('## Model result files')
lines.append(f'- `{MODEL_FP}`')
lines.append(f'- `{MODEL_SPECIES_FP}`')
lines.append('')

lines.append('## Figure inventory (_FIXED)')
for fp in FIG_FILES:
    lines.append(f'- `{fp}`')
lines.append('')

lines.append('## Limitations')
lines.append('- Trait missingness remains high for several RES dimensions, limiting precision and some multivariate analyses.')
lines.append('- Species-level traits are reused across samples; mixed-effects and species-level analyses reduce, but do not eliminate, dependence concerns.')
lines.append('- Plant phylogeny is not explicitly modeled; trait effects may partially proxy phylogenetic structure.')
lines.append('- EcoBank coverage and sequencing characteristics differ from GlobalAMFungi; source-stratified results are therefore critical context.')
lines.append('')

lines.append('## Notebook log')
lines.append('```text')
lines.extend(LOG)
lines.append('```')

SUMMARY_FP.write_text('\n'.join(lines), encoding='utf-8')
wrote(SUMMARY_FP)

emit('Final generated files:')
for f in [NOTEBOOK_FP, SUMMARY_FP, MODEL_FP, MODEL_SPECIES_FP] + [Path(x) for x in FIG_FILES]:
    emit(f'Wrote: {f}')
