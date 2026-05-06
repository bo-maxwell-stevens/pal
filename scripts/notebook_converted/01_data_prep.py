# Auto-converted from 01_data_prep.ipynb


# %% [cell 1] type=code
# --- Block 1: Imports + paths (run once) ---
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import geopandas as gpd

DATA_DIR = Path("../Data")
OUT_DIR  = DATA_DIR  # (R saved into the same folder; change if you want)

ECObank_FP = DATA_DIR / "ecobank_full.csv"
VT_FP      = DATA_DIR / "ecobank_cultured.csv"
FRED_FP    = DATA_DIR / "FRED3_fineRoots.csv"

assert ECObank_FP.exists(), ECObank_FP
assert VT_FP.exists(), VT_FP
assert FRED_FP.exists(), FRED_FP


# %% [cell 2] type=code
# --- Block 2: Load data (matches R read.csv + separators) ---
ecobank = pd.read_csv(ECObank_FP, sep=";")
vt      = pd.read_csv(VT_FP, sep=";")

# --- Block: Load FRED properly (skip first row, use second row as header) ---

# Read raw file without assigning header
fred_raw = pd.read_csv(FRED_FP, header=None)

# Use the SECOND row (index 1) as column names
fred_raw.columns = fred_raw.iloc[1].astype(str)

# Drop the first two rows (original row 1 + header row)
fred = fred_raw.iloc[2:].copy()

# Reset index
fred.reset_index(drop=True, inplace=True)

print("FRED shape:", fred.shape)
print("First columns:", fred.columns[:10])


# %% [cell 3] type=code
# --- Block 3: Sort + merge + filter roots (matches R order/merge/subset) ---
ecobank = ecobank.sort_values("sample")
vt      = vt.sort_values("sample")

roots = ecobank.merge(vt, on="sample", how="inner")

# R: roots <- subset(roots, isRoot == TRUE)
# In some exports this may be True/False, 0/1, or "TRUE"/"FALSE".
if roots["isRoot"].dtype == object:
    is_root = roots["isRoot"].astype(str).str.upper().isin(["TRUE", "T", "1", "YES"])
else:
    is_root = roots["isRoot"].astype(bool)

roots = roots.loc[is_root].copy()

print("roots:", roots.shape)


# %% [cell 4] type=code
# --- Block 4: Unique plant_species + write CSV (matches R unique + write.csv) ---
plant_species = (
    pd.Series(roots["plant_species"].dropna().unique(), name="plant_species")
    .sort_values()
    .reset_index(drop=True)
)

out_species_fp = DATA_DIR / "plant_species_Tartu.csv"
plant_species.to_csv(out_species_fp, index=False)

out_species_fp, plant_species.shape


# %% [cell 5] type=code
# --- Block: Load custom Natural Earth boundary shapefile ---
MAP_FP = Path("../../darkdivnet/Manuscript/Data/Maps/110m_cultural/ne_110m_admin_0_boundary_lines_land.shp")

assert MAP_FP.exists(), MAP_FP

basemap = gpd.read_file(MAP_FP).to_crs("EPSG:4326")

print(basemap.head())


# %% [cell 6] type=code
# --- Block 5: Basemap (map_units) + FRED points ---

MAP_DIR = Path("../../darkdivnet/Manuscript/Data/Maps/110m_cultural")

MAP_UNITS_FP = MAP_DIR / "ne_110m_admin_0_map_units.shp"
BORDERS_FP   = MAP_DIR / "ne_110m_admin_0_boundary_lines_land.shp"

assert MAP_UNITS_FP.exists(), MAP_UNITS_FP
assert BORDERS_FP.exists(), BORDERS_FP

# Basemap polygons (matches R: ne_countries(type="map_units"))
map_units = gpd.read_file(MAP_UNITS_FP).to_crs("EPSG:4326")

# Optional overlay: boundary linework for sharper borders
borders = gpd.read_file(BORDERS_FP).to_crs("EPSG:4326")

# Coerce FRED coordinates to numeric
fred["lat"] = pd.to_numeric(fred["Latitude_Main"], errors="coerce")
fred["lon"] = pd.to_numeric(fred["Longitude_Main"], errors="coerce")

mask_fred = np.isfinite(fred["lat"].to_numpy()) & np.isfinite(fred["lon"].to_numpy())
fred_xy = fred.loc[mask_fred, ["lat", "lon"]].copy()

print(f"FRED total rows: {len(fred):,}")
print(f"FRED rows with finite Latitude_Main/Longitude_Main: {len(fred_xy):,}")

fred_pts = gpd.GeoDataFrame(
    fred_xy,
    geometry=gpd.points_from_xy(fred_xy["lon"], fred_xy["lat"]),
    crs="EPSG:4326",
)

# Plot
fig, ax = plt.subplots(figsize=(11.81, 5.91))

# Polygons (no fill, just outlines)
map_units.boundary.plot(ax=ax, linewidth=0.4, color="black")

# Optional crisp borders overlay (can comment out if redundant)
borders.plot(ax=ax, linewidth=0.3, color="black")

# Avoid geopandas auto-aspect computation
ax.set_aspect("equal", adjustable="box")

if len(fred_pts) == 0:
    ax.text(
        0.5, 0.5,
        "No FRED points with finite Latitude_Main/Longitude_Main",
        ha="center", va="center",
        transform=ax.transAxes
    )
else:
    fred_pts.plot(ax=ax, markersize=5)

ax.set_axis_off()
plt.tight_layout()

out_map_fred = DATA_DIR / "map_fred.jpg"
plt.savefig(out_map_fred, dpi=300)
plt.show()

out_map_fred


# %% [cell 7] type=code
# --- Block 6: Basemap (map_units) + roots + FRED points ---

# Coerce roots coordinates to numeric
roots["lat"] = pd.to_numeric(roots["lat"], errors="coerce")
roots["lon"] = pd.to_numeric(roots["lon"], errors="coerce")

mask_roots = np.isfinite(roots["lat"].to_numpy()) & np.isfinite(roots["lon"].to_numpy())
roots_xy = roots.loc[mask_roots, ["lat", "lon"]].copy()

print(f"roots total rows: {len(roots):,}")
print(f"roots rows with finite lat/lon: {len(roots_xy):,}")

roots_pts = gpd.GeoDataFrame(
    roots_xy,
    geometry=gpd.points_from_xy(roots_xy["lon"], roots_xy["lat"]),
    crs="EPSG:4326",
)

# Plot
fig, ax = plt.subplots(figsize=(11.81, 5.91))

map_units.boundary.plot(ax=ax, linewidth=0.4, color="black")
borders.plot(ax=ax, linewidth=0.3, color="black")

ax.set_aspect("equal", adjustable="box")

if len(roots_pts) == 0 and len(fred_pts) == 0:
    ax.text(
        0.5, 0.5,
        "No roots or FRED points with finite coordinates",
        ha="center", va="center",
        transform=ax.transAxes
    )
else:
    if len(roots_pts) > 0:
        roots_pts.plot(ax=ax, markersize=5)
    if len(fred_pts) > 0:
        fred_pts.plot(ax=ax, markersize=5)

ax.set_axis_off()
plt.tight_layout()

out_map_both = DATA_DIR / "map_both.jpg"
plt.savefig(out_map_both, dpi=300)
plt.show()

out_map_both


# %% [cell 8] type=code
# --- Block 7: Build robust species names in FRED + overlap with roots ---

# Prefer accepted taxonomy if present; fallback to data-source taxonomy
genus_candidates = [
    "Plant taxonomy_Accepted genus_TPL",
    "Plant taxonomy_Genus_Data Source",
]
species_candidates = [
    "Plant Taxonomy_Accepted species_TPL",
    "Plant taxonomy_Species_Data source",
]

def first_existing_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of these columns found: {candidates}")

genus_col   = first_existing_col(fred, genus_candidates)
species_col = first_existing_col(fred, species_candidates)

print("Using genus column:  ", genus_col)
print("Using species column:", species_col)

# Clean and build binomial name
genus_clean = fred[genus_col].astype(str).str.strip()
species_clean = fred[species_col].astype(str).str.strip()

# Convert common "nan"/"None"/empty strings to NA before concatenation
genus_clean = genus_clean.replace({"nan": np.nan, "None": np.nan, "": np.nan})
species_clean = species_clean.replace({"nan": np.nan, "None": np.nan, "": np.nan})

fred["species"] = (genus_clean + " " + species_clean).str.replace(r"\s+", " ", regex=True).str.strip()

# Drop rows where we couldn't form a species name
fred = fred.loc[fred["species"].notna() & (fred["species"] != "")].copy()

unique_fred  = set(fred["species"].dropna().unique())
unique_roots = set(roots["plant_species"].dropna().astype(str).str.strip().unique())

overlap = sorted(unique_fred.intersection(unique_roots))

print(f"Unique FRED species:  {len(unique_fred):,}")
print(f"Unique roots species: {len(unique_roots):,}")
print(f"Overlap species:      {len(overlap):,}")


# %% [cell 9] type=code
# --- Block 8: Filter to overlap + coerce numeric columns in roots_vt ---

fred_roots = fred.loc[fred["species"].isin(overlap)].copy()
roots_vt   = roots.loc[roots["plant_species"].astype(str).str.strip().isin(overlap)].copy()

num_cols = [
    "prop.cult",
    "prop.ancestral",
    "prop.edaphophilic",
    "prop.rhizophilic",
    "spec.ances",
    "spec.edaph",
    "spec.rhizo",
    "richness",
]

for c in num_cols:
    if c in roots_vt.columns:
        roots_vt[c] = pd.to_numeric(roots_vt[c], errors="coerce")

print("fred_roots:", fred_roots.shape)
print("roots_vt:  ", roots_vt.shape)


# %% [cell 10] type=code
# --- Block 9: Species-level summaries + trait means merged ---

# VT/roots summaries
richness_vt = (
    roots_vt.groupby("plant_species", dropna=True)
    .agg(
        mean_richness=("richness", "mean"),
        mean_prop_cult=("prop.cult", "mean"),
        mean_prop_ancestral=("prop.ancestral", "mean"),
        mean_prop_edaphophilic=("prop.edaphophilic", "mean"),
        mean_prop_rhizophilic=("prop.rhizophilic", "mean"),
        mean_sp_ancestral=("spec.ances", "mean"),
        mean_sp_rhizophilic=("spec.rhizo", "mean"),
        mean_sp_edaphophilic=("spec.edaph", "mean"),
    )
    .reset_index()
    .rename(columns={"plant_species": "species"})
)

# Trait means from FRED
trait_map = {
    "Specific root length (SRL)": "mean_srl",
    "Root tissue density (RTD)": "mean_rtd",
    "Root diameter": "mean_rd",
    "Root N content": "mean_rnc",
}

trait_dfs = []
for trait_col, out_col in trait_map.items():
    if trait_col in fred_roots.columns:
        fred_roots[trait_col] = pd.to_numeric(fred_roots[trait_col], errors="coerce")
        tmp = (
            fred_roots.groupby("species", dropna=True)[trait_col]
            .mean()
            .rename(out_col)
            .reset_index()
        )
        trait_dfs.append(tmp)

# Merge: keep only species with VT summaries and trait data
df = richness_vt.copy()
for tdf in trait_dfs:
    df = df.merge(tdf, on="species", how="inner")

df = df.sort_values("mean_richness", ascending=False).reset_index(drop=True)

print("Merged df shape:", df.shape)
df.head()


# %% [cell 11] type=code
# --- Block 10: Scatter + regression helper + replicate plots ---

from scipy.stats import linregress

def scatter_lm(df, x, y, xlim=None, xlabel=None, ylabel=None, out_fp=None):
    d = df[[x, y]].dropna()
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(d[x], d[y], s=20)

    lr = linregress(d[x].values, d[y].values)
    xx = np.linspace(d[x].min(), d[x].max(), 100)
    ax.plot(xx, lr.intercept + lr.slope * xx)

    if xlim is not None:
        ax.set_xlim(xlim)

    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)

    ax.text(
        0.02, 0.98,
        f"r = {lr.rvalue:.2f}",
        transform=ax.transAxes,
        ha="left", va="top",
    )

    plt.tight_layout()
    if out_fp is not None:
        plt.savefig(out_fp, dpi=300)
    plt.show()

# Plots analogous to the R code
scatter_lm(df, "mean_richness", "mean_srl", xlim=(8, 37),
           xlabel="Mean VT richness", ylabel="Mean SRL",
           out_fp=DATA_DIR/"richness_vs_srl.jpg")

scatter_lm(df, "mean_richness", "mean_rtd", xlim=(8, 37),
           xlabel="Mean VT richness", ylabel="Mean RTD",
           out_fp=DATA_DIR/"richness_vs_rtd.jpg")

scatter_lm(df, "mean_richness", "mean_rd", xlim=(8, 37),
           xlabel="Mean VT richness", ylabel="Mean RD",
           out_fp=DATA_DIR/"richness_vs_rd.jpg")

scatter_lm(df, "mean_richness", "mean_rnc", xlim=(8, 37),
           xlabel="Mean VT richness", ylabel="Mean root N content",
           out_fp=DATA_DIR/"richness_vs_rnc.jpg")


# %% [cell 12] type=code
# --- Block: Report matching + non-matching (roots/ecobank vs FRED) genus+species ---

import re

# 1) Helper to split "Genus species" (keeps only first two tokens)
def split_binomial(name: str):
    if pd.isna(name):
        return (np.nan, np.nan)
    s = str(name).strip()
    if not s:
        return (np.nan, np.nan)
    parts = re.split(r"\s+", s)
    genus = parts[0] if len(parts) >= 1 else np.nan
    species = parts[1] if len(parts) >= 2 else np.nan
    return (genus, species)

# 2) Build ROOTS/ecobank binomials from plant_species
roots_species = (
    roots["plant_species"]
    .dropna()
    .astype(str)
    .str.replace(r"\s+", " ", regex=True)
    .str.strip()
    .drop_duplicates()
)
roots_binom = pd.DataFrame(
    roots_species.apply(split_binomial).tolist(),
    columns=["genus", "species"],
)
roots_binom["binomial"] = roots_binom["genus"].fillna("") + " " + roots_binom["species"].fillna("")
roots_binom["binomial"] = roots_binom["binomial"].str.strip()
roots_binom = roots_binom[(roots_binom["genus"].notna()) & (roots_binom["species"].notna())]

# 3) Build FRED binomials (assumes you already created fred["species"] as "Genus species")
fred_species = (
    fred["species"]
    .dropna()
    .astype(str)
    .str.replace(r"\s+", " ", regex=True)
    .str.strip()
    .drop_duplicates()
)
fred_binom = pd.DataFrame(
    fred_species.apply(split_binomial).tolist(),
    columns=["genus", "species"],
)
fred_binom["binomial"] = fred_binom["genus"].fillna("") + " " + fred_binom["species"].fillna("")
fred_binom["binomial"] = fred_binom["binomial"].str.strip()
fred_binom = fred_binom[(fred_binom["genus"].notna()) & (fred_binom["species"].notna())]

# 4) Compute sets
roots_set = set(roots_binom["binomial"].unique())
fred_set  = set(fred_binom["binomial"].unique())

match_set      = roots_set & fred_set
roots_only_set = roots_set - fred_set
fred_only_set  = fred_set - roots_set

# 5) Pretty print
def print_list(title, items, n=50):
    items = sorted(items)
    print("\n" + "=" * 80)
    print(f"{title}  (n={len(items):,})")
    print("=" * 80)
    for x in items[:n]:
        print(x)
    if len(items) > n:
        print(f"... ({len(items)-n:,} more)")

print_list("MATCHING binomials (roots/ecobank ∩ FRED)", match_set, n=50)
print_list("ROOTS/ECOBANK-ONLY binomials (roots/ecobank \\ FRED)", roots_only_set, n=50)
print_list("FRED-ONLY binomials (FRED \\ roots/ecobank)", fred_only_set, n=50)

# 6) Optional: genus-level overlap summary
roots_genera = set(roots_binom["genus"].unique())
fred_genera  = set(fred_binom["genus"].unique())

print("\n" + "=" * 80)
print("GENUS overlap summary")
print("=" * 80)
print(f"Roots/ecobank genera: {len(roots_genera):,}")
print(f"FRED genera:          {len(fred_genera):,}")
print(f"Overlap genera:       {len(roots_genera & fred_genera):,}")
print(f"Roots-only genera:    {len(roots_genera - fred_genera):,}")
print(f"FRED-only genera:     {len(fred_genera - roots_genera):,}")


# %% [cell 13] type=code
# --- Block: Find likely misspellings / near-matches between roots/ecobank and FRED (fixed) ---

import re
import difflib

import numpy as np
import pandas as pd

# ----------------------------
# 0) Build clean binomial sets from both datasets (self-contained)
# ----------------------------
def split_binomial(name: str):
    if pd.isna(name):
        return (np.nan, np.nan)
    s = str(name).strip()
    if not s:
        return (np.nan, np.nan)
    parts = re.split(r"\s+", s)
    genus = parts[0] if len(parts) >= 1 else np.nan
    species = parts[1] if len(parts) >= 2 else np.nan
    return (genus, species)

def clean_name_series(s: pd.Series) -> pd.Series:
    s = s.dropna().astype(str)
    s = s.str.replace(r"\s+", " ", regex=True).str.strip()
    # Normalize common artifacts
    s = s.str.replace(r"\u00d7", "x", regex=False)        # hybrid sign
    s = s.str.replace(r"[^\w\s\.-]", "", regex=True)      # drop punctuation except . and -
    s = s.str.replace(r"\s+", " ", regex=True).str.strip()
    return s.drop_duplicates()

roots_species = clean_name_series(roots["plant_species"])
fred_species  = clean_name_series(fred["species"])

roots_binom = pd.DataFrame(roots_species.apply(split_binomial).tolist(), columns=["genus", "species"])
fred_binom  = pd.DataFrame(fred_species.apply(split_binomial).tolist(),  columns=["genus", "species"])

# Keep only proper binomials (2 tokens)
roots_binom = roots_binom.dropna().copy()
fred_binom  = fred_binom.dropna().copy()

roots_binom["binomial"] = (roots_binom["genus"] + " " + roots_binom["species"]).str.strip()
fred_binom["binomial"]  = (fred_binom["genus"]  + " " + fred_binom["species"]).str.strip()

roots_set = set(roots_binom["binomial"].unique())
fred_set  = set(fred_binom["binomial"].unique())

roots_only = sorted(roots_set - fred_set)
fred_only  = sorted(fred_set - roots_set)

# ----------------------------
# 1) Similarity backend: RapidFuzz if available, else difflib
# ----------------------------
try:
    from rapidfuzz import fuzz, process
    HAVE_RAPIDFUZZ = True
except Exception:
    HAVE_RAPIDFUZZ = False

def best_match(query, choices):
    """Return (match, score_0_100)."""
    if HAVE_RAPIDFUZZ:
        m = process.extractOne(query, choices, scorer=fuzz.WRatio)
        return (m[0], float(m[1])) if m else (None, 0.0)
    else:
        matches = difflib.get_close_matches(query, choices, n=1, cutoff=0.0)
        if not matches:
            return (None, 0.0)
        score = difflib.SequenceMatcher(None, query, matches[0]).ratio() * 100.0
        return (matches[0], float(score))

# ----------------------------
# 2) Genus-constrained near-matches (fix: handle empty out cleanly)
# ----------------------------
def genus_constrained_pairs(roots_only, fred_only, min_score=90.0, max_pairs=200):
    roots_df = pd.DataFrame([split_binomial(x) for x in roots_only], columns=["genus", "species"])
    roots_df["binomial"] = roots_only
    fred_df = pd.DataFrame([split_binomial(x) for x in fred_only], columns=["genus", "species"])
    fred_df["binomial"] = fred_only

    out = []
    common_genera = sorted(set(roots_df["genus"]) & set(fred_df["genus"]))
    for g in common_genera:
        r_species = roots_df.loc[roots_df["genus"] == g, "species"].dropna().astype(str).unique().tolist()
        f_species = fred_df.loc[fred_df["genus"] == g, "species"].dropna().astype(str).unique().tolist()
        if not r_species or not f_species:
            continue

        if HAVE_RAPIDFUZZ:
            for rs in r_species:
                match = process.extractOne(rs, f_species, scorer=fuzz.WRatio)
                if match and float(match[1]) >= min_score:
                    fs, sc = match[0], float(match[1])
                    out.append({"genus": g, "roots": f"{g} {rs}", "fred": f"{g} {fs}", "score": sc, "method": "genus_constrained"})
        else:
            for rs in r_species:
                ms = difflib.get_close_matches(rs, f_species, n=1, cutoff=0.0)
                if ms:
                    fs = ms[0]
                    sc = difflib.SequenceMatcher(None, rs, fs).ratio() * 100.0
                    if sc >= min_score:
                        out.append({"genus": g, "roots": f"{g} {rs}", "fred": f"{g} {fs}", "score": float(sc), "method": "genus_constrained"})

        if len(out) >= max_pairs:
            break

    if not out:
        # Return an empty df with expected columns so downstream printing/saving is stable
        return pd.DataFrame(columns=["genus", "roots", "fred", "score", "method"])

    return (
        pd.DataFrame(out)
        .sort_values(["score", "genus"], ascending=[False, True])
        .head(max_pairs)
        .reset_index(drop=True)
    )

# ----------------------------
# 3) Global near-match search on full binomial (catches genus typos too)
# ----------------------------
def global_pairs(one_side, other_side, min_score=92.0, max_pairs=200, label_left="left", label_right="right"):
    out = []
    other_choices = list(other_side)
    for q in one_side:
        m, sc = best_match(q, other_choices)
        if m is not None and sc >= min_score:
            out.append({label_left: q, label_right: m, "score": sc, "method": "global_binomial"})
        if len(out) >= max_pairs:
            break

    if not out:
        return pd.DataFrame(columns=[label_left, label_right, "score", "method"])

    return pd.DataFrame(out).sort_values("score", ascending=False).head(max_pairs).reset_index(drop=True)

# ----------------------------
# 4) Run and print
# ----------------------------
MIN_SCORE_GENUS  = 80.0
MIN_SCORE_GLOBAL = 85.0

df_genus = genus_constrained_pairs(roots_only, fred_only, min_score=MIN_SCORE_GENUS, max_pairs=200)
df_r2f   = global_pairs(roots_only, fred_only, min_score=MIN_SCORE_GLOBAL, max_pairs=200, label_left="roots_only", label_right="fred_only_best")
df_f2r   = global_pairs(fred_only, roots_only, min_score=MIN_SCORE_GLOBAL, max_pairs=200, label_left="fred_only", label_right="roots_only_best")

print("\n" + "="*90)
print(f"LIKELY MISSPELLINGS / NEAR-MATCHES (GENUS-CONSTRAINED, score ≥ {MIN_SCORE_GENUS:g})")
print("="*90)
print("None found at this threshold." if df_genus.empty else df_genus.to_string(index=False))

print("\n" + "="*90)
print(f"LIKELY MISSPELLINGS / NEAR-MATCHES (GLOBAL BINOMIAL: roots_only → fred_only, score ≥ {MIN_SCORE_GLOBAL:g})")
print("="*90)
print("None found at this threshold." if df_r2f.empty else df_r2f.to_string(index=False))

print("\n" + "="*90)
print(f"LIKELY MISSPELLINGS / NEAR-MATCHES (GLOBAL BINOMIAL: fred_only → roots_only, score ≥ {MIN_SCORE_GLOBAL:g})")
print("="*90)
print("None found at this threshold." if df_f2r.empty else df_f2r.to_string(index=False))

# ----------------------------
# 5) Save (optional)
# ----------------------------
OUT_FP = Path("../Data/near_match_candidates.csv")

to_save = []

if not df_genus.empty:
    to_save.append(df_genus.rename(columns={"roots": "roots_name", "fred": "fred_name"}))

if not df_r2f.empty:
    to_save.append(df_r2f.rename(columns={"roots_only": "roots_name", "fred_only_best": "fred_name"}))

if not df_f2r.empty:
    to_save.append(df_f2r.rename(columns={"fred_only": "fred_name", "roots_only_best": "roots_name"}))

if to_save:
    pd.concat(to_save, ignore_index=True).drop_duplicates().to_csv(OUT_FP, index=False)
    print(f"\nSaved candidate near-matches to: {OUT_FP}")
else:
    print("\nNo candidates found; nothing saved.")


# %% [cell 14] type=code
# --- Block: Near-miss detection using canonicalized binomials (fixed for empty results) ---

import re
import difflib
import numpy as np
import pandas as pd
from pathlib import Path

# Optional faster/better fuzzy matching if installed
try:
    from rapidfuzz import fuzz, process
    HAVE_RAPIDFUZZ = True
except Exception:
    HAVE_RAPIDFUZZ = False

# ----------------------------
# 1) Canonicalize plant names to "Genus epithet"
# ----------------------------
QUALIFIERS = {
    "sp", "sp.", "spp", "spp.", "indet", "indet.", "cf", "cf.", "aff", "aff.",
    "nr", "nr.", "unidentified", "unknown"
}
INFRA = {"subsp", "subsp.", "ssp", "ssp.", "var", "var.", "forma", "forma.", "f", "f."}

def canonical_binomial(name: str):
    if pd.isna(name):
        return (np.nan, np.nan, np.nan)

    s = str(name).strip()
    if not s:
        return (np.nan, np.nan, np.nan)

    # Remove parenthetical authorship and normalize
    s = re.sub(r"\([^)]*\)", " ", s)
    s = s.replace("×", "x")
    s = re.sub(r"[^\w\s-]", " ", s)   # keep hyphens
    s = re.sub(r"\s+", " ", s).strip()

    parts = s.split(" ")
    if len(parts) < 2:
        return (np.nan, np.nan, np.nan)

    genus = parts[0]
    epithet = None

    for tok in parts[1:]:
        t = tok.lower().strip()

        if t in INFRA or t in QUALIFIERS:
            continue
        if len(t) == 1:
            continue
        if re.fullmatch(r"[a-z]+(-[a-z]+)?", t) is None:
            continue

        epithet = t
        break

    if epithet is None:
        return (np.nan, np.nan, np.nan)

    canon = f"{genus} {epithet}"
    return (canon, genus, epithet)

def build_name_table(series: pd.Series, source_label: str) -> pd.DataFrame:
    s = (
        series.dropna()
        .astype(str)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
        .drop_duplicates()
    )
    rows = [canonical_binomial(x) for x in s.tolist()]
    df = pd.DataFrame(rows, columns=["canon", "genus", "epithet"])
    df["raw"] = s.values
    df["source"] = source_label
    df = df.dropna(subset=["canon"]).drop_duplicates(subset=["canon", "raw"])
    return df

roots_tbl = build_name_table(roots["plant_species"], "roots/ecobank")
fred_tbl  = build_name_table(fred["species"], "FRED")

roots_canon = set(roots_tbl["canon"].unique())
fred_canon  = set(fred_tbl["canon"].unique())

roots_only_canon = sorted(roots_canon - fred_canon)
fred_only_canon  = sorted(fred_canon - roots_canon)

print(f"Canonical roots-only: {len(roots_only_canon):,}")
print(f"Canonical FRED-only:  {len(fred_only_canon):,}")

# ----------------------------
# 2) Fuzzy matching helpers
# ----------------------------
def best_match(query, choices):
    if HAVE_RAPIDFUZZ:
        m = process.extractOne(query, choices, scorer=fuzz.WRatio)
        return (m[0], float(m[1])) if m else (None, 0.0)
    else:
        matches = difflib.get_close_matches(query, choices, n=1, cutoff=0.0)
        if not matches:
            return (None, 0.0)
        sc = difflib.SequenceMatcher(None, query, matches[0]).ratio() * 100.0
        return (matches[0], float(sc))

def get_raw_examples(tbl: pd.DataFrame, canon: str, n=3):
    ex = tbl.loc[tbl["canon"] == canon, "raw"].head(n).tolist()
    return "; ".join(ex)

# ----------------------------
# 3) Genus-constrained near-matches among canonical names
# ----------------------------
MIN_SCORE_GENUS  = 80.0   # set low for exploration
MIN_SCORE_GLOBAL = 85.0

fred_by_genus = {
    g: sorted(fred_tbl.loc[fred_tbl["genus"] == g, "canon"].unique())
    for g in fred_tbl["genus"].dropna().unique()
}

genus_hits = []
for q in roots_only_canon:
    q_genus = q.split(" ", 1)[0]
    choices = fred_by_genus.get(q_genus, [])
    if not choices:
        continue
    m, sc = best_match(q, choices)
    if m and sc >= MIN_SCORE_GENUS:
        genus_hits.append({
            "roots_canon": q,
            "fred_canon": m,
            "score": sc,
            "roots_raw_examples": get_raw_examples(roots_tbl, q),
            "fred_raw_examples": get_raw_examples(fred_tbl, m),
            "method": "genus_constrained_canonical",
        })

# FIX: handle empty list before sorting
if genus_hits:
    df_genus = (
        pd.DataFrame(genus_hits)
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
else:
    df_genus = pd.DataFrame(columns=["roots_canon","fred_canon","score","roots_raw_examples","fred_raw_examples","method"])

# ----------------------------
# 4) Global near-matches among canonical names (roots_only_canon -> fred_only_canon)
# ----------------------------
global_hits = []
fred_choices = fred_only_canon[:]  # compare only among canonical non-overlaps

for q in roots_only_canon:
    m, sc = best_match(q, fred_choices)
    if m and sc >= MIN_SCORE_GLOBAL:
        global_hits.append({
            "roots_canon": q,
            "fred_canon": m,
            "score": sc,
            "roots_raw_examples": get_raw_examples(roots_tbl, q),
            "fred_raw_examples": get_raw_examples(fred_tbl, m),
            "method": "global_canonical",
        })

if global_hits:
    df_global = (
        pd.DataFrame(global_hits)
        .sort_values("score", ascending=False)
        .reset_index(drop=True)
    )
else:
    df_global = pd.DataFrame(columns=["roots_canon","fred_canon","score","roots_raw_examples","fred_raw_examples","method"])

# ----------------------------
# 5) Print results
# ----------------------------
print("\n" + "="*110)
print(f"GENUS-CONSTRAINED NEAR-MATCHES (canonical), score ≥ {MIN_SCORE_GENUS:g}")
print("="*110)
print("None found." if df_genus.empty else df_genus.head(50).to_string(index=False))

print("\n" + "="*110)
print(f"GLOBAL NEAR-MATCHES (canonical), score ≥ {MIN_SCORE_GLOBAL:g}")
print("="*110)
print("None found." if df_global.empty else df_global.head(50).to_string(index=False))

# ----------------------------
# 6) Optional: save
# ----------------------------
OUT_FP = Path("../Data/near_match_candidates_canonical.csv")
df_out = pd.concat([df_genus, df_global], ignore_index=True).drop_duplicates()

if not df_out.empty:
    df_out.to_csv(OUT_FP, index=False)
    print(f"\nSaved: {OUT_FP}")
else:
    print("\nNo near-match candidates found (even after canonicalization).")


# %% [cell 15] type=code
# how many 'sp.' / qualifiers are in roots-only vs fred-only?
pd.Series([x.split(" ",1)[1] if " " in x else "" for x in roots_only]).value_counts().head(20)
pd.Series([x.split(" ",1)[1] if " " in x else "" for x in fred_only]).value_counts().head(20)


# %% [cell 16] type=code
# --- Block: Keep only "likely misspellings" (not just nearest neighbor) + safe concat (no FutureWarning) ---
import os
import re
import difflib
import numpy as np
import pandas as pd
from pathlib import Path

def norm_for_compare(s: str) -> str:
    s = str(s).lower().strip()
    s = re.sub(r"\s+", " ", s)
    return s

def epithet(s: str) -> str:
    # expects "Genus epithet"
    parts = str(s).split(" ", 1)
    return parts[1].strip().lower() if len(parts) == 2 else ""

def is_likely_misspelling(a: str, b: str, max_edit=2, min_prefix=5) -> bool:
    """
    Heuristic:
      - same genus required (caller should ensure)
      - epithet differs by small edit distance OR shares long prefix
    """
    ea, eb = epithet(a), epithet(b)
    if not ea or not eb or ea == eb:
        return False

    # edit distance via SequenceMatcher "opcodes" (approx)
    # Convert similarity to an estimated edit distance bound:
    sm = difflib.SequenceMatcher(None, ea, eb)
    ratio = sm.ratio()
    # crude upper bound: edits <= max_edit tends to give high ratio; we also allow long shared prefix
    common_prefix = len(os.path.commonprefix([ea, eb])) if "os" in globals() else len(re.match(r"^([a-z]+)", ea).group(1)) if False else len(os.path.commonprefix([ea, eb]))

    # quick exact Levenshtein if rapidfuzz is available
    try:
        from rapidfuzz.distance import Levenshtein
        ed = Levenshtein.distance(ea, eb)
    except Exception:
        # fallback: approximate edit distance from ratio (rough)
        ed = int(round((1 - ratio) * max(len(ea), len(eb))))

    return (ed <= max_edit) or (common_prefix >= min_prefix)

# 1) Load previously saved candidate file (from your earlier block)
IN_FP = Path("../Data/near_match_candidates_canonical.csv")
assert IN_FP.exists(), IN_FP
cand = pd.read_csv(IN_FP)

# 2) Keep only genus-constrained rows (those are interpretable)
cand = cand.loc[cand["method"].eq("genus_constrained_canonical")].copy()

# 3) Require same genus explicitly + apply misspelling heuristic
cand["roots_genus"] = cand["roots_canon"].astype(str).str.split(" ", n=1).str[0]
cand["fred_genus"]  = cand["fred_canon"].astype(str).str.split(" ", n=1).str[0]
cand = cand.loc[cand["roots_genus"].eq(cand["fred_genus"])].copy()

cand["likely_misspelling"] = [
    is_likely_misspelling(r, f, max_edit=2, min_prefix=6)
    for r, f in zip(cand["roots_canon"], cand["fred_canon"])
]

cand_miss = cand.loc[cand["likely_misspelling"]].copy()

print("\n" + "="*110)
print("LIKELY MISSPELLINGS (same genus; small edit distance / long shared prefix)")
print("="*110)
if cand_miss.empty:
    print("None detected under this heuristic. (Differences are likely true taxonomic non-overlap.)")
else:
    print(cand_miss[["roots_canon","fred_canon","score","roots_raw_examples","fred_raw_examples"]].to_string(index=False))

# 4) Save filtered list
OUT_FP = Path("../Data/likely_misspellings_only.csv")
cand_miss.to_csv(OUT_FP, index=False)
print(f"\nSaved: {OUT_FP}")
