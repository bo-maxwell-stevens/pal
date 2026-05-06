from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .loaders import read_fred_table
from .text import canonical_binomial, canonicalize_species, genus_of


MANUAL_CORRECTIONS = {
    "deschampsia caespitosa": "deschampsia cespitosa",
    "hypochoeris radicata": "hypochaeris radicata",
}


FRED_SOURCES = {
    "FRED3": "Data/FRED3_fineRoots.csv",
    "FRED4_filtered_fineroot_lt2mm": "Data/FRED_4_20250921_filteredforMicrobeNet_AMF_FineRootsLess2mm.csv",
    "FRED4_filtered_1storder": "Data/FRED_4_20250921_filteredforMicrobeNet_AMF_1stOrderRoots.csv",
    "FRED4_full": "Data/FRED_4_full_20260312_2.csv",
}


def ensure_dirs(root: Path) -> dict[str, Path]:
    out = root / "Output"
    paths = {
        "output": out,
        "tables": out / "Tables",
        "figures": out / "Figures",
        "reports": out / "Reports",
        "intermediate": out / "Intermediate",
        "scripts": root / "scripts",
        "src": root / "src",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def _apply_corrections(series: pd.Series) -> pd.Series:
    s = series.map(canonical_binomial)
    return s.map(lambda x: MANUAL_CORRECTIONS.get(x, x) if x else x)


def load_ecobank(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    eco = pd.read_csv(root / "Data/ecobank_full.csv", sep=";", low_memory=False)
    cult = pd.read_csv(root / "Data/ecobank_cultured.csv", sep=";", low_memory=False)

    eco["canonical_species"] = _apply_corrections(eco.get("plant_species"))
    eco = eco[eco["canonical_species"].notna()].copy()

    merged = eco.merge(cult, on="sample", how="left")
    for c in ["richness", "prop.cult", "prop.ancestral", "prop.edaphophilic", "prop.rhizophilic"]:
        if c in merged.columns:
            merged[c] = pd.to_numeric(merged[c], errors="coerce")

    # sample-level metrics projected to species
    species = (
        merged.groupby("canonical_species", as_index=False)
        .agg(
            EcoBank_n_samples=("sample", "nunique"),
            EcoBank_amf_richness_mean=("richness", "mean"),
            EcoBank_prop_cult_mean=("prop.cult", "mean"),
            EcoBank_prop_ancestral_mean=("prop.ancestral", "mean"),
            EcoBank_prop_edaphophilic_mean=("prop.edaphophilic", "mean"),
            EcoBank_prop_rhizophilic_mean=("prop.rhizophilic", "mean"),
            EcoBank_latitude_mean=("lat", "mean"),
            EcoBank_longitude_mean=("lon", "mean"),
        )
    )
    return merged, species


def load_globalamf(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    fp = root / "Output/globalamfungi_sample_level_full.csv"
    if not fp.exists():
        fp = root / "Output/globalamfungi_sample_level.csv"
    g = pd.read_csv(fp)
    g["canonical_species"] = _apply_corrections(g.get("plant_species"))
    g = g[g["canonical_species"].notna()].copy()
    for c in ["amf_seq_richness", "amf_genus_richness", "latitude", "longitude"]:
        if c in g.columns:
            g[c] = pd.to_numeric(g[c], errors="coerce")
    species = (
        g.groupby("canonical_species", as_index=False)
        .agg(
            GlobalAMFungi_n_samples=("id", "nunique"),
            GlobalAMFungi_amf_richness_mean=("amf_seq_richness", "mean"),
            GlobalAMFungi_amf_genus_richness_mean=("amf_genus_richness", "mean"),
            GlobalAMFungi_latitude_mean=("latitude", "mean"),
            GlobalAMFungi_longitude_mean=("longitude", "mean"),
        )
    )
    return g, species


def load_fred_species(root: Path, path: Path) -> pd.DataFrame:
    df = read_fred_table(path)
    genus = df.get("Plant taxonomy_Accepted genus_WFO")
    if genus is None:
        genus = df.get("Plant taxonomy_Accepted genus_TPL")
    species = df.get("Plant Taxonomy_Accepted species_WFO")
    if species is None:
        species = df.get("Plant Taxonomy_Accepted species_TPL")
    if species is None:
        species = df.get("Plant taxonomy_Species_Data source")

    if genus is not None and species is not None:
        candidate = genus.fillna("").astype(str).str.strip() + " " + species.fillna("").astype(str).str.strip()
        candidate = candidate.str.strip()
    else:
        candidate = df.get("Name", pd.Series(index=df.index, dtype=str)).astype(str)

    df = df.copy()
    df["canonical_species"] = _apply_corrections(candidate)
    df["genus"] = df["canonical_species"].map(genus_of)
    return df


def overlap_summary(eco_species: set[str], global_species: set[str], fred_sources: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for name, fdf in fred_sources.items():
        fset = set(x for x in fdf["canonical_species"].dropna().unique().tolist() if x)
        fgen = set(x for x in fdf["genus"].dropna().unique().tolist() if x)
        for dataset, sset in [("EcoBank", eco_species), ("GlobalAMFungi", global_species)]:
            sgen = {genus_of(s) for s in sset if genus_of(s)}
            exact = len(sset & fset)
            rows.append(
                {
                    "dataset": dataset,
                    "fred_source": name,
                    "n_dataset_species": len(sset),
                    "n_fred_species": len(fset),
                    "exact_matches": exact,
                    "match_percent": 100 * exact / max(len(sset), 1),
                    "genus_overlap": len(sgen & fgen),
                }
            )
    return pd.DataFrame(rows)


def build_fred4_species_master(fred4: pd.DataFrame, eco_sp: pd.DataFrame, glo_sp: pd.DataFrame) -> pd.DataFrame:
    df = fred4.copy()
    numeric_cols = [
        "Root diameter",
        "Specific root length (SRL)",
        "Root tissue density (RTD)",
        "Root N content",
        "Root P content",
        "Mycorrhiza_Fraction of root length or tips colonized",
    ]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    grp_base = (
        df.groupby("canonical_species", as_index=False)
        .agg(
            woodiness=("Plant woodiness_TRY", lambda s: s.dropna().astype(str).mode().iloc[0] if s.dropna().size else np.nan),
            growth_form=("Plant growth form_TRY", lambda s: s.dropna().astype(str).mode().iloc[0] if s.dropna().size else np.nan),
            n_FRED_records=("Notes_Row ID", "count"),
        )
    )

    mean_frames = []
    count_frames = []
    for c in numeric_cols:
        if c in df.columns:
            mean_frames.append(df.groupby("canonical_species", as_index=False)[c].mean().rename(columns={c: c + "__mean"}))
            count_frames.append(
                df.assign(**{c: pd.to_numeric(df[c], errors="coerce")})
                .groupby("canonical_species", as_index=False)[c]
                .count()
                .rename(columns={c: c + "__count"})
            )

    grp = grp_base
    for fr in mean_frames + count_frames:
        grp = grp.merge(fr, on="canonical_species", how="left")

    grp = grp.rename(
        columns={
            "Root diameter__mean": "root_diameter_mean",
            "Specific root length (SRL)__mean": "SRL_mean",
            "Root tissue density (RTD)__mean": "RTD_mean",
            "Root N content__mean": "root_N_mean",
            "Root P content__mean": "root_P_mean",
            "Mycorrhiza_Fraction of root length or tips colonized__mean": "mycorrhizal_colonization_mean",
            "Root diameter__count": "n_root_diameter_records",
            "Specific root length (SRL)__count": "n_SRL_records",
            "Root tissue density (RTD)__count": "n_RTD_records",
            "Root N content__count": "n_root_N_records",
            "Root P content__count": "n_root_P_records",
            "Mycorrhiza_Fraction of root length or tips colonized__count": "n_mycorrhizal_colonization_records",
        }
    )

    master = grp.merge(eco_sp, on="canonical_species", how="outer").merge(glo_sp, on="canonical_species", how="outer")
    master["in_EcoBank"] = master["EcoBank_n_samples"].fillna(0).gt(0).astype(int)
    master["in_GlobalAMFungi"] = master["GlobalAMFungi_n_samples"].fillna(0).gt(0).astype(int)
    master["in_FRED4_full"] = master["n_FRED_records"].fillna(0).gt(0).astype(int)
    return master


def build_genus_master(species_master: pd.DataFrame) -> pd.DataFrame:
    s = species_master.copy()
    s["genus"] = s["canonical_species"].map(genus_of)
    numeric = [
        "root_diameter_mean",
        "SRL_mean",
        "RTD_mean",
        "root_N_mean",
        "root_P_mean",
        "mycorrhizal_colonization_mean",
        "EcoBank_amf_richness_mean",
        "GlobalAMFungi_amf_richness_mean",
        "EcoBank_prop_cult_mean",
        "EcoBank_prop_ancestral_mean",
        "EcoBank_prop_edaphophilic_mean",
        "EcoBank_prop_rhizophilic_mean",
    ]
    keep = [c for c in numeric if c in s.columns]
    agg = {c: "mean" for c in keep}
    agg.update(
        {
            "canonical_species": "count",
            "n_FRED_records": "sum",
            "in_EcoBank": "sum",
            "in_GlobalAMFungi": "sum",
        }
    )
    g = s[s["genus"].notna()].groupby("genus", as_index=False).agg(agg)
    g = g.rename(
        columns={
            "canonical_species": "n_species",
            "n_FRED_records": "n_fred_records",
            "in_EcoBank": "n_species_in_ecobank",
            "in_GlobalAMFungi": "n_species_in_globalamfungi",
        }
    )
    return g


def _fit_ols_rows(df: pd.DataFrame, response: str, predictors: list[str], dataset: str, level: str) -> list[dict[str, object]]:
    cols = [response] + predictors
    d = df[cols].dropna().copy()
    if len(d) < (len(predictors) + 5):
        return []
    X = sm.add_constant(d[predictors], has_constant="add")
    y = d[response]
    model = sm.OLS(y, X).fit()
    out = []
    for p in ["const"] + predictors:
        out.append(
            {
                "analysis_level": level,
                "dataset": dataset,
                "response": response,
                "predictor": p,
                "coef": model.params.get(p, np.nan),
                "se": model.bse.get(p, np.nan),
                "p_value": model.pvalues.get(p, np.nan),
                "r2": model.rsquared,
                "n": int(model.nobs),
            }
        )
    return out


def logit_transform(p: pd.Series, n: int) -> pd.Series:
    p_adj = (p * (n - 1) + 0.5) / n
    p_adj = p_adj.clip(1e-6, 1 - 1e-6)
    return np.log(p_adj / (1 - p_adj))


def run_models_species(master: pd.DataFrame) -> pd.DataFrame:
    predictors = ["root_diameter_mean", "SRL_mean", "RTD_mean", "root_N_mean"]
    if "woodiness" in master.columns:
        master = master.copy()
        master["woodiness_bin"] = master["woodiness"].astype(str).str.lower().str.contains("wood").astype(float)
        predictors = predictors + ["woodiness_bin"]

    rows: list[dict[str, object]] = []
    for dataset, resp in [
        ("EcoBank", "EcoBank_amf_richness_mean"),
        ("GlobalAMFungi", "GlobalAMFungi_amf_richness_mean"),
    ]:
        if resp in master.columns:
            rows.extend(_fit_ols_rows(master, resp, predictors, dataset, "species"))

    for resp in ["EcoBank_prop_cult_mean", "EcoBank_prop_ancestral_mean", "EcoBank_prop_edaphophilic_mean", "EcoBank_prop_rhizophilic_mean"]:
        if resp in master.columns:
            d = master.copy()
            d[resp + "_logit"] = logit_transform(pd.to_numeric(d[resp], errors="coerce"), max(len(d), 2))
            rows.extend(_fit_ols_rows(d, resp + "_logit", predictors, "EcoBank", "species"))

    return pd.DataFrame(rows)


def run_models_genus(genus_master: pd.DataFrame) -> pd.DataFrame:
    predictors = [c for c in ["root_diameter_mean", "SRL_mean", "RTD_mean", "root_N_mean"] if c in genus_master.columns]
    rows: list[dict[str, object]] = []
    for dataset, resp in [
        ("EcoBank", "EcoBank_amf_richness_mean"),
        ("GlobalAMFungi", "GlobalAMFungi_amf_richness_mean"),
    ]:
        if resp in genus_master.columns:
            rows.extend(_fit_ols_rows(genus_master, resp, predictors, dataset, "genus"))

    for resp in ["EcoBank_prop_cult_mean", "EcoBank_prop_ancestral_mean", "EcoBank_prop_edaphophilic_mean", "EcoBank_prop_rhizophilic_mean"]:
        if resp in genus_master.columns:
            d = genus_master.copy()
            d[resp + "_logit"] = logit_transform(pd.to_numeric(d[resp], errors="coerce"), max(len(d), 2))
            rows.extend(_fit_ols_rows(d, resp + "_logit", predictors, "EcoBank", "genus"))

    return pd.DataFrame(rows)


def add_pca(master: pd.DataFrame) -> pd.DataFrame:
    traits = ["root_diameter_mean", "SRL_mean", "RTD_mean", "root_N_mean"]
    d = master.copy()
    comp = d[traits].dropna()
    if len(comp) < 5:
        d["PC1"] = np.nan
        d["PC2"] = np.nan
        return d
    scaler = StandardScaler()
    X = scaler.fit_transform(comp)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X)
    d["PC1"] = np.nan
    d["PC2"] = np.nan
    d.loc[comp.index, "PC1"] = pcs[:, 0]
    d.loc[comp.index, "PC2"] = pcs[:, 1]
    return d
