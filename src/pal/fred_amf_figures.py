from __future__ import annotations

import string
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm


def _panel_labels(axs):
    for i, ax in enumerate(axs.ravel()):
        ax.text(0.01, 0.99, string.ascii_uppercase[i], transform=ax.transAxes, va="top", ha="left", fontweight="bold")


def save_png_pdf(fig, base: Path) -> None:
    fig.tight_layout()
    fig.savefig(base.with_suffix(".png"), dpi=300)
    fig.savefig(base.with_suffix(".pdf"))
    plt.close(fig)


def fig1_overlap(summary: pd.DataFrame, out_base: Path) -> None:
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.5))
    d = summary.copy()
    order = ["FRED3", "FRED4_filtered_fineroot_lt2mm", "FRED4_filtered_1storder", "FRED4_full"]
    for ax, y, title in zip(
        axs,
        ["exact_matches", "match_percent", "genus_overlap"],
        ["Exact species matches", "Match percent", "Genus overlap"],
    ):
        x = np.arange(len(order))
        w = 0.35
        eco = [d[(d.dataset == "EcoBank") & (d.fred_source == o)][y].mean() for o in order]
        glo = [d[(d.dataset == "GlobalAMFungi") & (d.fred_source == o)][y].mean() for o in order]
        ax.bar(x - w / 2, eco, width=w, label="EcoBank")
        ax.bar(x + w / 2, glo, width=w, label="GlobalAMFungi")
        ax.set_xticks(x)
        ax.set_xticklabels(order, rotation=20, ha="right")
        ax.set_title(title)
    axs[0].legend(frameon=False)
    _panel_labels(axs)
    save_png_pdf(fig, out_base)


def fig2_geo(fred4: pd.DataFrame, eco_sample: pd.DataFrame, glo_sample: pd.DataFrame, out_base: Path) -> None:
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    panels = [
        (fred4, "Longitude_Main", "Latitude_Main", "FRED4 full"),
        (eco_sample, "lon", "lat", "EcoBank"),
        (glo_sample, "longitude", "latitude", "GlobalAMFungi"),
    ]
    for ax, (df, x, y, title) in zip(axs.ravel()[:3], panels):
        if x in df.columns and y in df.columns:
            xx = pd.to_numeric(df[x], errors="coerce")
            yy = pd.to_numeric(df[y], errors="coerce")
            ax.scatter(xx, yy, s=2, alpha=0.3)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-60, 85)
        ax.set_title(title)
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

    ax = axs.ravel()[3]
    if all(c in fred4.columns for c in ["canonical_species", "Longitude_Main", "Latitude_Main"]):
        overlap = set(eco_sample.get("canonical_species", [])) | set(glo_sample.get("canonical_species", []))
        d = fred4[fred4["canonical_species"].isin(overlap)]
        ax.scatter(pd.to_numeric(d["Longitude_Main"], errors="coerce"), pd.to_numeric(d["Latitude_Main"], errors="coerce"), s=2, alpha=0.3)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 85)
    ax.set_title("Overlapping taxa locations")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    _panel_labels(axs)
    save_png_pdf(fig, out_base)


def fig3_trait_distributions(master: pd.DataFrame, out_base: Path) -> None:
    traits = ["root_diameter_mean", "SRL_mean", "RTD_mean", "root_N_mean", "root_P_mean", "mycorrhizal_colonization_mean"]
    traits = [t for t in traits if t in master.columns]
    n = len(traits)
    fig, axs = plt.subplots(2, int(np.ceil(n / 2)), figsize=(14, 7))
    axs = np.array(axs).reshape(2, -1)

    d_all = master[master["in_FRED4_full"] == 1]
    d_eco = master[(master["in_FRED4_full"] == 1) & (master["in_EcoBank"] == 1)]
    d_glo = master[(master["in_FRED4_full"] == 1) & (master["in_GlobalAMFungi"] == 1)]

    for i, t in enumerate(traits):
        ax = axs.ravel()[i]
        for data, label, color in [(d_all, "All FRED4", "#4c78a8"), (d_eco, "Overlap EcoBank", "#f58518"), (d_glo, "Overlap GlobalAMFungi", "#54a24b")]:
            v = pd.to_numeric(data[t], errors="coerce").dropna()
            if len(v):
                ax.hist(v, bins=25, alpha=0.35, density=True, label=label, color=color)
        ax.set_title(t)
    for j in range(i + 1, axs.size):
        axs.ravel()[j].axis("off")
    axs.ravel()[0].legend(frameon=False)
    _panel_labels(axs)
    save_png_pdf(fig, out_base)


def fig4_pca(master: pd.DataFrame, out_base: Path) -> None:
    d = master.dropna(subset=["PC1", "PC2"]).copy()
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))

    color_wood = d.get("woodiness", pd.Series(index=d.index, dtype=str)).astype(str).str.lower().str.contains("wood")
    axs[0, 0].scatter(d["PC1"], d["PC2"], c=np.where(color_wood, 1, 0), s=18, alpha=0.7)
    axs[0, 0].set_title("PCA by woodiness")

    mem = d["in_EcoBank"].astype(int) + 2 * d["in_GlobalAMFungi"].astype(int)
    axs[0, 1].scatter(d["PC1"], d["PC2"], c=mem, s=18, alpha=0.7)
    axs[0, 1].set_title("PCA by AMF dataset membership")

    for ax, x in [(axs[1, 0], "PC1"), (axs[1, 1], "PC2")]:
        for ycol, label in [("EcoBank_amf_richness_mean", "EcoBank"), ("GlobalAMFungi_amf_richness_mean", "GlobalAMFungi")]:
            if ycol in d.columns:
                dd = d[[x, ycol]].dropna()
                if len(dd):
                    ax.scatter(dd[x], dd[ycol], s=16, alpha=0.5, label=label)
        ax.set_title(f"{x} vs AMF richness")
        ax.set_xlabel(x)
        ax.set_ylabel("AMF richness")
    axs[1, 0].legend(frameon=False)
    _panel_labels(axs)
    save_png_pdf(fig, out_base)


def _scatter_with_fit(ax, d: pd.DataFrame, x: str, y: str, label: str) -> tuple[int, float, float]:
    dd = d[[x, y]].dropna()
    if len(dd) < 5:
        return len(dd), np.nan, np.nan
    ax.scatter(dd[x], dd[y], s=18, alpha=0.55, label=label)
    X = sm.add_constant(dd[[x]])
    model = sm.OLS(dd[y], X).fit()
    xx = np.linspace(dd[x].min(), dd[x].max(), 100)
    yy = model.params["const"] + model.params[x] * xx
    ax.plot(xx, yy, linewidth=1.5)
    return len(dd), model.rsquared, float(model.pvalues[x])


def fig5_richness_vs_traits(master: pd.DataFrame, out_base: Path) -> None:
    traits = ["root_diameter_mean", "SRL_mean", "RTD_mean", "root_N_mean", "root_P_mean"]
    fig, axs = plt.subplots(2, 3, figsize=(14, 8))
    for i, t in enumerate(traits):
        ax = axs.ravel()[i]
        text = []
        for y, label in [("EcoBank_amf_richness_mean", "EcoBank"), ("GlobalAMFungi_amf_richness_mean", "GlobalAMFungi")]:
            if y in master.columns:
                n, r2, p = _scatter_with_fit(ax, master, t, y, label)
                if n:
                    text.append(f"{label}: n={n}, R2={r2:.2f}, p={p:.3g}")
        ax.set_title(t)
        ax.set_xlabel(t)
        ax.set_ylabel("AMF richness")
        if text:
            ax.text(0.03, 0.97, "\n".join(text), transform=ax.transAxes, va="top", fontsize=8)
    axs.ravel()[0].legend(frameon=False)
    axs.ravel()[-1].axis("off")
    _panel_labels(axs)
    save_png_pdf(fig, out_base)


def fig6_guilds_vs_traits(master: pd.DataFrame, out_base: Path) -> None:
    guilds = ["EcoBank_prop_cult_mean", "EcoBank_prop_ancestral_mean", "EcoBank_prop_rhizophilic_mean", "EcoBank_prop_edaphophilic_mean"]
    traits = ["root_diameter_mean", "SRL_mean", "RTD_mean", "root_N_mean", "root_P_mean"]
    fig, axs = plt.subplots(len(guilds), len(traits), figsize=(16, 12), sharey=False)
    for i, g in enumerate(guilds):
        if g not in master.columns:
            continue
        for j, t in enumerate(traits):
            ax = axs[i, j]
            dd = master[[g, t]].dropna().copy()
            if len(dd) < 6:
                ax.set_axis_off()
                continue
            n = len(dd)
            p_adj = (dd[g] * (n - 1) + 0.5) / n
            dd["logit"] = np.log(p_adj.clip(1e-6, 1 - 1e-6) / (1 - p_adj.clip(1e-6, 1 - 1e-6)))
            nn, r2, p = _scatter_with_fit(ax, dd.rename(columns={"logit": g + "_logit"}), t, g + "_logit", "")
            ax.set_title(f"{g} vs {t}", fontsize=9)
            ax.text(0.03, 0.97, f"n={nn}, R2={r2:.2f}, p={p:.3g}", transform=ax.transAxes, va="top", fontsize=7)
    _panel_labels(axs)
    save_png_pdf(fig, out_base)


def fig7_woody_sensitivity(master: pd.DataFrame, out_base: Path) -> None:
    d = master.copy()
    d["group"] = np.where(d.get("woodiness", "").astype(str).str.lower().str.contains("wood"), "woody", "non_woody")
    traits = ["root_diameter_mean", "RTD_mean", "root_N_mean"]
    fig, axs = plt.subplots(2, len(traits), figsize=(14, 7))
    for r, resp in enumerate(["EcoBank_amf_richness_mean", "EcoBank_prop_rhizophilic_mean"]):
        for c, t in enumerate(traits):
            ax = axs[r, c]
            for g, dd in d.groupby("group"):
                col = resp
                data = dd[[t, col]].dropna().copy()
                if resp.endswith("prop_rhizophilic_mean") and len(data):
                    n = len(data)
                    p_adj = (data[col] * (n - 1) + 0.5) / n
                    data[col] = np.log(p_adj.clip(1e-6, 1 - 1e-6) / (1 - p_adj.clip(1e-6, 1 - 1e-6)))
                _scatter_with_fit(ax, data, t, col, g)
            ax.set_title(f"{resp} vs {t}")
    axs[0, 0].legend(frameon=False)
    _panel_labels(axs)
    save_png_pdf(fig, out_base)


def fig8_genus_sensitivity(genus_master: pd.DataFrame, out_base: Path) -> None:
    traits = ["root_diameter_mean", "SRL_mean", "RTD_mean", "root_N_mean"]
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    for ax, t in zip(axs.ravel(), traits):
        if t in genus_master.columns and "EcoBank_amf_richness_mean" in genus_master.columns:
            d = genus_master[[t, "EcoBank_amf_richness_mean", "n_species"]].dropna()
            if len(d):
                ax.scatter(d[t], d["EcoBank_amf_richness_mean"], s=10 + d["n_species"] * 2, alpha=0.5)
                _scatter_with_fit(ax, d, t, "EcoBank_amf_richness_mean", "")
            ax.set_title(f"Genus richness vs {t}")
    _panel_labels(axs)
    save_png_pdf(fig, out_base)
