from __future__ import annotations

from dataclasses import dataclass
from difflib import SequenceMatcher, get_close_matches
from typing import Iterable

from .text import genus_of


@dataclass
class MatchStats:
    dataset: str
    fred_source: str
    n_dataset_species: int
    n_fred_species: int
    exact_matches: int
    genus_overlap: int


def compute_match_stats(dataset: str, fred_source: str, dataset_species: set[str], fred_species: set[str]) -> MatchStats:
    exact = dataset_species & fred_species
    dataset_genera = {genus_of(s) for s in dataset_species if genus_of(s)}
    fred_genera = {genus_of(s) for s in fred_species if genus_of(s)}
    genus_overlap = len(dataset_genera & fred_genera)
    return MatchStats(
        dataset=dataset,
        fred_source=fred_source,
        n_dataset_species=len(dataset_species),
        n_fred_species=len(fred_species),
        exact_matches=len(exact),
        genus_overlap=genus_overlap,
    )


def likely_matches(unmatched: Iterable[str], fred_species: set[str], cutoff: float = 0.86, top_n: int = 50) -> list[dict[str, object]]:
    fred_list = sorted(fred_species)
    out: list[dict[str, object]] = []
    for name in sorted(set(unmatched)):
        g = genus_of(name)
        in_genus = [x for x in fred_list if genus_of(x) == g] if g else []
        pool = in_genus if in_genus else fred_list
        candidates = get_close_matches(name, pool, n=1, cutoff=cutoff)
        if not candidates:
            continue
        cand = candidates[0]
        ratio = SequenceMatcher(None, name, cand).ratio()
        out.append(
            {
                "unmatched_name": name,
                "candidate_fred_name": cand,
                "similarity_ratio": round(ratio, 3),
                "same_genus": int(genus_of(name) == genus_of(cand)),
            }
        )
    out.sort(key=lambda d: (d["similarity_ratio"], d["same_genus"]), reverse=True)
    return out[:top_n]
