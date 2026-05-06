from __future__ import annotations

import re
from typing import Optional


def canonicalize_species(value: object) -> Optional[str]:
    if value is None:
        return None
    s = str(value).strip().lower()
    if not s or s in {"nan", "none", "na", "n/a"}:
        return None
    s = re.sub(r"\s+", " ", s)
    s = re.sub(r"[\.;:,]+$", "", s)
    return s


def genus_of(species_name: Optional[str]) -> Optional[str]:
    if species_name is None:
        return None
    s = str(species_name).strip()
    if not s or s.lower() in {"nan", "none"}:
        return None
    return s.split(" ", 1)[0]


def canonical_binomial(value: object) -> Optional[str]:
    s = canonicalize_species(value)
    if not s:
        return None
    tokens = [t for t in s.split(" ") if t]
    if len(tokens) < 2:
        return None
    genus, epithet = tokens[0], tokens[1]
    bad = {"sp", "sp.", "spp", "cf", "cf.", "aff", "aff.", "x", "hybrid"}
    if epithet in bad:
        return None
    if not re.match(r"^[a-z-]+$", genus):
        return None
    if not re.match(r"^[a-z-]+$", epithet):
        return None
    return f"{genus} {epithet}"
