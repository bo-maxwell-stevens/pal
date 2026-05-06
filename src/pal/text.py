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
    if not species_name:
        return None
    return species_name.split(" ", 1)[0]
