from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd


def read_fred_table(path: Path) -> pd.DataFrame:
    """Read FRED CSV where row 2 holds headers and data starts row 3."""
    last_err = None
    for enc in ("utf-8", "latin-1", "cp1252"):
        try:
            raw = pd.read_csv(path, header=None, dtype=str, encoding=enc, low_memory=False)
            break
        except UnicodeDecodeError as exc:
            last_err = exc
    else:
        raise last_err  # type: ignore[misc]

    if raw.shape[0] < 3:
        raise ValueError(f"FRED file has too few rows: {path}")

    raw.columns = raw.iloc[1].astype(str)
    df = raw.iloc[2:].copy()
    df.reset_index(drop=True, inplace=True)

    row_id_col = "Notes_Row ID"
    if row_id_col in df.columns:
        row_id = df[row_id_col].astype(str).str.strip()
        numeric_mask = row_id.str.fullmatch(r"\d+")
        if numeric_mask.any():
            df = df[numeric_mask].copy()
            df.reset_index(drop=True, inplace=True)
    return df


def first_present(columns: Sequence[str], candidates: Iterable[str]) -> str | None:
    for name in candidates:
        if name in columns:
            return name
    return None
