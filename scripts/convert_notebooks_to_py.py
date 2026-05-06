#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path


def notebook_to_py(nb_path: Path, out_path: Path) -> None:
    data = json.loads(nb_path.read_text(encoding="utf-8"))
    lines: list[str] = []
    lines.append(f"# Auto-converted from {nb_path.name}")
    lines.append("\n")
    for i, cell in enumerate(data.get("cells", []), start=1):
        cell_type = cell.get("cell_type", "unknown")
        lines.append(f"# %% [cell {i}] type={cell_type}")
        src = "".join(cell.get("source", []))
        if cell_type == "code":
            lines.append(src.rstrip())
        else:
            for ln in src.splitlines():
                lines.append(f"# {ln}")
        lines.append("\n")
    out_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> None:
    code_dir = Path("Code")
    out_dir = Path("scripts/notebook_converted")
    out_dir.mkdir(parents=True, exist_ok=True)

    notebooks = sorted([p for p in code_dir.glob("*.ipynb") if not p.name.startswith(".")])
    for nb in notebooks:
        out = out_dir / f"{nb.stem}.py"
        notebook_to_py(nb, out)
        print(f"Wrote: {out}")


if __name__ == "__main__":
    main()
