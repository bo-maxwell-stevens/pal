# Auto-converted from 02_fetch_data.ipynb


# %% [cell 1] type=markdown
# # 02_fetch_data
# 
# Fetch a small root-only subset (100 samples) from the GlobalAMFungi archive.
# 
# This notebook uses the public Zenodo archive for `globalbioticinteractions/globalamfungi` and writes a filtered CSV into `../Data`.


# %% [cell 2] type=code
from pathlib import Path
from urllib.request import urlretrieve
import csv
import io
import zipfile

ZENODO_ZIP_URL = "https://zenodo.org/api/records/18424800/files/data.zip/content"
TARGET_N = 100

code_dir = Path.cwd()
data_dir = (code_dir / "../Data").resolve()
data_dir.mkdir(parents=True, exist_ok=True)

zip_path = data_dir / "globalamfungi_data.zip"
out_path = data_dir / "globalamfungi_roots_100.csv"

print(f"Data dir: {data_dir}")
print(f"Zip path: {zip_path}")
print(f"Output path: {out_path}")


# %% [cell 3] type=code
if not zip_path.exists():
    print("Downloading GlobalAMFungi archive (one-time)...")
    urlretrieve(ZENODO_ZIP_URL, zip_path)
    print("Download complete.")
else:
    print("Archive already present, skipping download.")

print(f"Archive size: {zip_path.stat().st_size / (1024**2):.1f} MB")


# %% [cell 4] type=code
def find_main_table_member(zf: zipfile.ZipFile) -> str:
    candidates = [n for n in zf.namelist() if not n.endswith('/')]

    for name in candidates:
        with zf.open(name, "r") as fh:
            header = fh.readline().decode("utf-8", errors="ignore").strip()

        if "sample_type" in header and "plants_dominant" in header and "id" in header:
            return name

    raise RuntimeError("Could not find the main GlobalAMFungi table in the ZIP.")

with zipfile.ZipFile(zip_path, "r") as zf:
    main_member = find_main_table_member(zf)

print(f"Main table member: {main_member}")


# %% [cell 5] type=code
selected = []

with zipfile.ZipFile(zip_path, "r") as zf:
    with zf.open(main_member, "r") as raw:
        text_stream = io.TextIOWrapper(raw, encoding="utf-8", errors="replace", newline="")
        reader = csv.DictReader(text_stream, delimiter='\t')

        for row in reader:
            row = {k: v for k, v in row.items() if k and k.strip()}
            sample_type = (row.get("sample_type") or "").strip().lower()
            if "root" in sample_type:
                selected.append(row)

            if len(selected) >= TARGET_N:
                break

if len(selected) < TARGET_N:
    raise RuntimeError(f"Only found {len(selected)} root rows, expected at least {TARGET_N}.")

fieldnames = list(selected[0].keys())
with out_path.open("w", encoding="utf-8", newline="") as out_fh:
    writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(selected)

print(f"Wrote {len(selected)} root rows to: {out_path}")


# %% [cell 6] type=code
with out_path.open("r", encoding="utf-8", newline="") as fh:
    reader = csv.DictReader(fh)
    rows = list(reader)

print("Rows in output:", len(rows))
assert len(rows) == TARGET_N, f"Expected {TARGET_N}, got {len(rows)}"

sample_types = sorted({(r.get('sample_type') or '').strip() for r in rows})
print("Sample types present:", sample_types)

print("\nPreview (first 5 rows):")
for i, r in enumerate(rows[:5], start=1):
    print({
        'id': r.get('id'),
        'sample_type': r.get('sample_type'),
        'plants_dominant': r.get('plants_dominant'),
        'latitude': r.get('latitude'),
        'longitude': r.get('longitude'),
        'doi': r.get('doi')
    })
