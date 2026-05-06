# Auto-converted from Untitled.ipynb


# %% [cell 1] type=code
import pandas as pd

url = "https://github.com/globalbioticinteractions/globalamfungi/raw/main/globalamfungi.tsv.gz"
df = pd.read_csv(url, sep="\t", compression="gzip", low_memory=False)
print(df.shape)
print(df.columns.tolist()[:50])


# %% [cell 2] type=code
df.to_csv('../Data/globalamf.csv')


# %% [cell 3] type=code
set(df['sample_type'])


# %% [cell 4] type=code
root = df[df['sample_type'] == 'root'].copy()

root['n_plants'] = root['plants_dominant'].fillna('').apply(
    lambda x: len([p for p in x.split(';') if p.strip()])
)

print(root['n_plants'].describe())

print("\nProportion single plant:", (root['n_plants'] == 1).mean())
print("Proportion multiple plants:", (root['n_plants'] > 1).mean())


# %% [cell 5] type=code
root = df[df['sample_type'] == 'root'].copy()

root['plants_dominant'] = root['plants_dominant'].str.strip()

root['n_plants'] = root['plants_dominant'].apply(
    lambda x: len(x.split(';')) if pd.notna(x) else 0
)

strict = root[
    (root['n_plants'] == 1) &
    (root['plants_dominant'].notna())
].copy()

print("Samples retained:", strict['id'].nunique())
