
# exp_real/01_prep.py
import json, sys
from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold

HERE = Path(__file__).parent
DATA = HERE.parent / "data" / "real.csv"
OUT  = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

print(f"[prep] Reading: {DATA}")
df = pd.read_csv(DATA)

# infer label column
LABEL_CANDIDATES = ["label", "class", "y", "target"]
label_col = None
for c in LABEL_CANDIDATES:
    if c in df.columns:
        label_col = c
        break
if label_col is None:
    raise RuntimeError(f"No label column found. Expected one of {LABEL_CANDIDATES}")

# rename to 'label' internally
if label_col != "label":
    df = df.rename(columns={label_col: "label"})

# ensure label is categorical string (HC/MDD-like if 0/1)
if pd.api.types.is_numeric_dtype(df["label"]):
    df["label"] = df["label"].map({0: "HC", 1: "MDD"}).astype("category")
else:
    df["label"] = df["label"].astype("category")

X = df.drop(columns=["label"])
y = df["label"]

# drop constant features
selector = VarianceThreshold(threshold=0.0)
X_sel = pd.DataFrame(selector.fit_transform(X), columns=X.columns[selector.get_support(indices=True)])

# stratified split
Xtr, Xte, ytr, yte = train_test_split(X_sel, y, test_size=0.25, random_state=42, stratify=y)

# save CSV splits
(Xtr).to_csv(OUT / "X_train.csv", index=False)
(Xte).to_csv(OUT / "X_test.csv",  index=False)
ytr.to_frame("label").to_csv(OUT / "y_train.csv", index=False)
yte.to_frame("label").to_csv(OUT / "y_test.csv",  index=False)

# manifest
manifest = {
    "data": str(DATA),
    "out_dir": str(OUT),
    "n_train": int(Xtr.shape[0]),
    "n_test": int(Xte.shape[0]),
    "n_features": int(Xtr.shape[1]),
    "label_name_original": label_col
}
with open(OUT / "prep_manifest.json", "w") as f:
    json.dump(manifest, f, indent=2)

print("[prep] Done. Wrote splits to:", OUT)
