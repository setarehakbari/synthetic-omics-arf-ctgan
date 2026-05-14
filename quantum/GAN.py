import os, time
import pandas as pd
from sdv.metadata import Metadata
from sdv.single_table import CTGANSynthesizer

# ---- Paths ----
base_dir = os.path.expanduser("~/Desktop/My paper")
in_csv   = os.path.join(base_dir, "gene_with_label.csv")
out_csv  = os.path.join(base_dir, "synthetic_cgan.csv")  # for R
meta_json = os.path.join(base_dir, "cgan_metadata.json")


# ---- 1) Load data ----
df = pd.read_csv(in_csv)

# ---- 2) Build metadata ----
print("Building metadata ...")
meta = Metadata.detect_from_dataframe(df)

# Explicitly set 'class' as categorical
if 'class' in df.columns:
    meta.update_column('class', sdtype='categorical')

# Save metadata JSON
if os.path.exists(meta_json):
    os.remove(meta_json)
meta.save_to_json(meta_json)
print(f"Metadata saved to: {meta_json}")

# ---- 3) Initialize CTGAN ----
print("Initializing CTGANSynthesizer ...")
syn = CTGANSynthesizer(
    metadata=meta,
    epochs=50,        # for testing; increase to ~300 for real runs
    batch_size=256,
    pac=1,
    verbose=True,
    # device="cpu"   # or "cuda" if you have GPU
)

# ---- 4) Fit ----
print("Fitting (this can take a while) ...")
t0 = time.time()
syn.fit(df)
print(f"Fit done in {time.time()-t0:.1f}s")

# ---- 5) Sample ----
print("Sampling synthetic data ...")
syn_all = syn.sample(num_rows=200)
print("Sample shape =", syn_all.shape)

# ---- 6) Save output ----
syn_all.to_csv(out_csv, index=False)
print(" Synthetic data saved to:", out_csv)

# ---- 7) Check class balance ----
print("Real class balance:")
print(df['class'].value_counts(normalize=True).round(3))

print("Synthetic class balance:")
print(syn_all['class'].value_counts(normalize=True).round(3))
