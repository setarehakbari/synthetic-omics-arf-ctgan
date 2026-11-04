# exp_transfer/01_align_data.py
# Align synthetic data to real Top-K features and z-score with real μ,σ.

import pandas as pd, numpy as np
from pathlib import Path

OUT = Path("../outputs"); (OUT/"tables").mkdir(parents=True, exist_ok=True)
mu = np.load(OUT/"real_mu.npy"); sd = np.load(OUT/"real_sd.npy")
topk = pd.read_csv(OUT/"tables/real_topK_features.csv")["feature"].tolist()

s_arf   = pd.read_csv("../data/synth_arf.csv")
s_ctgan = pd.read_csv("../data/synth_ctgan.csv")

def align_and_z(df):
    cols = [c for c in topk if c in df.columns]
    Z = (df[cols] - mu[:len(cols)]) / sd[:len(cols)]
    return Z.values

np.save(OUT/"X_synth_arf.npy",   align_and_z(s_arf))
np.save(OUT/"X_synth_ctgan.npy", align_and_z(s_ctgan))
print("Aligned ARF/CTGAN to real Top-K and z-scored.")
