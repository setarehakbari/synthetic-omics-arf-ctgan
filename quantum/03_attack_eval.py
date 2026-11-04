# exp_real/03_attack_eval.py
# Evaluate adversarial robustness on REAL test set:
# - loads X_test / y_test (CSV or NPY)
# - applies the same StandardScaler used at training (if available)
# - rebuilds MLP EXACTLY to match the checkpoint (asymmetric hidden widths)
# - loads PQCWithProj (proj/head [+qclf if present in state])
# - crafts FGSM adversarial examples on MLP and measures transfer to PQC
# - saves clean accuracies + FGSM sweep table and an accuracy-vs-epsilon plot

import warnings, json, re, pickle
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
import joblib

HERE = Path(__file__).resolve().parent
OUT  = HERE / "outputs"
(OUT/"figs").mkdir(parents=True, exist_ok=True)
(OUT/"tables").mkdir(parents=True, exist_ok=True)

# ---------------- I/O helpers ----------------
def load_array(stem: str):
    """
    Load array for 'X_train/X_test/y_train/y_test' from NPY if present, else CSV.
    For y*, supports headered ('label','y','target') or single-column CSV.
    """
    npy = OUT / f"{stem}.npy"
    csv = OUT / f"{stem}.csv"
    if npy.exists():
        return np.load(npy, allow_pickle=False)
    if csv.exists():
        if stem.startswith("X_"):
            return pd.read_csv(csv).to_numpy()
        df = pd.read_csv(csv)
        if df.shape[1] == 1:
            return df.iloc[:, 0].astype(str).to_numpy()
        for c in ["label", "y", "target"]:
            if c in df.columns:
                return df[c].astype(str).to_numpy()
        return df.iloc[:, 0].astype(str).to_numpy()
    raise FileNotFoundError(f"Missing both {npy} and {csv}")

def encode_labels_like_train(y_str, mapping=None):
    """
    Map string labels to ints. If HC/MDD present, uses {'HC':0,'MDD':1},
    else factorizes in sorted order. Returns (y_int, mapping).
    """
    y_str = pd.Series(y_str).astype(str)
    if mapping is None:
        if {"HC","MDD"}.issubset(set(y_str.unique())):
            mapping = {"HC": 0, "MDD": 1}
        else:
            labs = sorted(y_str.unique())
            mapping = {lab: i for i, lab in enumerate(labs)}
    y = y_str.map(mapping)
    if y.isna().any():
        cats = sorted(y_str[y.isna()].unique())
        warnings.warn(f"[warn] Unrecognized labels {cats}; factorizing in sorted order.")
        labs = sorted(y_str.unique())
        mapping = {lab: i for i, lab in enumerate(labs)}
        y = y_str.map(mapping)
    return y.astype(int).to_numpy(), mapping

# ---------------- Models ----------------
class AsymMLP(nn.Module):
    """
    MLP with two possibly different hidden widths (h1 != h2).
    Matches a typical Sequential:
      net.0 Linear(in_dim, h1)
      net.3 Linear(h1,    h2)
      net.6 Linear(h2,    out_dim)
    """
    def __init__(self, in_dim, h1, h2, out_dim, p=0.2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, h1), nn.ReLU(), nn.Dropout(p),
            nn.Linear(h1,   h2),   nn.ReLU(), nn.Dropout(p),
            nn.Linear(h2, out_dim)
        )
    def forward(self, x): 
        return self.net(x)

def infer_mlp_dims_from_state(ckpt: dict):
    """
    Infer (in_dim, h1, h2, out_dim) from a Sequential with keys like:
      net.0.weight -> (h1, in_dim)
      net.3.weight -> (h2, h1)
      net.6.weight -> (out_dim, h2)
    Returns None if shapes not found.
    """
    try:
        w0 = ckpt["net.0.weight"]  # (h1, in_dim)
        w1 = ckpt["net.3.weight"]  # (h2, h1)
        w2 = ckpt["net.6.weight"]  # (out_dim, h2)
        h1, in_dim  = int(w0.shape[0]), int(w0.shape[1])
        h2, _h1c    = int(w1.shape[0]), int(w1.shape[1])
        out_dim, h2c= int(w2.shape[0]), int(w2.shape[1])
        # quick consistency check
        assert _h1c == h1 and h2c == h2
        return in_dim, h1, h2, out_dim
    except Exception:
        return None

class PQCWithProj(nn.Module):
    """
    Projection + (optional QClassifier) + head, mirroring the training wrapper:
      proj: R^D -> R^(n_qubits)
      head: R^(n_qubits) -> R^C
    Inference path is robust: head(proj(x)).
    """
    def __init__(self, in_dim, out_dim, n_qubits):
        super().__init__()
        self.proj = nn.Linear(in_dim, n_qubits)
        self.head = nn.Linear(n_qubits, out_dim)
        self.qclf = None  # left None unless custom quantum restore is added
    def forward(self, x):
        z = self.proj(x)
        return self.head(z)

# ---------------- Load data ----------------
Xte_np = load_array("X_test")
yte_str = load_array("y_test")

# Apply the same StandardScaler as in training (if present)
scaler_path = OUT / "scaler.pkl"
if scaler_path.exists():
    try:
        scaler = joblib.load(scaler_path)
    except Exception:
        with open(scaler_path, "rb") as f:
            scaler = pickle.load(f)
    Xte_np = scaler.transform(Xte_np)
else:
    warnings.warn("[warn] scaler.pkl not found; evaluating on unscaled features.")

# Discover sizes from manifest if present
n_feat = int(Xte_np.shape[1])
manifest = OUT / "train_manifest.json"
n_classes = 2
n_qubits  = 6
if manifest.exists():
    mf = json.loads(manifest.read_text())
    n_classes = int(mf.get("n_classes", 2))
    n_qubits  = int(mf.get("n_qubits", 6))

yte, lab_map = encode_labels_like_train(yte_str)

# Basic shape/type prep
Xte = torch.tensor(np.asarray(Xte_np, dtype=np.float32))
yte = torch.tensor(np.asarray(yte,    dtype=np.int64))

# Safety: ensure non-empty & matching lengths
n = min(len(Xte), len(yte))
if n != len(Xte) or n != len(yte):
    warnings.warn(f"[warn] length mismatch X:{len(Xte)} vs y:{len(yte)} → truncating to {n}.")
Xte, yte = Xte[:n], yte[:n]

print(f"[info] X_test shape: {tuple(Xte.shape)} | y_test shape: {tuple(yte.shape)} | y classes: {sorted(set(yte.tolist()))}")

# ---------------- Load trained models ----------------
# MLP — rebuild EXACT architecture from checkpoint (asymmetric hidden widths)
mlp_path = OUT / "mlp.pt"
if mlp_path.exists():
    ckpt = torch.load(mlp_path, map_location="cpu")
    dims = infer_mlp_dims_from_state(ckpt)
    if dims is None:
        # Fallback: try symmetric 128-128 (won't match your checkpoint but keeps code running)
        warnings.warn("[warn] Could not infer MLP dims from checkpoint; using fallback 128-128.")
        mlp = AsymMLP(n_feat, 128, 128, n_classes)
        mlp.load_state_dict(ckpt, strict=False)
    else:
        in_dim, h1, h2, out_dim = dims
        # Sanity print to make debugging easy
        print(f"[mlp] inferred dims from ckpt: in={in_dim}, h1={h1}, h2={h2}, out={out_dim}")
        # Build EXACT same architecture then strict-load
        mlp = AsymMLP(in_dim, h1, h2, out_dim)
        mlp.load_state_dict(ckpt, strict=True)
        # If in/out differ from current test dims/classes, we still can run:
        # - If in_dim != n_feat: a scaler/projection mismatch exists; ideally use the same features.
        if in_dim != n_feat:
            warnings.warn(f"[warn] Test features ({n_feat}) differ from train MLP input ({in_dim}). "
                          f"Ensure you evaluate on the same selected features as training.")
        if out_dim != n_classes:
            warnings.warn(f"[warn] Test n_classes ({n_classes}) differs from MLP out_dim ({out_dim}). "
                          f"Predictions will still work but class mapping may differ.")
else:
    warnings.warn("[warn] Missing mlp.pt; MLP will be randomly initialized (not recommended).")
    mlp = AsymMLP(n_feat, 256, 128, n_classes)
mlp.eval()

# PQC wrapper — restore proj/head if present
pqc_path = OUT / "pqc.pt"
pqc = PQCWithProj(in_dim=n_feat, out_dim=n_classes, n_qubits=n_qubits)
if pqc_path.exists():
    state = torch.load(pqc_path, map_location="cpu")
    if isinstance(state, dict):
        if "proj" in state: pqc.proj.load_state_dict(state["proj"])
        if "head" in state: pqc.head.load_state_dict(state["head"])
else:
    warnings.warn("[warn] Missing pqc.pt; PQC will be randomly initialized (not recommended).")
pqc.eval()

# ---------------- Metrics helpers ----------------
@torch.no_grad()
def acc(model, x, y):
    """Compute standard accuracy; returns NaN if x is empty."""
    if x.nelement() == 0:
        return float("nan")
    return (model(x).argmax(1) == y).float().mean().item()

def fgsm(model, x, y, eps=0.1):
    """
    Single-step FGSM crafted on 'model' with label 'y'.
    Gradient sign is taken on input; result is clipped to a reasonable range.
    """
    x = x.clone().detach().requires_grad_(True)
    loss = nn.CrossEntropyLoss()(model(x), y)
    loss.backward()
    with torch.no_grad():
        x_adv = x + eps * x.grad.sign()
    return torch.clamp(x_adv, -5, 5)

# ---------------- Clean accuracies (no attack) ----------------
clean_mlp = acc(mlp, Xte, yte)
clean_pqc = acc(pqc, Xte, yte)
pd.DataFrame({"mlp_clean": [clean_mlp], "pqc_clean": [clean_pqc]}).to_csv(
    OUT / "tables" / "real_clean_acc.csv", index=False
)
print(f"[clean] MLP={clean_mlp:.3f} | PQC={clean_pqc:.3f}")

# ---------------- FGSM sweep (crafted on MLP; transfer to PQC) ----------------
eps_list = [0.00, 0.03, 0.06, 0.10, 0.20]
rows = []
for eps in eps_list:
    x_adv = fgsm(mlp, Xte, yte, eps=eps)  # craft on MLP
    rows.append({
        "eps": eps,
        "acc_mlp": acc(mlp, x_adv, yte),
        "acc_pqc": acc(pqc, x_adv, yte)   # transfer to PQCWithProj
    })

tab = pd.DataFrame(rows)
tab_path = OUT / "tables" / "real_fgsm_transfer.csv"
tab.to_csv(tab_path, index=False)
print(tab)

# ---------------- Plot: accuracy vs epsilon ----------------
# # =========================================================
# Plot: adversarial accuracy vs FGSM perturbation magnitude
# =========================================================

import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 15

fig, ax = plt.subplots(figsize=(6.2, 4.5))

# MLP: solid orange line
ax.plot(
    tab["eps"],
    tab["acc_mlp"],
    linestyle="-",
    linewidth=2.2,
    color="#E69F00",
    label="MLP (attacked)"
)

# PQC: dashed green line
ax.plot(
    tab["eps"],
    tab["acc_pqc"],
    linestyle="--",
    linewidth=2.2,
    color="#009E73",
    label="PQC (transfer)"
)

# Axis labels
ax.set_xlabel("FGSM perturbation magnitude (ε)", fontsize=15, fontname="Times New Roman")
ax.set_ylabel("Accuracy", fontsize=15, fontname="Times New Roman")

# Axis limits and ticks
ax.set_xlim(-0.002, 0.205)
ax.set_ylim(0, 1.0)
ax.tick_params(axis="both", labelsize=15)

# Legend (upper right)
leg = ax.legend(
    loc="upper right",
    frameon=True,
    fontsize=15
)
for text in leg.get_texts():
    text.set_fontname("Times New Roman")

# No title
ax.set_title("")

plt.tight_layout()
plt.savefig("adv_accuracy_plot.png", dpi=300, bbox_inches="tight")
plt.show()
