

# exp_real/02_train_models.py  (stronger baseline + scaler + clean metrics)
import os, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import torch, torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader

HERE = Path(__file__).resolve().parent
OUT  = HERE / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

# -------------------- I/O helpers --------------------
def load_array(stem):
    npy = OUT / f"{stem}.npy"
    csv = OUT / f"{stem}.csv"
    if npy.exists():
        return np.load(npy, allow_pickle=False)
    if csv.exists():
        if stem.startswith("X_"):
            return pd.read_csv(csv).to_numpy()
        df = pd.read_csv(csv)
        if df.shape[1] == 1:
            return df.iloc[:,0].astype(str).to_numpy()
        for c in ["label","y","target"]:
            if c in df.columns:
                return df[c].astype(str).to_numpy()
        return df.iloc[:,0].astype(str).to_numpy()
    raise FileNotFoundError(f"[load_array] Missing both {npy} and {csv}")

def encode_labels(y_str):
    y_str = pd.Series(y_str).astype(str)
    if set(["HC","MDD"]).issubset(set(y_str.unique())):
        mapping = {"HC":0, "MDD":1}
    else:
        labs = sorted(y_str.unique())
        mapping = {lab:i for i,lab in enumerate(labs)}
    print("[y-map]", mapping)
    return y_str.map(mapping).astype(int).to_numpy(), len(mapping)

print("[train] Loading CSV/NPY splits from exp_real/outputs/ ...")
Xtr = load_array("X_train"); ytr_str = load_array("y_train")
Xte = load_array("X_test");  yte_str = load_array("y_test")

ytr, n_cls_tr = encode_labels(ytr_str)
yte, n_cls_te = encode_labels(yte_str)
n_cls = max(n_cls_tr, n_cls_te)

Xtr = np.asarray(Xtr, dtype=np.float32)
Xte = np.asarray(Xte, dtype=np.float32)
ytr = np.asarray(ytr, dtype=np.int64)
yte = np.asarray(yte, dtype=np.int64)
n_feat = Xtr.shape[1]
print(f"[train] Xtr={Xtr.shape} Xte={Xte.shape} | ytr={ytr.shape} yte={yte.shape} | n_cls={n_cls}")

# -------------------- Standardize features --------------------
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler().fit(Xtr)
Xtr_s = scaler.transform(Xtr)
Xte_s = scaler.transform(Xte)

# persist scaler for eval/attacks
try:
    import joblib, pickle
    joblib.dump(scaler, OUT/"scaler.pkl")
except Exception:
    # fallback to pickle
    import pickle
    with open(OUT/"scaler.pkl", "wb") as f:
        pickle.dump(scaler, f)
print(f"[train] Saved scaler → {OUT/'scaler.pkl'}")

# -------------------- Random Forest (baseline) --------------------
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

rf = RandomForestClassifier(n_estimators=1200, max_features="sqrt",
                            min_samples_leaf=2, n_jobs=-1, random_state=42)
rf.fit(Xtr_s, ytr)
rf_acc = accuracy_score(yte, rf.predict(Xte_s))
pd.DataFrame({"rf_test_accuracy":[rf_acc]}).to_csv(OUT/"rf_test_metrics.csv", index=False)
print(f"[RF] test acc = {rf_acc:.3f} | saved → {OUT/'rf_test_metrics.csv'}")

# -------------------- Stronger MLP --------------------
class MLP(nn.Module):
    def __init__(self, in_dim, out_dim, h1=256, h2=128, p=0.25):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, h1), nn.ReLU(), nn.Dropout(p),
            nn.Linear(h1, h2),    nn.ReLU(), nn.Dropout(p),
            nn.Linear(h2, out_dim)
        )
    def forward(self, x): return self.net(x)

def train_mlp(X, y, epochs=120, bs=32, lr=1e-3, wd=1e-4):
    mlp = MLP(n_feat, n_cls, 256, 128, 0.25)
    opt = torch.optim.Adam(mlp.parameters(), lr=lr, weight_decay=wd)
    crit= nn.CrossEntropyLoss()
    dl  = DataLoader(TensorDataset(torch.tensor(X), torch.tensor(y)),
                     batch_size=bs, shuffle=True)
    mlp.train()
    for ep in range(1, epochs+1):
        s=0.0
        for xb, yb in dl:
            xb=xb.float(); yb=yb.long()
            opt.zero_grad()
            loss = crit(mlp(xb), yb); loss.backward(); opt.step()
            s += loss.item()*len(xb)
        if ep==1 or ep%20==0 or ep==epochs:
            print(f"[MLP] epoch {ep:3d}/{epochs} | loss={s/len(dl.dataset):.4f}")
    torch.save(mlp.state_dict(), OUT/"mlp.pt")
    print(f"[MLP] Saved → {OUT/'mlp.pt'}")
    return mlp

mlp = train_mlp(Xtr_s, ytr)

# -------------------- PQC wrapper (with projection) --------------------
HAS_QML = True
try:
    from exp_real_02_utils import QClassifier
except Exception as e:
    HAS_QML = False
    warnings.warn(f"[warn] QClassifier not found ({e}); using classical-only head.")

def get_n_qubits_default():
    if HAS_QML and hasattr(QClassifier, "N_QUBITS"):
        try: return int(getattr(QClassifier, "N_QUBITS"))
        except Exception: pass
    return 6

N_QUBITS = get_n_qubits_default()
print(f"[train] Using N_QUBITS={N_QUBITS}")

class PQCWithProj(nn.Module):
    def __init__(self, in_dim, out_dim, n_qubits):
        super().__init__()
        self.proj = nn.Linear(in_dim, n_qubits)
        self.head = nn.Linear(n_qubits, out_dim)
        self.qclf = None
    def forward(self, x):
        z = self.proj(x)
        return self.head(z)

def train_pqc(X, y, epochs=25, bs=32, lr=1e-3):
    model = PQCWithProj(n_feat, n_cls, N_QUBITS)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    crit= nn.CrossEntropyLoss()
    dl  = DataLoader(TensorDataset(torch.tensor(X), torch.tensor(y)),
                     batch_size=bs, shuffle=True)
    model.train()
    for ep in range(1, epochs+1):
        s=0.0
        for xb, yb in dl:
            xb=xb.float(); yb=yb.long()
            opt.zero_grad(); loss = crit(model(xb), yb)
            loss.backward(); opt.step()
            s += loss.item()*len(xb)
        if ep==1 or ep%5==0 or ep==epochs:
            print(f"[PQC] epoch {ep:2d}/{epochs} | loss={s/len(dl.dataset):.4f}")
    state = {
        "proj": model.proj.state_dict(),
        "head": model.head.state_dict(),
        "n_qubits": N_QUBITS,
        "n_feat": n_feat,
        "n_cls": n_cls
    }
    torch.save(state, OUT/"pqc.pt")
    print(f"[PQC] Saved → {OUT/'pqc.pt'}")
    return model

pqc = train_pqc(Xtr_s, ytr)

# -------------------- Clean test metrics (RF/MLP/PQC) --------------------
@torch.no_grad()
def acc_nn(model, X, y):
    model.eval()
    X = torch.tensor(X, dtype=torch.float32)
    y = torch.tensor(y, dtype=torch.long)
    return (model(X).argmax(1) == y).float().mean().item()

clean_mlp = acc_nn(mlp, Xte_s, yte)
clean_pqc = acc_nn(pqc, Xte_s, yte)

clean_tbl = pd.DataFrame({
    "rf_test_acc":[rf_acc],
    "mlp_test_acc":[clean_mlp],
    "pqc_test_acc":[clean_pqc]
})
clean_tbl.to_csv(OUT/"clean_test_metrics.csv", index=False)
print(f"[clean] RF={rf_acc:.3f} | MLP={clean_mlp:.3f} | PQC={clean_pqc:.3f}")
print(f"[train] Wrote manifest + metrics in: {OUT}")
with open(OUT/"train_manifest.json","w") as f:
    json.dump({
        "n_train": int(Xtr.shape[0]),
        "n_test" : int(Xte.shape[0]),
        "n_features": int(n_feat),
        "n_classes": int(n_cls),
        "n_qubits": int(N_QUBITS),
        "files": {
            "rf_metrics": str(OUT/"rf_test_metrics.csv"),
            "mlp": str(OUT/"mlp.pt"),
            "pqc": str(OUT/"pqc.pt"),
            "scaler": str(OUT/"scaler.pkl"),
            "clean_metrics": str(OUT/"clean_test_metrics.csv"),
        }
    }, f, indent=2)


