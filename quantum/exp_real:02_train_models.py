# exp_real/02_train_models.py
# Train classical baselines on real-data splits produced by 01_prep.py.
# Reads CSVs (not NPY), saves models & metrics into exp_real/outputs/.

import os, json, sys
import pandas as pd
import numpy as np

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT  = os.path.join(os.path.dirname(__file__), "outputs")
os.makedirs(OUT, exist_ok=True)

# --------- Helper: safe imports ----------
have_torch = False
try:
    import torch
    have_torch = True
except Exception:
    pass

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score, f1_score
import joblib

def load_csv_split():
    xtr = os.path.join(OUT, "X_train.csv")
    ytr = os.path.join(OUT, "y_train.csv")
    xte = os.path.join(OUT, "X_test.csv")
    yte = os.path.join(OUT, "y_test.csv")

    for p in [xtr, ytr, xte, yte]:
        if not os.path.exists(p):
            print(f"[train][ERROR] Missing file: {p}", file=sys.stderr)
            sys.exit(1)

    X_train = pd.read_csv(xtr)
    y_train = pd.read_csv(ytr).iloc[:,0]
    X_test  = pd.read_csv(xte)
    y_test  = pd.read_csv(yte).iloc[:,0]

    # Ensure numeric dtypes
    X_train = X_train.apply(pd.to_numeric, errors="coerce")
    X_test  = X_test.apply(pd.to_numeric,  errors="coerce")

    return X_train, y_train, X_test, y_test

def eval_metrics(y_true, y_prob_or_label, is_proba):
    if is_proba:
        y_prob = y_prob_or_label
        y_pred = (y_prob >= 0.5).astype(int)
    else:
        y_pred = y_prob_or_label
        # best-effort probas for metrics uniform:
        y_prob = None

    acc = accuracy_score(y_true, y_pred)
    f1  = f1_score(y_true, y_pred, average="macro")
    try:
        auc = roc_auc_score(y_true, y_pred if y_prob is None else y_prob)
    except Exception:
        auc = float("nan")
    return {"accuracy": acc, "macro_f1": f1, "auc": auc}

def train_sklearn_models(X_train, y_train, X_test, y_test):
    results = {}

    # Random Forest
    rf = RandomForestClassifier(
        n_estimators=1000,
        max_features="sqrt",
        n_jobs=-1,
        random_state=42
    )
    rf.fit(X_train, y_train)
    rf_proba = rf.predict_proba(X_test)[:,1] if len(np.unique(y_train))==2 else None
    rf_pred  = rf.predict(X_test)
    results["random_forest"] = eval_metrics(y_test, rf_proba if rf_proba is not None else rf_pred, is_proba=(rf_proba is not None))
    joblib.dump(rf, os.path.join(OUT, "model_rf.joblib"))

    # Logistic Regression (baseline)
    lr = LogisticRegression(max_iter=2000, n_jobs=-1)
    lr.fit(X_train, y_train)
    lr_proba = lr.predict_proba(X_test)[:,1] if len(np.unique(y_train))==2 else None
    lr_pred  = lr.predict(X_test)
    results["log_reg"] = eval_metrics(y_test, lr_proba if lr_proba is not None else lr_pred, is_proba=(lr_proba is not None))
    joblib.dump(lr, os.path.join(OUT, "model_lr.joblib"))

    return results

def maybe_torch_mlp(X_train, y_train, X_test, y_test):
    # Optional tiny MLP if torch is available
    if not have_torch:
        return None

    import torch.nn as nn
    import torch.optim as optim

    Xtr = torch.tensor(X_train.values, dtype=torch.float32)
    ytr = torch.tensor(y_train.values, dtype=torch.long)
    Xte = torch.tensor(X_test.values,  dtype=torch.float32)
    yte = torch.tensor(y_test.values,  dtype=torch.long)

    in_dim = Xtr.shape[1]
    hidden = 64
    out_dim = len(np.unique(y_train))

    class MLP(nn.Module):
        def __init__(self):
            super().__init__()
            self.net = nn.Sequential(
                nn.Linear(in_dim, hidden),
                nn.ReLU(),
                nn.Linear(hidden, out_dim)
            )
        def forward(self, x):
            return self.net(x)

    model = MLP()
    opt = optim.Adam(model.parameters(), lr=1e-3)
    crit = nn.CrossEntropyLoss()

    model.train()
    for epoch in range(50):
        opt.zero_grad()
        logits = model(Xtr)
        loss = crit(logits, ytr)
        loss.backward()
        opt.step()

    model.eval()
    with torch.no_grad():
        logits = model(Xte)
        prob = torch.softmax(logits, dim=1)[:,1].cpu().numpy() if out_dim==2 else None
        pred = torch.argmax(logits, dim=1).cpu().numpy()

    metrics = eval_metrics(y_test.values, prob if prob is not None else pred, is_proba=(prob is not None))
    # Save torch model
    torch.save(model.state_dict(), os.path.join(OUT, "model_mlp.pt"))
    return {"mlp_torch": metrics}

def main():
    print("[train] Loading CSV splits from exp_real/outputs/ ...")
    X_train, y_train, X_test, y_test = load_csv_split()
    print(f"[train] Shapes: X_train={X_train.shape}, X_test={X_test.shape}")

    print("[train] Training sklearn models ...")
    results = train_sklearn_models(X_train, y_train, X_test, y_test)

    print("[train] (Optional) Training small torch MLP ...")
    mlp_res = maybe_torch_mlp(X_train, y_train, X_test, y_test)
    if mlp_res is not None:
        results.update(mlp_res)

    # Save metrics
    metrics_path = os.path.join(OUT, "metrics_train.json")
    with open(metrics_path, "w") as f:
        json.dump(results, f, indent=2)

    print("[train] Saved models & metrics to exp_real/outputs/")
    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
