# exp_transfer/02_eval_transfer.py
# Compare PQC/MLP confidence on real-adv vs synthetic (clean) inputs.

import numpy as np, torch, pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from exp_real_02_utils import MLP, QClassifier
import torch.nn as nn

OUT = Path("../outputs"); (OUT/"figs").mkdir(parents=True, exist_ok=True)
Xte = torch.tensor(np.load(OUT/"X_test.npy"),  dtype=torch.float32)
yte = torch.tensor(np.load(OUT/"y_test.npy"),  dtype=torch.long)
X_arf   = torch.tensor(np.load(OUT/"X_synth_arf.npy"),   dtype=torch.float32)
X_ctgan = torch.tensor(np.load(OUT/"X_synth_ctgan.npy"), dtype=torch.float32)

mlp = MLP(Xte.shape[1]); mlp.load_state_dict(torch.load(OUT/"mlp.pt")); mlp.eval()
pqc = QClassifier(); pqc.load_state_dict(torch.load(OUT/"pqc.pt")); pqc.eval()

def fgsm(model, x, y, eps=0.1):
    x_adv = x.clone().detach().requires_grad_(True)
    loss = nn.CrossEntropyLoss()(model(x_adv), y)
    loss.backward()
    with torch.no_grad():
        x_adv = x_adv + eps * x_adv.grad.sign()
    return torch.clamp(x_adv, -5, 5)

def mean_confidence(m, X):
    with torch.no_grad():
        p = torch.softmax(m(X), dim=1).max(1).values
    return float(p.mean().item())

eps_list = [0.0, 0.03, 0.06, 0.1, 0.2]
rows = []
for eps in eps_list:
    X_adv_real = fgsm(mlp, Xte, yte, eps=eps)
    rows.append({
        "eps": eps,
        "conf_pqc_real": mean_confidence(pqc, X_adv_real),
        "conf_pqc_arf":  mean_confidence(pqc, X_arf),
        "conf_pqc_ctg":  mean_confidence(pqc, X_ctgan),
        "conf_mlp_real": mean_confidence(mlp, X_adv_real),
        "conf_mlp_arf":  mean_confidence(mlp, X_arf),
        "conf_mlp_ctg":  mean_confidence(mlp, X_ctgan),
    })

tab = pd.DataFrame(rows)
tab.to_csv(OUT/"tables/transfer_confidence_vs_eps.csv", index=False)
print(tab)

plt.figure(figsize=(7,4))
plt.plot(tab["eps"], tab["conf_pqc_real"], marker="o", label="PQC on Real (adv)")
plt.plot(tab["eps"], tab["conf_pqc_arf"],  marker="o", label="PQC on ARF (clean)")
plt.plot(tab["eps"], tab["conf_pqc_ctg"],  marker="o", label="PQC on CTGAN (clean)")
plt.xlabel("FGSM ε (crafted on real)"); plt.ylabel("Mean max-probability"); plt.ylim(0,1); plt.legend()
plt.title("Transfer — PQC confidence under FGSM crafted on Real")
plt.tight_layout(); plt.savefig(OUT/"figs/Transfer_PQC_confidence.png", dpi=300)
