# exp_real_02_utils.py
# Shared model definitions used by exp_real and exp_transfer.

import torch
import torch.nn as nn
import pennylane as qml

# ----- Baseline MLP -----
class MLP(nn.Module):
    def __init__(self, d):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d, 32), nn.ReLU(),
            nn.Linear(32, 16), nn.ReLU(),
            nn.Linear(16, 2)
        )
    def forward(self, x): 
        return self.net(x)

# ----- PQC Classifier (6-qubit re-uploading, 2 layers) -----
N_QUBITS = 6
N_LAYERS = 2
DEV = qml.device("default.qubit", wires=N_QUBITS, shots=None)

def angle_block(x):
    qml.templates.AngleEmbedding(x, wires=range(N_QUBITS), rotation="Y")

def var_block(weights):
    for w in weights:
        for i in range(N_QUBITS):
            qml.RY(w[i], wires=i)
        for i in range(N_QUBITS-1):
            qml.CZ(wires=[i, i+1])

@qml.qnode(DEV, interface="torch")
def qcircuit(x, weights):
    angle_block(x)
    var_block(weights[0])
    angle_block(x)      # simple re-upload
    var_block(weights[1])
    return [qml.expval(qml.PauliZ(i)) for i in range(N_QUBITS)]

class QClassifier(nn.Module):
    def __init__(self):
        super().__init__()
        self.weights = nn.Parameter(0.01*torch.randn(N_LAYERS, N_QUBITS))
        self.W = nn.Linear(N_QUBITS, 2)   # classical head
    def forward(self, x):
        # pad to N_QUBITS if needed
        if x.shape[1] < N_QUBITS:
            pad = torch.zeros(x.shape[0], N_QUBITS - x.shape[1], dtype=x.dtype, device=x.device)
            x = torch.cat([x, pad], dim=1)
        outs = []
        for i in range(x.shape[0]):
            out = qcircuit(x[i], self.weights)
            outs.append(torch.stack(out))
        Z = torch.stack(outs)               # (B, N_QUBITS)
        return self.W(Z)
