import torch
import os, sys
from rdkit import Chem

ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "graph"))
sys.path.append(os.path.join(ROOT, "model"))

from mol_to_graph import mol_to_graph
from gnn_model import SimpleGNN

# -------- Load model --------
node_dim = 5
edge_dim = 6

model = SimpleGNN(node_dim, edge_dim)
model.load_state_dict(
    torch.load(os.path.join(ROOT, "checkpoints", "rc_gnn.pt"),
               map_location="cpu")
)
model.eval()

# -------- Predict reactive bonds --------
def predict_reactive_bonds(smiles, top_k=3):
    g = mol_to_graph(smiles)

    with torch.no_grad():
        scores = model(
            g["x"],
            g["edge_index"],
            g["edge_attr"]
        ).squeeze(1)

    # top-k edges
    topk = torch.topk(scores, min(top_k, scores.numel()))

    bonds = []
    edge_index = g["edge_index"]

    for idx in topk.indices.tolist():
        a1 = edge_index[0, idx].item()
        a2 = edge_index[1, idx].item()
        bonds.append((a1, a2, scores[idx].item()))

    return bonds


# -------- Test --------
if __name__ == "__main__":
    smi = "CC(=O)NC(C)C"
    preds = predict_reactive_bonds(smi)

    print("Predicted reactive bonds:")
    for b in preds:
        print(f"Atoms {b[0]}–{b[1]} | score={b[2]:.3f}")
