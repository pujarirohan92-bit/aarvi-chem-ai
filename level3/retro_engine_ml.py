import os, sys
from rdkit import Chem

# -------- PATHS --------
ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "graph"))
sys.path.append(os.path.join(ROOT, "model"))
sys.path.append(os.path.join(ROOT, "..", "v2_retro"))

from mol_to_graph import mol_to_graph
from gnn_model import SimpleGNN
from engine import predict_reactions   # LEVEL-2 rules

import torch

# -------- Load trained ML model --------
node_dim = 5
edge_dim = 6

model = SimpleGNN(node_dim, edge_dim)
model.load_state_dict(
    torch.load(os.path.join(ROOT, "checkpoints", "rc_gnn.pt"),
               map_location="cpu")
)
model.eval()


# -------- ML: predict reactive bonds --------
def ml_predict_bonds(smiles, top_k=2):
    g = mol_to_graph(smiles)

    with torch.no_grad():
        scores = model(
            g["x"],
            g["edge_index"],
            g["edge_attr"]
        ).squeeze(1)

    topk = torch.topk(scores, min(top_k, scores.numel()))

    bonds = []
    for idx in topk.indices.tolist():
        a1 = g["edge_index"][0, idx].item()
        a2 = g["edge_index"][1, idx].item()
        bonds.append({
            "atoms": (a1, a2),
            "score": scores[idx].item()
        })

    return bonds


# -------- Apply bond disconnection --------
def disconnect_bond(smiles, atom_pair):
    mol = Chem.MolFromSmiles(smiles)
    rw = Chem.RWMol(mol)

    rw.RemoveBond(atom_pair[0], atom_pair[1])
    frags = Chem.GetMolFrags(rw.GetMol(), asMols=True)

    return [Chem.MolToSmiles(f) for f in frags]


# -------- MAIN retrosynthesis --------
def retrosynthesis_ml(smiles):
    routes = []

    # 1️⃣ ML predicts reactive bonds
    ml_bonds = ml_predict_bonds(smiles)

    # 2️⃣ Level-2 rule suggestions
    lvl2_rxns, fgs = predict_reactions(smiles)

    for b in ml_bonds:
        precursors = disconnect_bond(smiles, b["atoms"])

        route = {
            "ml_score": round(b["score"], 3),
            "broken_bond": b["atoms"],
            "precursors": precursors,
            "level2_reactions": [r["name"] for r in lvl2_rxns]
        }

        routes.append(route)

    # 3️⃣ Rank routes (ML score + chemistry)
    routes = sorted(routes, key=lambda x: x["ml_score"], reverse=True)

    return routes


# -------- TEST --------
if __name__ == "__main__":
    smi = "CC(=O)NC(C)C"
    routes = retrosynthesis_ml(smi)

    print("RETROSYNTHESIS RESULTS\n")
    for i, r in enumerate(routes, 1):
        print(f"ROUTE {i}")
        print(" ML score :", r["ml_score"])
        print(" Broken bond :", r["broken_bond"])
        print(" Precursors :", r["precursors"])
        print(" Possible reactions :", r["level2_reactions"])
        print("-" * 40)
