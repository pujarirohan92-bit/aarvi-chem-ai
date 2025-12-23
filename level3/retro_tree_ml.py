import os
import sys
import torch
from rdkit import Chem

# ================= PATH SETUP =================
ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "graph"))
sys.path.append(os.path.join(ROOT, "model"))
sys.path.append(os.path.join(ROOT, "..", "v2_retro"))

from mol_to_graph import mol_to_graph
from gnn_model import SimpleGNN
from engine import predict_reactions   # Level-2 rules

# ================= LOAD MODEL =================
node_dim = 5
edge_dim = 6

model = SimpleGNN(node_dim, edge_dim)
model.load_state_dict(
    torch.load(os.path.join(ROOT, "checkpoints", "rc_gnn.pt"),
               map_location="cpu")
)
model.eval()

# ================= ML BOND PREDICTION =================
def ml_predict_bonds(smiles, top_k=2):
    g = mol_to_graph(smiles)

    # 🔐 GUARD: no bonds → no ML
    if g["edge_attr"].shape[0] == 0:
        return []

    with torch.no_grad():
        scores = model(
            g["x"],
            g["edge_index"],
            g["edge_attr"]
        ).squeeze(1)

    k = min(top_k, scores.numel())
    topk = torch.topk(scores, k)

    bonds = []
    for idx in topk.indices.tolist():
        a1 = g["edge_index"][0, idx].item()
        a2 = g["edge_index"][1, idx].item()
        bonds.append({
            "atoms": (a1, a2),
            "score": float(scores[idx].item())
        })

    return bonds

# ================= DISCONNECT BOND =================
def disconnect_bond(smiles, atom_pair):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    rw = Chem.RWMol(mol)
    try:
        rw.RemoveBond(atom_pair[0], atom_pair[1])
    except Exception:
        return []

    frags = Chem.GetMolFrags(rw.GetMol(), asMols=True)
    return [Chem.MolToSmiles(f, canonical=True) for f in frags]

# ================= UTILITY =================
def canonical(smiles):
    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m, canonical=True) if m else smiles

# ================= TREE BUILDER =================
def build_tree(smiles, depth, top_k, visited):
    smiles = canonical(smiles)

    if depth == 0 or smiles in visited:
        return {
            "smiles": smiles,
            "stop": True
        }

    visited.add(smiles)

    node = {
        "smiles": smiles,
        "routes": []
    }

    # Level-2 chemistry context
    lvl2_rxns, _ = predict_reactions(smiles)

    # ML bond prediction
    ml_bonds = ml_predict_bonds(smiles, top_k=top_k)
    if not ml_bonds:
        return {
            "smiles": smiles,
            "stop": True
        }

    for b in ml_bonds:
        precursors = disconnect_bond(smiles, b["atoms"])
        if not precursors:
            continue

        children = []
        for p in precursors:
            child = build_tree(
                p,
                depth=depth - 1,
                top_k=top_k,
                visited=visited.copy()
            )
            children.append(child)

        route = {
            "ml_score": round(b["score"], 3),
            "broken_bond": b["atoms"],
            "precursors": precursors,
            "level2_reactions": [r["name"] for r in lvl2_rxns],
            "children": children
        }

        node["routes"].append(route)

    node["routes"] = sorted(
        node["routes"],
        key=lambda r: r["ml_score"],
        reverse=True
    )

    return node

# ================= PUBLIC API =================
def retrosynthesis_tree(smiles, depth=2, top_k=2):
    return build_tree(
        smiles,
        depth=depth,
        top_k=top_k,
        visited=set()
    )

# ================= TEST =================
if __name__ == "__main__":
    import json
    target = "CC(=O)NC(C)C"
    tree = retrosynthesis_tree(target, depth=2, top_k=2)
    print(json.dumps(tree, indent=2))
