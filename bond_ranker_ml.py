"""
STEP-24
Trained ML bond ranking + explainability
"""

import torch
import torch.nn as nn


class BondRankModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(3, 8),
            nn.ReLU(),
            nn.Linear(8, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.net(x)


def load_model(path="bond_ranker.pt"):
    model = BondRankModel()
    model.load_state_dict(torch.load(path))
    model.eval()
    return model


_model = load_model()


def featurize_bond(bond: dict):
    return torch.tensor([
        bond.get("priority", 0.5),
        len(bond.get("type", "")),
        len(bond.get("reason", ""))
    ], dtype=torch.float32)


def explain_bond(bond: dict):
    reasons = []
    if bond.get("priority", 0) > 0.8:
        reasons.append("High rule-based priority")
    if "amide" in bond.get("type", ""):
        reasons.append("Amide bond is retrosynthetically strategic")
    if len(bond.get("reason", "")) > 30:
        reasons.append("Strong chemical rationale")

    return "; ".join(reasons)


def rank_bonds(bonds: list):
    ranked = []

    with torch.no_grad():
        for b in bonds:
            x = featurize_bond(b)
            ml_score = float(_model(x).item())

            final_score = round(
                0.6 * b.get("priority", 0.5) + 0.4 * ml_score, 2
            )

            ranked.append({
                **b,
                "ml_score": round(ml_score, 2),
                "final_score": final_score,
                "explanation": explain_bond(b)
            })

    ranked.sort(key=lambda x: x["final_score"], reverse=True)
    return ranked
