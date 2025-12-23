"""
STEP-22
Bond disconnection intelligence (rule-based)
"""

from rdkit import Chem


def find_disconnections(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    suggestions = []

    if mol is None:
        return suggestions

    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        # Rule 1: Amide bond (C–N next to C=O)
        if (
            a1.GetSymbol() == "C"
            and a2.GetSymbol() == "N"
            and any(nb.GetSymbol() == "O" for nb in a1.GetNeighbors())
        ):
            suggestions.append({
                "bond": f"{a1.GetIdx()}-{a2.GetIdx()}",
                "type": "amide_disconnection",
                "reason": "Amide bond – common retrosynthetic cut",
                "priority": 0.9
            })

        # Rule 2: Ester bond (C–O next to C=O)
        if (
            a1.GetSymbol() == "C"
            and a2.GetSymbol() == "O"
            and any(nb.GetSymbol() == "O" for nb in a1.GetNeighbors())
        ):
            suggestions.append({
                "bond": f"{a1.GetIdx()}-{a2.GetIdx()}",
                "type": "ester_disconnection",
                "reason": "Ester bond – breaks into acid + alcohol",
                "priority": 0.8
            })

    # Sort by priority
    suggestions.sort(key=lambda x: x["priority"], reverse=True)
    return suggestions
