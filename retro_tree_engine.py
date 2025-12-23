from rdkit import Chem


def retrosynthesis_tree(smiles: str, depth: int = 2, top_k: int = 2):
    """
    Rule-based retrosynthesis tree with confidence scoring
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    routes = []

    # -------- Rule 1: Amide disconnection --------
    if "C(=O)N" in smiles:
        routes.append({
            "disconnection": "amide bond",
            "precursors": ["CC(=O)Cl", "CN"],
            "confidence": 0.82,
            "reason": "Standard amide coupling"
        })

    # -------- Rule 2: Ester-like fallback --------
    if "C(=O)O" in smiles:
        routes.append({
            "disconnection": "ester bond",
            "precursors": ["CC(=O)Cl", "CO"],
            "confidence": 0.65,
            "reason": "Common esterification"
        })

    # -------- Rule 3: Generic fallback --------
    routes.append({
        "disconnection": "generic cleavage",
        "precursors": ["fragment_A", "fragment_B"],
        "confidence": 0.35,
        "reason": "Low confidence generic split"
    })

    # -------- Sort by confidence --------
    routes = sorted(routes, key=lambda x: x["confidence"], reverse=True)

    # -------- Apply top_k --------
    return routes[:top_k]
