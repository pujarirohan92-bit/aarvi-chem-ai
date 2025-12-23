from rdkit import Chem

def retro_amide(smiles):
    routes = []

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return routes

    # simple amide check
    if "C(=O)N" not in smiles:
        return routes

    # generate acid
    acid = smiles.replace("NC", "O").replace("N(", "O(")

    # generate amine
    parts = smiles.split("N")
    if len(parts) < 2:
        return routes

    amine = "N" + parts[1]

    routes.append({
        "reaction": "Amide formation (retro)",
        "acid": acid,
        "amine": amine,
        "conditions": "Acid chloride or EDC/HOBt"
    })

    return routes
