from rdkit import Chem

# ---------- Functional group detection ----------
def detect_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)
    groups = []

    patterns = {
        "acid": "C(=O)[OX2H1]",
        "base": "[OH-]",
        "amine": "[NX3;H2,H1]",
        "alcohol": "[OX2H]"
    }

    for name, smarts in patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            groups.append(name)

    return groups


# ---------- Reaction type ----------
def classify_reaction(groups1, groups2):
    if "acid" in groups1 and "base" in groups2:
        return "Acid-Base Reaction"
    if "acid" in groups1 and "amine" in groups2:
        return "Acid-Amine Salt Formation"
    if "acid" in groups1 and "alcohol" in groups2:
        return "Esterification"
    return "Unknown Reaction"
