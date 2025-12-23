from rdkit import Chem

# Functional group SMARTS
PATTERNS = {
    "amine": "[NX3;H2,H1;!$(NC=O)]",
    "phenol": "c[OX2H]",
    "alcohol": "[OX2H][CX4]",
    "acid_anhydride": "C(=O)OC(=O)",
    "amide": "NC=O"
}

def detect_functional_groups(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    found = []
    for name, smarts in PATTERNS.items():
        patt = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(patt):
            found.append(name)

    return found


# -------- TEST CASES ----------
if __name__ == "__main__":

    print("4-aminophenol groups:")
    print(detect_functional_groups("Nc1ccc(O)cc1"))

    print("\nAcetic anhydride groups:")
    print(detect_functional_groups("CC(=O)OC(=O)C"))