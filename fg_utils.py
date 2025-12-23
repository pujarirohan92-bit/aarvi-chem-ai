from rdkit import Chem

FG_SMARTS = {
    "amine": "[NX3;H2,H1;!$(NC=O)]",
    "acid_chloride": "C(=O)Cl",
    "alcohol": "[OX2H]",
    "carboxylic_acid": "C(=O)[OX2H1]",
    "nitro": "[NX3](=O)=O"
}

def detect_functional_groups(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    found = []
    for fg, smarts in FG_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(patt):
            found.append(fg)
    return found
