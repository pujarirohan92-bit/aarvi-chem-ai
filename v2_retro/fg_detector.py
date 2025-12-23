from rdkit import Chem

FG_SMARTS = {
    "acid": "C(=O)[OH]",
    "ester": "C(=O)O",
    "amide": "C(=O)N",
    "amine": "[NX3;H2,H1,H0]",
    "alcohol": "[OX2H]",
    "aldehyde": "[CH]=O",
    "aromatic": "a",
    "halide": "[Cl,Br,I,F]"
}

def detect_fg(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    found = []
    for name, smarts in FG_SMARTS.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            found.append(name)

    return list(set(found))
