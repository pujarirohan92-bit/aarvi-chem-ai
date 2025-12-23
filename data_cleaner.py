from rdkit import Chem


def clean_smiles(smiles: str):
    """
    Remove atom mapping, validate SMILES
    """
    if not smiles or not isinstance(smiles, str):
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Remove atom mapping numbers
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)

        return Chem.MolToSmiles(mol, canonical=True)

    except Exception:
        return None


def clean_reaction_record(record: dict):
    """
    Takes raw reaction record and returns cleaned record or None
    """
    reactants = []
    for smi in record.get("reactants", []):
        cs = clean_smiles(smi)
        if cs:
            reactants.append(cs)

    products = []
    for smi in record.get("products", []):
        cs = clean_smiles(smi)
        if cs:
            products.append(cs)

    if not reactants or not products:
        return None

    cleaned = {
        "reactants": reactants,
        "products": products,
        "reaction_class": record.get("reaction_class"),
        "conditions": record.get("conditions", {}),
        "outcome": record.get("outcome", {})
    }

    return cleaned
