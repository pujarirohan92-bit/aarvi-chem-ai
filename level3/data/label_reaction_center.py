from rdkit import Chem

# -------- bond set extractor --------
def get_bond_set(mol):
    bonds = set()
    for bond in mol.GetBonds():
        a = bond.GetBeginAtom().GetAtomicNum()
        b = bond.GetEndAtom().GetAtomicNum()
        bonds.add(tuple(sorted((a, b))))
    return bonds


# -------- reaction center labeling --------
def label_reaction_center(reactant_smiles, product_smiles):
    r_mol = Chem.MolFromSmiles(reactant_smiles)
    p_mol = Chem.MolFromSmiles(product_smiles)

    if r_mol is None or p_mol is None:
        raise ValueError("Invalid SMILES")

    r_bonds = get_bond_set(r_mol)
    p_bonds = get_bond_set(p_mol)

    formed = p_bonds - r_bonds
    broken = r_bonds - p_bonds

    return {
        "formed_bonds": list(formed),
        "broken_bonds": list(broken)
    }


# -------- test --------
if __name__ == "__main__":
    reactants = "CC(=O)O.NC(C)C"
    product = "CC(=O)NC(C)C"

    labels = label_reaction_center(reactants, product)

    print("Formed bonds :", labels["formed_bonds"])
    print("Broken bonds :", labels["broken_bonds"])
