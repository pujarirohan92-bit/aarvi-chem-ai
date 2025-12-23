from rdkit import Chem
import torch

# -------- Atom feature encoder --------
def atom_features(atom):
    return [
        atom.GetAtomicNum(),          # atomic number
        atom.GetDegree(),             # number of neighbors
        atom.GetTotalNumHs(),         # attached hydrogens
        int(atom.GetIsAromatic()),    # aromatic or not
        atom.GetFormalCharge()        # formal charge
    ]

# -------- Bond feature encoder --------
def bond_features(bond):
    bt = bond.GetBondType()
    return [
        int(bt == Chem.BondType.SINGLE),
        int(bt == Chem.BondType.DOUBLE),
        int(bt == Chem.BondType.TRIPLE),
        int(bt == Chem.BondType.AROMATIC),
        int(bond.GetIsConjugated()),
        int(bond.IsInRing())
    ]

# -------- Molecule → Graph --------
def mol_to_graph(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    # Node features
    x = []
    for atom in mol.GetAtoms():
        x.append(atom_features(atom))
    x = torch.tensor(x, dtype=torch.float)

    # Edge index + edge attributes
    edge_index = []
    edge_attr = []

    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()

        bf = bond_features(bond)

        # undirected graph → add both directions
        edge_index.append([i, j])
        edge_attr.append(bf)

        edge_index.append([j, i])
        edge_attr.append(bf)

    edge_index = torch.tensor(edge_index, dtype=torch.long).t()
    edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    return {
        "x": x,
        "edge_index": edge_index,
        "edge_attr": edge_attr
    }

# -------- Simple test --------
if __name__ == "__main__":
    g = mol_to_graph("CC(=O)NC(C)C")

    print("Node feature shape :", g["x"].shape)
    print("Edge index shape  :", g["edge_index"].shape)
    print("Edge attr shape   :", g["edge_attr"].shape)
