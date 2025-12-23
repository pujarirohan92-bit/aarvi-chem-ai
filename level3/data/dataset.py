import csv
import os, sys
from rdkit import Chem
import torch

# path fixes
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(ROOT, "graph"))

from mol_to_graph import mol_to_graph


def get_mapped_bonds(mol):
    bonds = set()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetAtomMapNum()
        a2 = bond.GetEndAtom().GetAtomMapNum()
        if a1 > 0 and a2 > 0:
            bonds.add(tuple(sorted((a1, a2))))
    return bonds


def label_reaction(reactants, product):
    r_mol = Chem.MolFromSmiles(reactants)
    p_mol = Chem.MolFromSmiles(product)

    r_bonds = get_mapped_bonds(r_mol)
    p_bonds = get_mapped_bonds(p_mol)

    formed = p_bonds - r_bonds
    return formed


class ReactionDataset(torch.utils.data.Dataset):
    def __init__(self, csv_path):
        self.rows = []
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                self.rows.append(r)

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, idx):
        row = self.rows[idx]
        reactants = row["reactants"]
        product = row["product"]

        # graph from product (common practice)
        g = mol_to_graph(
            Chem.MolToSmiles(Chem.MolFromSmiles(product), canonical=True)
        )

        formed = label_reaction(reactants, product)

        # edge labels (demo: mark first edge positive if any formed bond exists)
        num_edges = g["edge_index"].shape[1]
        y = torch.zeros((num_edges, 1))
        if len(formed) > 0 and num_edges > 0:
            y[0] = 1.0

        return g, y
