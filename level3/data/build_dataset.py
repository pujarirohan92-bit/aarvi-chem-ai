import sys
import os
from rdkit import Chem

# -------- path fix (IMPORTANT) --------
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
GRAPH_DIR = os.path.join(ROOT_DIR, "graph")
sys.path.append(GRAPH_DIR)

from mol_to_graph import mol_to_graph


def get_mapped_bonds(mol):
    bonds = set()
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetAtomMapNum()
        a2 = bond.GetEndAtom().GetAtomMapNum()
        if a1 > 0 and a2 > 0:
            bonds.add(tuple(sorted((a1, a2))))
    return bonds


def label_reaction(mapped_reactants, mapped_product):
    r_mol = Chem.MolFromSmiles(mapped_reactants)
    p_mol = Chem.MolFromSmiles(mapped_product)

    r_bonds = get_mapped_bonds(r_mol)
    p_bonds = get_mapped_bonds(p_mol)

    formed = p_bonds - r_bonds
    broken = r_bonds - p_bonds

    return formed, broken


def build_sample():
    reactants = "[CH3:1][C:2](=[O:3])[O:4].[NH:5][CH:6]([CH3:7])[CH3:8]"
    product   = "[CH3:1][C:2](=[O:3])[NH:5][CH:6]([CH3:7])[CH3:8]"

    formed, broken = label_reaction(reactants, product)

    graph = mol_to_graph("CC(=O)NC(C)C")

    return {
        "formed_bonds": list(formed),
        "broken_bonds": list(broken),
        "graph": graph
    }


if __name__ == "__main__":
    sample = build_sample()
    print("Formed bonds :", sample["formed_bonds"])
    print("Broken bonds :", sample["broken_bonds"])
    print("Graph keys   :", sample["graph"].keys())
