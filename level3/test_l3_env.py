from rdkit import Chem
import torch

print("=== LEVEL-3 ENVIRONMENT TEST ===")

mol = Chem.MolFromSmiles("CC(=O)NC(C)C")
if mol is None:
    raise ValueError("RDKit failed to read SMILES")

print("RDKit OK")

x = torch.tensor([1.0, 2.0, 3.0])
print("Torch OK | Tensor sum =", x.sum().item())

atom_features = []
for atom in mol.GetAtoms():
    atom_features.append({
        "atomic_number": atom.GetAtomicNum(),
        "degree": atom.GetDegree(),
        "is_aromatic": atom.GetIsAromatic(),
        "formal_charge": atom.GetFormalCharge()
    })

print("Atom features sample:")
for a in atom_features[:3]:
    print(a)

print("LEVEL-3 ENV READY ✅")
