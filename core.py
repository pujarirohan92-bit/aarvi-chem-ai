# core.py
from rdkit import Chem

def classify_reaction(smiles1, smiles2):
    """
    Input: two reactant SMILES
    Output: list of possible reaction types
    """

    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        return ["Invalid SMILES"]

    reactions = []

    acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    base = Chem.MolFromSmarts("[OH-]")
    amine = Chem.MolFromSmarts("[NX3;H2,H1,H0]")
    acyl = Chem.MolFromSmarts("[CX3](=O)[OX1]")

    # Acid–Base
    if (mol1.HasSubstructMatch(acid) and mol2.HasSubstructMatch(base)) or \
       (mol2.HasSubstructMatch(acid) and mol1.HasSubstructMatch(base)):
        reactions.append("Acid–Base Neutralization")

    # Amide formation
    if (mol1.HasSubstructMatch(amine) and mol2.HasSubstructMatch(acyl)) or \
       (mol2.HasSubstructMatch(amine) and mol1.HasSubstructMatch(acyl)):
        reactions.append("Amide Formation")

    # Esterification
    alcohol = Chem.MolFromSmarts("[OX2H]")
    if (mol1.HasSubstructMatch(acid) and mol2.HasSubstructMatch(alcohol)) or \
       (mol2.HasSubstructMatch(acid) and mol1.HasSubstructMatch(alcohol)):
        reactions.append("Esterification")

    if not reactions:
        reactions.append("No known reaction (rule-based)")

    return reactions
