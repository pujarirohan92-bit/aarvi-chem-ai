from rdkit import Chem


def label_reaction(reactants, products):
    """
    Very first baseline reaction class labeler.
    Rule-based, safe, interpretable.
    """

    rxn_smiles = ".".join(reactants) + ">>" + ".".join(products)

    # Simple heuristics (baseline)
    if "C(=O)N" in rxn_smiles or "NC(=O)" in rxn_smiles:
        return "amide formation"

    if "Br" in rxn_smiles or "Cl" in rxn_smiles:
        return "sn2 substitution"

    if "B(" in rxn_smiles and "Pd" in rxn_smiles:
        return "suzuki coupling"

    return "other"
