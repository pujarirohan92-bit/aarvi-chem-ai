# engine.py

def predict_reaction(reactants):
    """
    reactants: string (user input)
    return: predicted product (string)
    """

    reactants = reactants.lower()

    if "acid" in reactants and "base" in reactants:
        return "Salt + Water (Acid-Base Neutralization)"

    elif "amine" in reactants and "acid chloride" in reactants:
        return "Amide formation"

    elif "alcohol" in reactants and "acid" in reactants:
        return "Ester formation"

    else:
        return "Reaction not found in database"
