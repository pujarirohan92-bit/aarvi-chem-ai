from model import detect_groups, classify_reaction

def predict_reaction(smiles1, smiles2):
    g1 = detect_groups(smiles1)
    g2 = detect_groups(smiles2)

    reaction = classify_reaction(g1, g2)

    return {
        "Reactant 1 groups": g1,
        "Reactant 2 groups": g2,
        "Predicted Reaction": reaction
    }


# ----- test -----
if __name__ == "__main__":
    r = predict_reaction("CC(=O)O", "N")  # acetic acid + ammonia
    print(r)
