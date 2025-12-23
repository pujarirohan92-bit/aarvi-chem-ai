from fg_detector import detect_fg
from reactions_db import REACTIONS

def predict_reactions(smiles):
    fgs = detect_fg(smiles)
    matches = []

    for r in REACTIONS:
        if all(n in fgs for n in r["needs"]):
            matches.append(r)

    return matches, fgs
