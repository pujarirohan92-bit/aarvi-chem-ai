def get_reagents_for_reaction(disconnection: str):
    """
    Returns reagents + conditions based on reaction type
    """

    # ---------- Amide coupling ----------
    if disconnection == "amide bond":
        return {
            "reaction_type": "Amide coupling",
            "reagents": ["EDC-HCl", "HOBt"],
            "base": "DIPEA",
            "solvent": "DMF",
            "temperature": "0–25 °C",
            "time": "4–6 h",
            "notes": "Standard peptide-style amide formation"
        }

    # ---------- Ester formation ----------
    if disconnection == "ester bond":
        return {
            "reaction_type": "Esterification",
            "reagents": ["SOCl2"],
            "solvent": "DCM",
            "temperature": "0 °C → RT",
            "time": "2–3 h",
            "notes": "Acid chloride mediated esterification"
        }

    # ---------- Generic fallback ----------
    return {
        "reaction_type": "Unknown / generic",
        "reagents": ["To be determined"],
        "solvent": "Unknown",
        "temperature": "Unknown",
        "time": "Unknown",
        "notes": "No standard conditions available"
    }
