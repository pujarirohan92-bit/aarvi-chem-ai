from condition_priors import CONDITION_PRIORS


def predict_conditions(reaction_class: str):
    prior = CONDITION_PRIORS.get(reaction_class)

    if not prior:
        return {
            "solvent": "Unknown",
            "base": "Unknown",
            "temperature": "Unknown"
        }

    return {
        "solvent": prior["solvent"][0],
        "base": prior["base"][0],
        "temperature": prior["temperature"]
    }
