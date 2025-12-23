# chem_synonyms.py

SYNONYMS = {
    "acid derivative": ["acid chloride", "acyl chloride", "ester", "anhydride"],
    "amine": ["amine", "primary amine", "secondary amine"],
    "alcohol": ["alcohol", "oh", "hydroxyl"],
    "carboxylic acid": ["acid", "cooh", "carboxylic acid"]
}

def normalize_groups(user_groups):
    normalized = set()

    for g in user_groups:
        for key, values in SYNONYMS.items():
            if g == key or g in values:
                normalized.add(key)

    return list(normalized)
