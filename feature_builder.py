def build_features(records):
    """
    Convert reaction records into ML features.
    """
    X = []
    y = []

    for r in records:
        reaction_class = r.get("reaction_class", "other")

        # Simple encoding
        class_map = {
            "amide formation": 0,
            "sn2 substitution": 1,
            "suzuki coupling": 2,
            "other": 3
        }

        class_id = class_map.get(reaction_class, 3)

        n_reactants = len(r.get("reactants", []))
        n_products = len(r.get("products", []))

        X.append([class_id, n_reactants, n_products])

        outcome = r.get("outcome", {})
        y.append(outcome.get("success", 0))

    return X, y
