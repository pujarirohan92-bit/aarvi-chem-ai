def build_yield_features(records):
    X = []
    y = []

    class_map = {
        "amide formation": 0,
        "sn2 substitution": 1,
        "suzuki coupling": 2,
        "other": 3
    }

    temp_map = {
        "0–25 °C": 0,
        "40–80 °C": 1,
        "70–100 °C": 2
    }

    for r in records:
        rc = r.get("reaction_class", "other")
        class_id = class_map.get(rc, 3)

        n_reactants = len(r.get("reactants", []))
        n_products = len(r.get("products", []))

        cond = r.get("conditions", {})
        temp_bucket = temp_map.get(cond.get("temperature"), 1)

        X.append([class_id, n_reactants, n_products, temp_bucket])

        # Yield proxy: map bins to continuous
        outcome = r.get("outcome", {})
        ybin = outcome.get("yield_bin", "medium")
        yval = {"low": 0.3, "medium": 0.6, "high": 0.85}.get(ybin, 0.5)
        y.append(yval)

    return X, y
