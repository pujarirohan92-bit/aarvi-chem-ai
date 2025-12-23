# ml_engine.py — STEP-47
MODEL_FEAS = None
MODEL_YIELD = None


def load_models(feas_model, yield_model):
    global MODEL_FEAS, MODEL_YIELD
    MODEL_FEAS = feas_model
    MODEL_YIELD = yield_model


def reaction_feasibility(step: dict) -> float:
    if MODEL_FEAS is None:
        return 0.5

    class_map = {
        "amide formation": 0,
        "sn2 substitution": 1,
        "suzuki coupling": 2,
        "other": 3
    }

    class_id = class_map.get(step.get("reaction", "other"), 3)
    n_reactants = len(step.get("reactants", []))
    n_products = 1

    X = [[class_id, n_reactants, n_products]]
    return round(float(MODEL_FEAS.predict_proba(X)[0][1]), 2)


def reaction_yield(step: dict) -> float:
    if MODEL_YIELD is None:
        return 0.5

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

    class_id = class_map.get(step.get("reaction", "other"), 3)
    n_reactants = len(step.get("reactants", []))
    n_products = 1

    temp_bucket = temp_map.get(step.get("conditions", {}).get("temperature"), 1)

    X = [[class_id, n_reactants, n_products, temp_bucket]]
    return round(float(MODEL_YIELD.predict(X)[0]), 2)
