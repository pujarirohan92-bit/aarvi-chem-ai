from inorganic_templates import INORGANIC_TEMPLATES

def plan_inorganic(smiles: str, max_steps: int = 2):
    """
    Baseline planner for coordination compounds.
    Returns safe, explainable steps (no organic nonsense).
    """

    steps = []
    step_no = 1

    for tpl in INORGANIC_TEMPLATES:
        if step_no > max_steps:
            break

        step = {
            "step": step_no,
            "reaction_class": tpl["reaction_class"],
            "description": tpl["description"],
            "conditions": tpl["conditions"],
            "notes": tpl["notes"],
            "confidence": 0.65  # conservative baseline
        }

        steps.append(step)
        step_no += 1

    return {
        "route_id": 1,
        "route_score": 0.6,
        "logic": "Coordination chemistry baseline route",
        "steps": steps,
        "overall_confidence": 0.6
    }
