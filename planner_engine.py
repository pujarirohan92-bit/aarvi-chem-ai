from reaction_templates import REACTION_TEMPLATES
from ml_engine import reaction_feasibility, reaction_yield
from condition_engine import predict_conditions


def get_template_by_reaction(reaction_name: str):
    for t in REACTION_TEMPLATES:
        if t.get("reaction") == reaction_name:
            return t
    return None


def plan_multistep(smiles: str, max_steps: int):
    steps = []

    for i in range(1, max_steps + 1):
        reaction_name = "amide formation" if i == 1 else "sn2 substitution"
        template = get_template_by_reaction(reaction_name)
        conditions = predict_conditions(reaction_name)

        step = {
            "step": i,
            "reaction": reaction_name,
            "bond_disconnection": template.get("bond") if template else "Unknown",
            "reactants": template.get("reactants", []) if template else [],
            "conditions": conditions,
            "confidence": 0.0,
            "yield": 0.0,
            "warnings": []
        }

        feas = reaction_feasibility(step)
        yld = reaction_yield(step)

        step["confidence"] = feas
        step["yield"] = yld

        if feas < 0.4:
            step["warnings"].append("Low feasibility")
        if yld < 0.4:
            step["warnings"].append("Low expected yield")

        steps.append(step)

    overall_confidence = round(
        sum(s["confidence"] for s in steps) / len(steps), 2
    )

    return {
        "target_smiles": smiles,
        "steps": steps,
        "overall_confidence": overall_confidence
    }
