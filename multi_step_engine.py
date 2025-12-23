def build_multistep_chain(smiles: str, max_steps: int):
    steps = []

    for i in range(1, max_steps + 1):
        step = {
            "step": i,
            "reaction": f"Generic reaction step {i}",
            "reactants": [f"Reactant_{i}A", f"Reactant_{i}B"],
            "products": [f"Intermediate_{i}"],
            "conditions": f"Solvent, heat, catalyst (step {i})",
            "warning": "Check temperature control"
        }
        steps.append(step)

    return steps
