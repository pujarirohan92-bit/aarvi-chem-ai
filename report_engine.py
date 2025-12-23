def generate_synthesis_report(target_smiles: str, plan: dict) -> str:
    lines = []

    lines.append("AARVI CHEM AI – SYNTHESIS REPORT")
    lines.append("=" * 60)
    lines.append(f"Target Molecule: {target_smiles}")
    lines.append(f"Overall Confidence: {plan.get('overall_confidence')}")
    lines.append("=" * 60)
    lines.append("")

    for step in plan.get("steps", []):
        lines.append(f"STEP {step.get('step')}")
        lines.append("-" * 40)
        lines.append(f"Reaction: {step.get('reaction')}")
        lines.append(f"Bond Disconnection: {step.get('bond_disconnection')}")
        lines.append(
            "Reactants: " + ", ".join(step.get("reactants", []))
        )

        cond = step.get("conditions", {})
        lines.append(
            f"Conditions: Solvent={cond.get('solvent')}, "
            f"Temp={cond.get('temperature')}, "
            f"Catalyst={cond.get('catalyst')}"
        )

        lines.append(f"Step Confidence: {step.get('confidence')}")

        warnings = step.get("warnings", [])
        if warnings:
            lines.append("Warnings:")
            for w in warnings:
                lines.append(f" - {w}")
        else:
            lines.append("Warnings: None")

        lines.append("")

    return "<pre>" + "\n".join(lines) + "</pre>"
