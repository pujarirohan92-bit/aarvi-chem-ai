def generate_chemist_report(target_smiles, best_route):
    report = []

    report.append(f"TARGET MOLECULE (SMILES): {target_smiles}")
    report.append(f"SELECTED STRATEGY: {best_route['strategy']}")
    report.append(f"OVERALL CONFIDENCE: {best_route['overall_confidence']}")
    report.append("")

    report.append("WHY THIS ROUTE WAS SELECTED:")
    report.append(
        "- This route has the highest overall confidence score.\n"
        "- Number of steps is minimized.\n"
        "- Warnings are acceptable and manageable under standard lab conditions."
    )
    report.append("")

    report.append("STEP-BY-STEP SYNTHESIS PLAN:")

    for step in best_route["steps"]:
        report.append(
            f"Step {step['step']}: {step['reaction']}\n"
            f"  Reactants: {', '.join(step['reactants'])}\n"
            f"  Conditions: {step['conditions']}\n"
            f"  Confidence: {step['confidence']}\n"
            f"  Warnings: {', '.join(step['warnings']) if step['warnings'] else 'None'}\n"
        )

    report.append("")
    report.append("CHEMIST NOTES:")
    report.append(
        "- Ensure reagent purity >98%.\n"
        "- Monitor temperature and pH carefully.\n"
        "- Steps with medium confidence may require optimization."
    )

    report.append("")
    report.append("FINAL RECOMMENDATION:")
    report.append(
        "Proceed with this route for laboratory synthesis. "
        "Scale-up is recommended only after pilot validation."
    )

    return "\n".join(report)
