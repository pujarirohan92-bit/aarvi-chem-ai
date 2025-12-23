from planner_engine import plan_multistep


def run_benchmark(benchmarks, max_steps=3):
    results = []

    for b in benchmarks:
        plan = plan_multistep(b["smiles"], max_steps)
        predicted_rxns = [s["reaction"] for s in plan["steps"]]

        match = any(
            rxn in predicted_rxns for rxn in b["expected_reactions"]
        )

        results.append({
            "name": b["name"],
            "expected": b["expected_reactions"],
            "predicted": predicted_rxns,
            "match": match,
            "overall_confidence": plan["overall_confidence"]
        })

    return results
