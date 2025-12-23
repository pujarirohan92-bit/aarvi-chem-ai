def rank_routes(routes):
    ranked = []

    for r in routes:
        total_warnings = sum(len(s["warnings"]) for s in r["steps"])
        n_steps = len(r["steps"])

        score = (
            r["overall_confidence"]
            - 0.1 * total_warnings
            - 0.05 * n_steps
        )

        ranked.append({
            "route_id": r["route_id"],
            "strategy": r["strategy"],
            "score": round(score, 3),
            "overall_confidence": r["overall_confidence"],
            "steps": r["steps"]
        })

    ranked.sort(key=lambda x: x["score"], reverse=True)
    return ranked
