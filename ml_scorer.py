def score_step(step: dict) -> float:
    """
    Returns a confidence score (0–1) for a single reaction step
    """

    score = 1.0

    # Penalize harsh conditions
    if "heat" in step["conditions"].lower():
        score -= 0.1
    if "strong acid" in step["conditions"].lower():
        score -= 0.2

    # Penalize warnings
    if step["warning"]:
        score -= 0.15

    # Normalize
    return max(round(score, 2), 0.1)


def score_route(steps: list) -> float:
    """
    Score full synthetic route
    """

    if not steps:
        return 0.0

    total = sum(score_step(step) for step in steps)
    avg_score = total / len(steps)

    return round(avg_score, 2)
