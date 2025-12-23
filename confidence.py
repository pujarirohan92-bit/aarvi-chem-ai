def calculate_confidence(rule_strength, fg_score, commonality, penalty):
    score = (
        0.45 * rule_strength +
        0.30 * fg_score +
        0.15 * commonality -
        0.10 * penalty
    )
    return round(max(0, min(score, 1)), 2)


def explain_confidence(rule_strength, fg_score, commonality, penalty):
    reasons = []

    reasons.append(
        "Strong known reaction rule"
        if rule_strength > 0.7 else
        "Weak or uncommon rule"
    )

    reasons.append(
        "Functional groups compatible"
        if fg_score > 0.7 else
        "Possible FG interference"
    )

    reasons.append(
        "Common laboratory reaction"
        if commonality > 0.6 else
        "Rare or specialized reaction"
    )

    if penalty > 0.4:
        reasons.append("Sensitive / incompatible groups present")

    return reasons
