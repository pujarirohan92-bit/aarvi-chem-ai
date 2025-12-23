def analyze_risk(route: dict):
    """
    Analyze chemical risks and return warnings + penalty
    """

    warnings = []
    penalty = 0.0

    disconnection = route.get("disconnection", "")
    conditions = route.get("conditions", {})

    reagents = conditions.get("reagents", [])
    temperature = conditions.get("temperature", "")

    # -------- Rule 1: Harsh reagents --------
    if "SOCl2" in reagents:
        warnings.append("SOCl2 is moisture sensitive and harsh")
        penalty += 0.15

    # -------- Rule 2: Amide coupling risks --------
    if disconnection == "amide bond":
        warnings.append("Check for steric hindrance near amine")
        penalty += 0.05

    # -------- Rule 3: High temperature --------
    if ">" in temperature or "reflux" in temperature.lower():
        warnings.append("High temperature may cause side reactions")
        penalty += 0.10

    # -------- Rule 4: Low base confidence --------
    base_conf = route.get("confidence", 0.5)
    if base_conf < 0.5:
        warnings.append("Low confidence route – not recommended")
        penalty += 0.20

    return {
        "warnings": warnings,
        "risk_penalty": round(min(penalty, 0.5), 2)
    }
