"""
STEP-26.2
Bond disconnection → Reaction template mapper
"""

from reaction_templates import get_reaction_template


def map_bond_to_reaction(bond: dict):
    """
    Given a bond disconnection, return reaction details
    """
    if not bond:
        return None

    bond_type = bond.get("type")
    template = get_reaction_template(bond_type)

    if not template:
        return {
            "reaction": "unknown_reaction",
            "reagents": [],
            "conditions": "Not specified",
            "confidence": 0.3
        }

    return {
        "reaction": template["name"],
        "reagents": template["typical_reagents"],
        "conditions": template["conditions"],
        "confidence": template["confidence"]
    }
