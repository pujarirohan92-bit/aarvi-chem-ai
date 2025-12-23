# reaction_templates.py

REACTION_TEMPLATES = [
    {
        "reaction": "amide formation",
        "bond": "C–N",
        "reactants": ["Acid chloride", "Amine"],
        "conditions": {
            "solvent": "DCM",
            "temperature": "0–5 °C",
            "base": "Et3N"
        }
    },
    {
        "reaction": "sn2 substitution",
        "bond": "C–X",
        "reactants": ["Alkyl halide", "Nucleophile"],
        "conditions": {
            "solvent": "DMF",
            "temperature": "60 °C"
        }
    }
]
