INORGANIC_TEMPLATES = [
    {
        "reaction_class": "ligand_exchange",
        "description": "Halide to ammine substitution",
        "typical_metals": ["Pt", "Pd"],
        "example": "K2PtCl4 → cis-[Pt(NH3)2Cl2]",
        "conditions": {
            "solvent": "Water",
            "reagent": "NH3 (aq)",
            "temperature": "25–60 °C"
        },
        "notes": [
            "Stepwise substitution",
            "cis product favored under controlled conditions"
        ]
    },
    {
        "reaction_class": "aquation",
        "description": "Halide replacement by water",
        "typical_metals": ["Pt"],
        "example": "[PtCl2(NH3)2] → [Pt(H2O)Cl(NH3)2]+",
        "conditions": {
            "solvent": "Water",
            "reagent": "AgNO3 (optional)",
            "temperature": "Room temperature"
        },
        "notes": [
            "Used to activate platinum complexes",
            "Precedes ligand substitution"
        ]
    }
]
