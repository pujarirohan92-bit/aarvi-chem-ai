REACTIONS = [

    # ===================== ACID DERIVATIVES =====================
    {
        "name": "Amide formation",
        "needs": ["acid", "amine"],
        "retro": "acid + amine",
        "conditions": "EDC/HOBt or acid chloride"
    },
    {
        "name": "Esterification",
        "needs": ["acid", "alcohol"],
        "retro": "acid + alcohol",
        "conditions": "H2SO4 reflux"
    },
    {
        "name": "Acid chloride formation",
        "needs": ["acid"],
        "retro": "acid chloride",
        "conditions": "SOCl2 / oxalyl chloride"
    },

    # ===================== N-CHEMISTRY =====================
    {
        "name": "Sulfonamide formation",
        "needs": ["amine"],
        "retro": "sulfonyl chloride + amine",
        "conditions": "Et3N / pyridine"
    },
    {
        "name": "Urea formation",
        "needs": ["amine"],
        "retro": "isocyanate + amine",
        "conditions": "Phosgene equivalent"
    },
    {
        "name": "Reductive amination",
        "needs": ["amine"],
        "retro": "carbonyl + amine",
        "conditions": "NaBH3CN"
    },

    # ===================== SUBSTITUTION =====================
    {
        "name": "SN2 substitution",
        "needs": ["halide"],
        "retro": "alkyl halide + nucleophile",
        "conditions": "Polar aprotic solvent"
    },
    {
        "name": "Williamson ether synthesis",
        "needs": ["alcohol", "halide"],
        "retro": "alkyl halide + alkoxide",
        "conditions": "NaH / DMF"
    },

    # ===================== CARBONYL REACTIONS =====================
    {
        "name": "Reduction (carbonyl)",
        "needs": ["aldehyde"],
        "retro": "aldehyde",
        "conditions": "NaBH4"
    },
    {
        "name": "Reduction (ester/acid)",
        "needs": ["ester"],
        "retro": "ester",
        "conditions": "LiAlH4"
    },

    # ===================== AROMATIC =====================
    {
        "name": "Friedel–Crafts acylation",
        "needs": ["aromatic"],
        "retro": "arene + acyl chloride",
        "conditions": "AlCl3"
    },
    {
        "name": "Friedel–Crafts alkylation",
        "needs": ["aromatic"],
        "retro": "arene + alkyl halide",
        "conditions": "AlCl3"
    },
    {
        "name": "Nitration",
        "needs": ["aromatic"],
        "retro": "arene",
        "conditions": "HNO3 / H2SO4"
    },
    {
        "name": "Sulfonation",
        "needs": ["aromatic"],
        "retro": "arene",
        "conditions": "SO3 / H2SO4"
    },

    # ===================== PROTECTION / DEPROTECTION =====================
    {
        "name": "Boc protection",
        "needs": ["amine"],
        "retro": "amine",
        "conditions": "Boc2O"
    },
    {
        "name": "Boc deprotection",
        "needs": ["amine"],
        "retro": "protected amine",
        "conditions": "TFA"
    }
]
