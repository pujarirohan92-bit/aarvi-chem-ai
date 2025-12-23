from retro_engine import retro_amide

def retro_tree(smiles, depth=2):
    tree = {}

    # level 1
    routes = retro_amide(smiles)
    if not routes:
        return {}

    tree["product"] = smiles
    tree["routes"] = []

    for r in routes:
        node = {
            "reaction": r["reaction"],
            "precursors": {
                "acid": r["acid"],
                "amine": r["amine"]
            },
            "next_steps": {}
        }

        # level 2 – acid breakdown
        acid = r["acid"]
        if "C(=O)O" in acid:
            node["next_steps"]["acid"] = [
                f"{acid.replace('O', 'Cl')}  (acid chloride)"
            ]

        # level 2 – amine breakdown
        amine = r["amine"]
        if amine.startswith("N"):
            node["next_steps"]["amine"] = [
                "Alkyl halide + NH3"
            ]

        tree["routes"].append(node)

    return tree
