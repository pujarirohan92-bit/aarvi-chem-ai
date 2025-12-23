from chemistry_type_detector import detect_chemistry_type

def route_by_chemistry(smiles: str) -> dict:
    info = detect_chemistry_type(smiles)
    ctype = info["type"]

    if ctype == "organic_small":
        return {"engine": "organic_engine", "info": info}

    if ctype == "organic_large":
        return {"engine": "organic_large_engine", "info": info}

    if ctype == "inorganic_coordination":
        return {"engine": "inorganic_engine", "info": info}

    if ctype == "organometallic":
        return {"engine": "organometallic_engine", "info": info}

    return {"engine": "unsupported", "info": info}
