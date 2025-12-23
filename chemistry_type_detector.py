import re

METALS = [
    "Pt", "Pd", "Ni", "Cu", "Ru", "Rh", "Ir", "Fe", "Co", "Zn",
    "Ag", "Au", "Sn", "Al", "Mg"
]

def detect_chemistry_type(smiles: str) -> dict:
    """
    Detect chemistry type from SMILES / formula-like strings.
    Returns a dict with type + notes.
    """

    if not smiles or not isinstance(smiles, str):
        return {"type": "unknown", "note": "Invalid input"}

    # Normalize
    s = smiles.strip()

    # 1) Inorganic / coordination: presence of metal tokens
    for m in METALS:
        # metal token as standalone or bracketed
        if re.search(rf"(\[{m}|\b{m}\b)", s):
            return {
                "type": "inorganic_coordination",
                "note": f"Metal center detected: {m}"
            }

    # 2) Organometallic: metal + carbon fragments (rough heuristic)
    if any(m in s for m in METALS) and ("C" in s or "c" in s):
        return {
            "type": "organometallic",
            "note": "Metal + carbon framework detected"
        }

    # 3) Organic: count heavy atoms (rough size proxy)
    heavy_atoms = len(re.findall(r"[A-Z][a-z]?", s))
    ring_count = s.count("1") + s.count("2") + s.count("3")

    if heavy_atoms >= 60 or ring_count >= 6:
        return {
            "type": "organic_large",
            "note": f"Large organic molecule (heavy_atoms={heavy_atoms})"
        }

    if "C" in s or "c" in s:
        return {
            "type": "organic_small",
            "note": "Carbon-based organic molecule"
        }

    # Fallback
    return {"type": "unknown", "note": "Could not classify"}
