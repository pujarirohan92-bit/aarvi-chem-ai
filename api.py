from fastapi import FastAPI
from pydantic import BaseModel
from typing import List, Dict

from report_engine import generate_synthesis_report

app = FastAPI(
    title="Aarvi Chem AI",
    description="Organic & Inorganic Multi-Step Retrosynthesis Intelligence",
    version="2.0"
)

# ---------- INPUT SCHEMA ----------

class SynthesisRequest(BaseModel):
    smiles: str

# ---------- ROOT ----------

@app.get("/")
def home():
    return {
        "message": "Aarvi Chem AI API is running",
        "endpoints": {
            "docs": "/docs",
            "generate": "/generate"
        }
    }

# ---------- CORE LOGIC (TEMP MOCK PLANNER) ----------

def mock_retrosynthesis_engine(smiles: str) -> List[Dict]:
    """
    Temporary planner.
    Later this will be replaced by ML / GNN / template engine.
    """

    return [
        {
            "reaction": "Amide / Sulfonylurea Formation",
            "reactants": ["Amine fragment", "Sulfonyl isocyanate"],
            "products": ["Target molecule"],
            "conditions": "Base, aprotic solvent, 0–25°C",
            "notes": "Heuristic AI suggestion"
        }
    ]

# ---------- API ENDPOINT ----------

@app.post("/generate")
def generate_synthesis(req: SynthesisRequest):
    smiles = req.smiles

    plan = mock_retrosynthesis_engine(smiles)

    report_text = generate_synthesis_report(smiles, plan)

    return {
        "status": "success",
        "input_smiles": smiles,
        "report": report_text,
        "plan": plan
    }
