from fastapi import FastAPI
from pydantic import BaseModel
from fastapi.responses import FileResponse
from typing import List, Dict

from report_engine import generate_synthesis_report
from pdf_generator import generate_pdf

app = FastAPI(
    title="Aarvi Chem AI",
    version="2.1"
)

class SynthesisRequest(BaseModel):
    smiles: str

@app.get("/")
def home():
    return {"message": "Aarvi Chem AI running"}

def mock_retrosynthesis_engine(smiles: str) -> List[Dict]:
    return [
        {
            "reaction": "Sulfonylurea / Amide Formation",
            "reactants": ["Amine fragment", "Sulfonyl isocyanate"],
            "products": ["Target molecule"],
            "conditions": "Base, aprotic solvent, 0–25°C",
            "notes": "AI heuristic route (non-patent)"
        }
    ]

@app.post("/generate-pdf")
def generate_pdf_report(req: SynthesisRequest):
    smiles = req.smiles
    plan = mock_retrosynthesis_engine(smiles)

    report_text = generate_synthesis_report(smiles, plan)

    pdf_path, filename = generate_pdf(report_text)

    return FileResponse(
        pdf_path,
        media_type="application/pdf",
        filename=filename
    )
