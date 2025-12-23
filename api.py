from fastapi import FastAPI
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import List, Dict

from pdf_generator import create_pdf

app = FastAPI(
    title="Aarvi Chem AI",
    version="2.2"
)

class SynthesisRequest(BaseModel):
    smiles: str

@app.get("/")
def root():
    return {"status": "Aarvi Chem AI running"}

def mock_ai_engine(smiles: str) -> Dict:
    report = f"""
AARVI CHEM AI – SYNTHESIS REPORT
================================

Target Molecule (SMILES):
{smiles}

Step 1:
Reaction Type: Amide / Sulfonylurea Formation
Reactants:
- Amine fragment
- Sulfonyl isocyanate

Products:
- Target molecule

Conditions:
- Base
- Aprotic solvent
- 0–25 °C

Notes:
AI heuristic retrosynthesis suggestion.

DISCLAIMER:
This synthesis route is AI-generated and may differ from
industrial or patented routes. Experimental validation required.
"""
    return {"report": report}

@app.post("/generate-report")
def generate_report(req: SynthesisRequest):
    data = mock_ai_engine(req.smiles)
    return {
        "status": "success",
        "smiles": req.smiles,
        "report": data["report"]
    }

@app.post("/download-pdf")
def download_pdf(req: SynthesisRequest):
    data = mock_ai_engine(req.smiles)
    pdf_path, filename = create_pdf(data["report"])

    return FileResponse(
        pdf_path,
        media_type="application/pdf",
        filename=filename
    )
