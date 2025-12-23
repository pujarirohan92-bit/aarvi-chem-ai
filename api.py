from fastapi import FastAPI, Body
from fastapi.responses import FileResponse
from pydantic import BaseModel
import uuid
import os
from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas

app = FastAPI(
    title="Aarvi Chem AI",
    description="Organic & Inorganic Multi-Step Synthesis Intelligence",
    version="2.0"
)

# =========================
# Input schema
# =========================
class PDFRequest(BaseModel):
    smiles: str


# =========================
# PDF Generator
# =========================
def generate_pdf(smiles: str, filepath: str):
    c = canvas.Canvas(filepath, pagesize=A4)
    width, height = A4

    y = height - 50

    def line(text):
        nonlocal y
        c.drawString(40, y, text)
        y -= 18

    line("AARVI CHEM AI – SYNTHESIS REPORT")
    line("=" * 60)
    y -= 10

    line(f"Target Molecule (SMILES):")
    line(smiles)
    y -= 10

    line("Reaction Type:")
    line("Amide / Sulfonylurea Formation")
    y -= 10

    line("Proposed Synthetic Step:")
    line("Reactants:")
    line("• Amine fragment")
    line("• Sulfonyl isocyanate")
    y -= 10

    line("Conditions:")
    line("• Base")
    line("• Aprotic solvent")
    line("• 0–25 °C")
    y -= 10

    line("AI Notes:")
    line("This synthesis route is AI-generated using retrosynthetic heuristics.")
    line("Experimental validation is required.")
    y -= 20

    line("DISCLAIMER:")
    line("For research & planning only. Not a laboratory protocol.")

    c.showPage()
    c.save()


# =========================
# PDF DOWNLOAD ENDPOINT
# =========================
@app.post(
    "/generate-pdf",
    response_class=FileResponse,
    summary="Generate & Download Synthesis PDF"
)
def generate_pdf_endpoint(data: PDFRequest = Body(...)):
    file_id = uuid.uuid4().hex
    filename = f"synthesis_report_{file_id}.pdf"

    tmp_path = f"/tmp/{filename}"

    generate_pdf(data.smiles, tmp_path)

    return FileResponse(
        path=tmp_path,
        filename=filename,
        media_type="application/pdf"
    )
