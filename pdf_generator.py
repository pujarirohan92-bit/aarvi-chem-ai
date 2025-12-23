from fpdf import FPDF
import os
import uuid

def generate_pdf(report_text: str) -> str:
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_font("Arial", size=10)

    for line in report_text.split("\n"):
        pdf.multi_cell(0, 6, line)

    filename = f"synthesis_report_{uuid.uuid4().hex}.pdf"
    file_path = os.path.join("/tmp", filename)

    pdf.output(file_path)
    return file_path, filename
