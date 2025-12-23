from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4


def generate_pdf_report(output: dict, filename: str):
    doc = SimpleDocTemplate(filename, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []

    story.append(Paragraph("<b>AARVI CHEM AI — SYNTHESIS REPORT</b>", styles["Title"]))
    story.append(Spacer(1, 12))

    story.append(Paragraph(f"<b>Target SMILES:</b> {output['target_smiles']}", styles["Normal"]))
    story.append(Paragraph(f"<b>Overall Confidence:</b> {output['overall_confidence']}", styles["Normal"]))
    story.append(Spacer(1, 12))

    best_route = output["routes"][0]
    story.append(Paragraph("<b>Selected Best Route</b>", styles["Heading2"]))
    story.append(Paragraph(f"Route Score: {best_route['score']}", styles["Normal"]))
    story.append(Spacer(1, 8))

    for step in best_route["steps"]:
        story.append(Paragraph(f"<b>Step {step['step']} — {step['reaction']}</b>", styles["Heading3"]))

        table_data = [
            ["Reactants", ", ".join(step["reactants"])],
            ["Solvent", step["conditions"].get("solvent")],
            ["Base", step["conditions"].get("base")],
            ["Temperature", step["conditions"].get("temperature")],
            ["Feasibility", str(step["confidence"])],
            ["Expected Yield", str(step.get("expected_yield", "N/A"))],
            ["Warnings", ", ".join(step["warnings"]) if step["warnings"] else "None"]
        ]

        table = Table(table_data, colWidths=[120, 350])
        story.append(table)
        story.append(Spacer(1, 12))

    story.append(Paragraph("<b>Final Recommendation</b>", styles["Heading2"]))
    story.append(Paragraph(
        "Proceed with the selected route for laboratory synthesis. "
        "Steps with low confidence or yield may require optimization.",
        styles["Normal"]
    ))

    doc.build(story)
