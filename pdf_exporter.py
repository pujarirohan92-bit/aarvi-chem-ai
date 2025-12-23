from reportlab.platypus import SimpleDocTemplate, Paragraph
from reportlab.lib.styles import getSampleStyleSheet


def export_pdf(text: str, filename: str = "synthesis_report.pdf"):
    doc = SimpleDocTemplate(filename)
    styles = getSampleStyleSheet()
    story = []

    for line in text.split("\n"):
        story.append(Paragraph(line.replace("&", "&amp;"), styles["Normal"]))

    doc.build(story)
    return filename
