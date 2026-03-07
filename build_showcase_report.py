#-----------------------------------------------------------------------------
# build_showcase_report.py
#-----------------------------------------------------------------------------
# Builds a single PDF report for collaborators: figures, top-20 table,
# and concise interpretation text. Run after 5_report.py.
#-----------------------------------------------------------------------------
# Input:   data/output/report_figures.png, report_top20.csv, config.json
# Output:  data/output/Showcase_Report.pdf
#-----------------------------------------------------------------------------

import os
import csv
import json
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak,
)
from reportlab.lib.enums import TA_LEFT, TA_CENTER

_script_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_script_dir, "config.json")) as f:
    cfg = json.load(f)
OUTPUT_DIR = cfg["output_dir"]

fig_path = os.path.join(OUTPUT_DIR, "report_figures.png")
fig_pdf_path = os.path.join(OUTPUT_DIR, "report_figures.pdf")
top20_path = os.path.join(OUTPUT_DIR, "report_top20.csv")
out_pdf = os.path.join(OUTPUT_DIR, "Showcase_Report.pdf")

# Use PNG if present; else try to create it from PDF (so we don't require re-running 5_report)
if not os.path.isfile(fig_path) and os.path.isfile(fig_pdf_path):
    try:
        import fitz  # PyMuPDF
        doc = fitz.open(fig_pdf_path)
        page = doc[0]
        pix = page.get_pixmap(dpi=150)
        pix.save(fig_path)
        doc.close()
    except Exception:
        pass

styles = getSampleStyleSheet()
title_style = ParagraphStyle(
    name="ReportTitle",
    parent=styles["Heading1"],
    fontSize=16,
    spaceAfter=12,
    alignment=TA_CENTER,
)
heading_style = ParagraphStyle(
    name="Section",
    parent=styles["Heading2"],
    fontSize=12,
    spaceBefore=14,
    spaceAfter=6,
)
body_style = ParagraphStyle(
    name="Body",
    parent=styles["Normal"],
    fontSize=10,
    spaceAfter=8,
)

doc = SimpleDocTemplate(
    out_pdf,
    pagesize=letter,
    rightMargin=0.75 * inch,
    leftMargin=0.75 * inch,
    topMargin=0.75 * inch,
    bottomMargin=0.75 * inch,
)
story = []

# ---- Title and introduction ----
story.append(Paragraph(
    "mRNA-Display Peptide Selection — NGS QC Report",
    title_style,
))
story.append(Paragraph(
    "Here’s a quick pass over the NGS data from the mRNA-display selection (rounds 1, 5, 7, 8). "
    "The plots give a sanity check on read counts and point to the top hits by abundance and enrichment; "
    "the table at the end has the full peptide sequences for those.",
    body_style,
))
story.append(Spacer(1, 0.2 * inch))

# ---- Figure ----
if os.path.isfile(fig_path):
    img = Image(fig_path, width=6.5 * inch, height=2.8 * inch)
    story.append(img)
    story.append(Spacer(1, 0.15 * inch))
else:
    story.append(Paragraph(
        "<i>Figure not found. Run 5_report.py first to generate report_figures.png.</i>",
        body_style,
    ))

# ---- Interpretation ----
story.append(Paragraph("What the figures show", heading_style))

story.append(Paragraph("<b>1. Reads per round (left)</b>", body_style))
story.append(Paragraph(
    "Just total reads (millions) per round. We see fewer reads in later rounds because we’re "
    "sequencing what came through selection, not the whole library—so the drop is normal. "
    "Depth looks fine at each step.",
    body_style,
))

story.append(Paragraph("<b>2. Top 20 by count, round 1 (middle)</b>", body_style))
story.append(Paragraph(
    "The 20 most common peptides in the starting library. Useful as a baseline: if something "
    "was rare here but shows up strong in the enrichment panel, that’s a good candidate to chase.",
    body_style,
))

story.append(Paragraph("<b>3. Top 20 by enrichment, R8 vs R1 (right)</b>", body_style))
story.append(Paragraph(
    "Enrichment = (frequency in round 8) / (frequency in round 1). These are the sequences that "
    "gained the most over the selection—the ones to prioritize for synthesis and binding checks.",
    body_style,
))

story.append(Spacer(1, 0.2 * inch))

# ---- Summary table ----
story.append(Paragraph("Top 20 peptides — full sequences", heading_style))

if os.path.isfile(top20_path):
    rows = []
    with open(top20_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            panel = row.get("panel", "")
            rank = row.get("rank", "")
            peptide = row.get("peptide_full", "")
            c1 = row.get("count_r1", "")
            enr = row.get("enrichment_r8_vs_r1", "")
            if enr and float(enr) != 0:
                try:
                    enr = f"{float(enr):,.0f}"
                except ValueError:
                    pass
            rows.append([panel, rank, peptide, c1, enr])
    if rows:
        head = ["Panel", "Rank", "Peptide (full)", "Count R1", "Enrichment R8/R1"]
        data = [head] + rows
        t = Table(data, colWidths=[1.4 * inch, 0.5 * inch, 2.2 * inch, 0.7 * inch, 1.0 * inch])
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#4472C4")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 9),
            ("FONTSIZE", (0, 1), (-1, -1), 8),
            ("ALIGN", (0, 0), (1, -1), "LEFT"),
            ("ALIGN", (2, 0), (2, -1), "LEFT"),
            ("ALIGN", (3, 0), (-1, -1), "RIGHT"),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
            ("TOPPADDING", (0, 0), (-1, 0), 8),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F2F2F2")]),
        ]))
        story.append(t)
    else:
        story.append(Paragraph("No data in report_top20.csv.", body_style))
else:
    story.append(Paragraph(
        "<i>report_top20.csv not found. Run 5_report.py first.</i>",
        body_style,
    ))

story.append(Spacer(1, 0.25 * inch))

# ---- Notation ----
story.append(Paragraph("Notation", heading_style))
story.append(Paragraph(
    "A * in a sequence is a stop (or ambiguous call) from the NNK codon.",
    body_style,
))

doc.build(story)
print(f"Showcase report saved -> {out_pdf}")
