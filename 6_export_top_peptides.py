#-----------------------------------------------------------------------------
# 6_export_top_peptides.py
#-----------------------------------------------------------------------------
# Reads merged_counts.csv, filters out stop-codon peptides, ranks by final-round
# count (e.g. R8), exports top N peptides to CSV and Excel.
# Run after 4_merge_rounds.py.
#-----------------------------------------------------------------------------
# Input:   data/output/merged_counts.csv, config.json
# Output:  data/output/top100_peptides.csv, data/output/top100_peptides.xlsx
#-----------------------------------------------------------------------------

import os
import json
import csv

_script_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_script_dir, "config.json")) as f:
    cfg = json.load(f)

OUTPUT_DIR = cfg["output_dir"]
ROUNDS = cfg["rounds"]
TOP_N = cfg.get("top_n_peptides", 100)

merged_path = os.path.join(OUTPUT_DIR, "merged_counts.csv")
out_csv = os.path.join(OUTPUT_DIR, f"top{TOP_N}_peptides.csv")
out_xlsx = os.path.join(OUTPUT_DIR, f"top{TOP_N}_peptides.xlsx")

if not os.path.isfile(merged_path):
    print(f"Error: {merged_path} not found. Run 4_merge_rounds.py first.")
    exit(1)

# Build output column names from config rounds
count_cols = [f"count_r{r}" for r in ROUNDS]
freq_cols = [f"freq_r{r}" for r in ROUNDS]
last_round = ROUNDS[-1]
count_last = f"count_r{last_round}"

# Read all rows, skip stop-codon peptides
rows = []
with open(merged_path) as f:
    reader = csv.DictReader(f)
    for row in reader:
        pep = (row.get("peptide") or "").strip()
        if "*" in pep:
            continue
        rows.append(row)

# Sort by final-round count descending, take top N
def get_count(r, row):
    try:
        return int(row.get(f"count_r{r}", 0) or 0)
    except (ValueError, TypeError):
        return 0

rows.sort(key=lambda r: get_count(last_round, r), reverse=True)
top = rows[:TOP_N]

# Output column order: peptide, R1_count, R5_count, ..., freq_R1, ..., enrichment, rank
out_count_names = [f"R{r}_count" for r in ROUNDS]
out_freq_names = [f"freq_R{r}" for r in ROUNDS]
out_header = ["peptide"] + out_count_names + out_freq_names + ["enrichment_R8_R1", f"rank_R{last_round}"]

# Write CSV
with open(out_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(out_header)
    for rank, row in enumerate(top, 1):
        pep = (row.get("peptide") or "").strip()
        counts = [str(row.get(c, "")) for c in count_cols]
        freqs = [str(row.get(c, "")) for c in freq_cols]
        enr = (row.get("enrichment_r8_vs_r1") or "").strip()
        w.writerow([pep] + counts + freqs + [enr, rank])
print(f"Wrote {out_csv}")

# Write Excel if openpyxl available
try:
    import openpyxl
    from openpyxl.styles import Font, Alignment
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Top peptides"
    for c, h in enumerate(out_header, 1):
        ws.cell(row=1, column=c, value=h)
        ws.cell(row=1, column=c).font = Font(bold=True)
    for r, row in enumerate(top, 2):
        pep = (row.get("peptide") or "").strip()
        counts = [row.get(c, "") for c in count_cols]
        freqs = [row.get(c, "") for c in freq_cols]
        enr = (row.get("enrichment_r8_vs_r1") or "").strip()
        for c, val in enumerate([pep] + counts + freqs + [enr, r - 1], 1):
            ws.cell(row=r, column=c, value=val)
    wb.save(out_xlsx)
    print(f"Wrote {out_xlsx}")
except ImportError:
    print("openpyxl not installed; skipping Excel. Install with: pip install openpyxl")