#-----------------------------------------------------------------------------
# 5_report.py
#-----------------------------------------------------------------------------
# Author:  reza.mozafari
# Date:    2025-03-06
#
# Summary: Reads pipeline outputs and generates a short visual QC report
#          (reads per round, top 20 by count R1, top 20 by enrichment R8/R1).
#          Writes report_top20.csv with full peptide sequences for the report.
#-----------------------------------------------------------------------------
# Input:   data/output/QC_summary.txt, round1_counts.csv, merged_counts.csv
# Output:  data/output/report_figures.pdf, report_figures.png, report_top20.csv
#-----------------------------------------------------------------------------

import os
import re
import json
import csv
import heapq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_script_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_script_dir, "config.json")) as f:
    cfg = json.load(f)
OUTPUT_DIR = cfg["output_dir"]
ROUNDS = cfg["rounds"]

qc_path = os.path.join(OUTPUT_DIR, "QC_summary.txt")
round1_counts_path = os.path.join(OUTPUT_DIR, "round1_counts.csv")
merged_path = os.path.join(OUTPUT_DIR, "merged_counts.csv")
out_pdf = os.path.join(OUTPUT_DIR, "report_figures.pdf")
out_top20_csv = os.path.join(OUTPUT_DIR, "report_top20.csv")

# Parse QC summary for reads per round
round_reads = {}
round_unique = {}
if os.path.isfile(qc_path):
    with open(qc_path) as f:
        for line in f:
            m = re.search(r"Round (\d+): (\d+) total reads, (\d+) unique", line)
            if m:
                r, reads, uniq = int(m.group(1)), int(m.group(2)), int(m.group(3))
                round_reads[r] = reads
                round_unique[r] = uniq

top_n = 20
fig, axes = plt.subplots(1, 3, figsize=(14, 6))

# 1) Bar chart: total reads per round
if round_reads:
    ax = axes[0]
    rounds = sorted(round_reads.keys())
    reads = [round_reads[r] / 1e6 for r in rounds]
    ax.bar([str(r) for r in rounds], reads, color="steelblue", edgecolor="black")
    ax.set_xlabel("Round")
    ax.set_ylabel("Total reads (millions)")
    ax.set_title("Reads per round")
else:
    axes[0].text(0.5, 0.5, "No QC data", ha="center", va="center")
    axes[0].set_title("Reads per round")

# 2) Top 20 by count R1 (from round1_counts.csv); collect full peptide + count for report file
counts_r1 = []
labels_count = []
top20_count_rows = []  # (full_peptide, count) for CSV
if os.path.isfile(round1_counts_path):
    with open(round1_counts_path) as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= top_n:
                break
            try:
                c = int(row.get("count", 0))
            except (ValueError, TypeError):
                c = 0
            pep = (row.get("peptide", "") or "").strip()
            counts_r1.append(c)
            labels_count.append(pep)  # full peptide (13 aa)
            top20_count_rows.append((pep, c))
if counts_r1:
    ax = axes[1]
    ax.barh(range(len(counts_r1)), counts_r1, color="coral", edgecolor="black")
    ax.set_yticks(range(len(labels_count)))
    ax.set_yticklabels(labels_count, fontsize=7)
    ax.set_xlabel("Count (round 1)")
    ax.set_title("Top 20 sequences by count (R1)")
else:
    axes[1].text(0.5, 0.5, "round1_counts.csv not found or empty", ha="center", va="center")
    axes[1].set_title("Top 20 by count (R1)")

# 3) Top 20 by enrichment R8/R1; collect full peptide + enrichment for report file
top_enrich = []
top20_enrich_rows = []  # (full_peptide, enrichment) for CSV
if os.path.isfile(merged_path):
    with open(merged_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw = row.get("enrichment_r8_vs_r1", "").strip()
            try:
                enr = float(raw) if raw else 0.0
            except (ValueError, TypeError):
                enr = 0.0
            if enr <= 0:
                continue
            pep = (row.get("peptide", "") or "").strip()
            if len(top_enrich) < top_n:
                heapq.heappush(top_enrich, (enr, pep, enr))
            elif enr > top_enrich[0][0]:
                heapq.heapreplace(top_enrich, (enr, pep, enr))
    top_enrich.sort(key=lambda x: -x[0])
    top20_enrich_rows = [(x[1], x[2]) for x in top_enrich]
if top_enrich:
    ax = axes[2]
    vals = [x[2] for x in top_enrich]
    lbls = [x[1] for x in top_enrich]  # full peptide
    ax.barh(range(len(vals)), vals, color="seagreen", edgecolor="black")
    ax.set_yticks(range(len(lbls)))
    ax.set_yticklabels(lbls, fontsize=7)
    ax.set_xlabel("Enrichment (freq R8 / freq R1)")
    ax.set_title("Top 20 by enrichment (R8 vs R1)")
else:
    axes[2].text(0.5, 0.5, "No enrichment data or merged_counts.csv not found", ha="center", va="center")
    axes[2].set_title("Top 20 by enrichment (R8 vs R1)")

plt.tight_layout()
out_png = os.path.join(OUTPUT_DIR, "report_figures.png")
plt.savefig(out_pdf, format="pdf", bbox_inches="tight")
plt.savefig(out_png, format="png", bbox_inches="tight", dpi=150)
plt.close()
print(f"Report saved -> {out_pdf}")

# Write report_top20.csv: full peptide + count (top 20 R1) and full peptide + enrichment (top 20 R8 vs R1)
with open(out_top20_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["panel", "rank", "peptide_full", "count_r1", "enrichment_r8_vs_r1"])
    for rank, (pep, c) in enumerate(top20_count_rows, 1):
        w.writerow(["top20_count_R1", rank, pep, c, ""])
    for rank, (pep, enr) in enumerate(top20_enrich_rows, 1):
        w.writerow(["top20_enrichment_R8_vs_R1", rank, pep, "", enr])
print(f"Top-20 table (full peptides) -> {out_top20_csv}")