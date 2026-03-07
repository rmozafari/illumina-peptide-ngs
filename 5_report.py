#-----------------------------------------------------------------------------
# 5_report.py
#-----------------------------------------------------------------------------
# Author:  reza.mozafari
# Date:    2025-03-06
#
# Summary: Reads pipeline outputs and generates a short visual QC report
#          (reads per round, top sequences by count/enrichment).
#-----------------------------------------------------------------------------
# Input:   data/output/QC_summary.txt, data/output/merged_counts.csv (or round*_counts.csv)
# Output:  data/output/report_figures.pdf
#-----------------------------------------------------------------------------

import os
import re
import json
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_script_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_script_dir, "config.json")) as f:
    cfg = json.load(f)
OUTPUT_DIR = cfg["output_dir"]
ROUNDS = cfg["rounds"]

qc_path = os.path.join(OUTPUT_DIR, "QC_summary.txt")
merged_path = os.path.join(OUTPUT_DIR, "merged_counts.csv")
out_pdf = os.path.join(OUTPUT_DIR, "report_figures.pdf")

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

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

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

# 2) Top 20 by count from merged (stream first N rows)
top_n = 20
counts_r1 = []
labels = []
if os.path.isfile(merged_path):
    with open(merged_path) as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i >= top_n:
                break
            try:
                c = int(row.get("count_r1", 0))
            except (ValueError, TypeError):
                c = 0
            counts_r1.append(c)
            pep = (row.get("peptide", "") or "")[:12]
            labels.append(pep + ("..." if len(row.get("peptide", "")) > 12 else ""))
    if counts_r1:
        ax = axes[1]
        ax.barh(range(len(counts_r1)), counts_r1, color="coral", edgecolor="black")
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Count (round 1)")
        ax.set_title("Top 20 sequences by count (R1)")
    else:
        axes[1].text(0.5, 0.5, "No merged data", ha="center", va="center")
else:
    axes[1].text(0.5, 0.5, "merged_counts.csv not found", ha="center", va="center")

plt.tight_layout()
plt.savefig(out_pdf, format="pdf", bbox_inches="tight")
plt.close()
print(f"Report saved -> {out_pdf}")