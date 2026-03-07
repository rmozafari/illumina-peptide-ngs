#-----------------------------------------------------------------------------
# 4_merge_rounds.py
#-----------------------------------------------------------------------------
# Author:  reza.mozafari
# Date:    2025-03-06
#
# Summary: Reads per-round count CSVs from Script 3, merges by dna_seq/peptide,
#          adds count and frequency per round and enrichment (e.g. freq_R8/freq_R1).
#          Streams when only one round is present to avoid memory issues.
#-----------------------------------------------------------------------------
# Input:   data/output/round{N}_counts.csv (from 3_count_peptides.py)
# Output:  data/output/merged_counts.csv, data/output/QC_summary.txt
#-----------------------------------------------------------------------------

import os
import json
import csv

# Load config
_script_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_script_dir, "config.json")) as f:
    cfg = json.load(f)

OUTPUT_DIR = cfg["output_dir"]
ROUNDS = cfg["rounds"]

# Check which round files exist
rounds_present = []
for r in ROUNDS:
    path = os.path.join(OUTPUT_DIR, f"round{r}_counts.csv")
    if os.path.isfile(path):
        rounds_present.append(r)
    else:
        print(f"Warning: {path} not found, skipping round {r}")

if not rounds_present:
    print("No round count files found. Exiting.")
    exit(1)

out_path = os.path.join(OUTPUT_DIR, "merged_counts.csv")
qc_path = os.path.join(OUTPUT_DIR, "QC_summary.txt")

header = ["dna_seq", "peptide"]
for r in ROUNDS:
    header.append(f"count_r{r}")
    header.append(f"freq_r{r}")
header.append("enrichment_r8_vs_r1")

def safe_enrichment(freq_last, freq_first):
    if freq_first and freq_first > 0:
        return f"{freq_last / freq_first:.4f}"
    return ""

# --- Single round: stream row-by-row (no big dict in memory) ---
if len(rounds_present) == 1:
    r1 = rounds_present[0]
    path = os.path.join(OUTPUT_DIR, f"round{r1}_counts.csv")
    total_reads = 0
    n_rows = 0
    with open(path) as infile, open(out_path, "w", newline="") as outfile:
        reader = csv.DictReader(infile)
        writer = csv.writer(outfile)
        writer.writerow(header)
        for row in reader:
            dna = row.get("dna_seq", "").strip()
            pep = row.get("peptide", "").strip()
            c1 = row.get("count", "0")
            f1 = row.get("frequency", "0")
            try:
                total_reads += int(c1)
                freq_r1_val = float(f1)
            except (ValueError, TypeError):
                freq_r1_val = 0.0
            out_row = [dna, pep]
            for r in ROUNDS:
                if r == r1:
                    out_row.append(c1)
                    out_row.append(f1 if f1 else "0")
                else:
                    out_row.append(0)
                    out_row.append("0")
            out_row.append("")  # no enrichment with one round
            writer.writerow(out_row)
            n_rows += 1
    print(f"Merged {n_rows} sequences (streamed from round {r1}) -> {out_path}")
    with open(qc_path, "w") as qc:
        qc.write("QC summary\n")
        qc.write("==========\n")
        qc.write(f"Round {r1}: {total_reads} total reads, {n_rows} unique sequences\n")
        qc.write(f"\nMerged output: {out_path}\n")
    print(f"QC summary -> {qc_path}")
    exit(0)

# --- Multiple rounds: load round data and merge (original logic) ---
round_data = {}
for r in rounds_present:
    path = os.path.join(OUTPUT_DIR, f"round{r}_counts.csv")
    round_data[r] = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            dna = row.get("dna_seq", "").strip()
            pep = row.get("peptide", "").strip()
            try:
                c = int(row.get("count", 0))
                freq = float(row.get("frequency", 0))
            except (ValueError, TypeError):
                c, freq = 0, 0.0
            round_data[r][dna] = {"peptide": pep, "count": c, "frequency": freq}

def get_data(r, dna):
    d = round_data.get(r, {}).get(dna, {"peptide": "", "count": 0, "frequency": 0.0})
    return d["peptide"], d["count"], d["frequency"]

all_dna = set()
for r in round_data:
    all_dna.update(round_data[r].keys())

with open(out_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    for dna in all_dna:
        peptide = ""
        for r in ROUNDS:
            peptide = get_data(r, dna)[0]
            if peptide:
                break
        row = [dna, peptide]
        for r in ROUNDS:
            _, c, freq = get_data(r, dna)
            row.append(c)
            row.append(f"{freq:.6e}" if freq else "0")
        freq_r1 = get_data(ROUNDS[0], dna)[2] if ROUNDS else 0
        freq_r8 = get_data(ROUNDS[-1], dna)[2] if ROUNDS else 0
        row.append(safe_enrichment(freq_r8, freq_r1))
        writer.writerow(row)
print(f"Merged {len(all_dna)} sequences -> {out_path}")

with open(qc_path, "w") as qc:
    qc.write("QC summary\n")
    qc.write("==========\n")
    for r in ROUNDS:
        if r not in round_data:
            continue
        total_reads = sum(d["count"] for d in round_data[r].values())
        n_unique = len(round_data[r])
        qc.write(f"Round {r}: {total_reads} total reads, {n_unique} unique sequences\n")
    qc.write(f"\nMerged unique sequences: {len(all_dna)}\n")
    qc.write(f"Output: {out_path}\n")
print(f"QC summary -> {qc_path}")