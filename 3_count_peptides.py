#-----------------------------------------------------------------------------
# 3_count_peptides.py
#-----------------------------------------------------------------------------
# Author:  reza.mozafari
# Date:    2025-03-06
#
# Summary: Reads translated per-round files from Script 2, counts each unique
#          DNA sequence (and its peptide), computes count, frequency, and rank.
#          Writes one CSV per round for Script 4 (merge and enrichment).
#-----------------------------------------------------------------------------
# Input:   data/output/round{N}_translated.txt (from 2_extract_translate.py)
# Output:  data/output/round{N}_counts.csv (dna_seq, peptide, count, frequency, rank)
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
PROGRESS_EVERY = cfg.get("progress_every", 0)


for round_num in ROUNDS:
    trans_path = os.path.join(OUTPUT_DIR, f"round{round_num}_translated.txt")
    out_path = os.path.join(OUTPUT_DIR, f"round{round_num}_counts.csv")

    if not os.path.isfile(trans_path):
        print(f"Skip round {round_num}: not found {trans_path}")
        continue

    print(f"Processing round {round_num}: {trans_path} -> {out_path}")

    # count[dna_seq] = number of reads with this dna_seq; keep one peptide per dna
    count_by_dna = {}
    dna_to_peptide = {}
    n_lines = 0

    with open(trans_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            dna = row.get("dna_seq", "").strip()
            pep = row.get("peptide", "").strip()
            if not dna:
                continue
            n_lines += 1
            count_by_dna[dna] = count_by_dna.get(dna, 0) + 1
            dna_to_peptide[dna] = pep
            if PROGRESS_EVERY and n_lines % PROGRESS_EVERY == 0:
                print(f"  round {round_num}: {n_lines} lines read")

    total = sum(count_by_dna.values())
    # Sort by count descending, then assign rank
    sorted_dna = sorted(count_by_dna.keys(), key=lambda d: -count_by_dna[d])

    with open(out_path, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(["dna_seq", "peptide", "count", "frequency", "rank"])
        for rank, dna in enumerate(sorted_dna, start=1):
            c = count_by_dna[dna]
            freq = c / total if total else 0
            writer.writerow([dna, dna_to_peptide[dna], c, f"{freq:.6e}", rank])

    print(f"  round {round_num}: done, {len(count_by_dna)} unique sequences -> {out_path}")