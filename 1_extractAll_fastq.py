#-----------------------------------------------------------------------------
# 1_extractAll_fastq.py
#-----------------------------------------------------------------------------
# Author:  reza.mozafari
# Date:    2025-03-06
#
# Summary: Reads gzipped FASTQ files (one per selection round), streams R1
#          sequence lines, and writes one DNA sequence per line to data/output/.
#          Used as step 1 of the mRNA-display NGS pipeline (rounds 1, 5, 7, 8).
#-----------------------------------------------------------------------------
# Input:   data/FASTQ/{round}_R1_001.fastq.gz (e.g. 1_R1_001.fastq.gz)
# Output:  data/output/round{N}_DNA.txt (one sequence per line)
#-----------------------------------------------------------------------------

import gzip
import os

#-------------------------------------------------------------------------------------
# CONFIG
#-------------------------------------------------------------------------------------
FASTQ_DIR = "data/FASTQ"
OUTPUT_DIR = "data/output"
ROUNDS = [1, 5, 7, 8]
READ = "R1"                    # which read has the construct (R1 or R2)
LANE = "001"                   # lane number in filename
PROGRESS_EVERY = 1_000_000     # print progress every N reads (set to 0 to disable)

#-------------------------------------------------------------------------------------

os.makedirs(OUTPUT_DIR, exist_ok=True)

for round_num in ROUNDS:
    fastq_name = f"{round_num}_{READ}_{LANE}.fastq.gz"
    fastq_path = os.path.join(FASTQ_DIR, fastq_name)
    out_path = os.path.join(OUTPUT_DIR, f"round{round_num}_DNA.txt")

    if not os.path.isfile(fastq_path):
        print(f"Skip round {round_num}: not found {fastq_path}")
        continue

    print(f"Processing round {round_num}: {fastq_path} -> {out_path}")

    open_fn = gzip.open if fastq_path.endswith(".gz") else open
    mode = "rt" if fastq_path.endswith(".gz") else "r"

    with open_fn(fastq_path, mode) as infile, open(out_path, "w") as outfile:
        i = -1
        n_reads = 0
        for line in infile:
            i += 1
            if i % 4 == 0:
                outfile.write(line.strip() + "\n")
                n_reads += 1
                if PROGRESS_EVERY and n_reads % PROGRESS_EVERY == 0:
                    print(f"  round {round_num}: {n_reads} reads")

    print(f"  round {round_num}: done, {n_reads} reads written to {out_path}")