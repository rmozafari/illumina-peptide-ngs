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
import json

# Load config from config.json (same folder as this script)
_script_dir = os.path.dirname(os.path.abspath(__file__))
_config_path = os.path.join(_script_dir, "config.json")
with open(_config_path) as f:
    cfg = json.load(f)


os.makedirs(cfg["output_dir"], exist_ok=True)
progress_every = cfg.get("progress_every", 0)

for round_num in cfg["rounds"]:
    fastq_name = f"{round_num}_{cfg['read']}_{cfg['lane']}.fastq.gz"
    fastq_path = os.path.join(cfg["fastq_dir"], fastq_name)
    out_path = os.path.join(cfg["output_dir"], f"round{round_num}_DNA.txt")

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
            if i % 4 == 1:
                outfile.write(line.strip() + "\n")
                n_reads += 1
                if progress_every and n_reads % progress_every == 0:
                    print(f"  round {round_num}: {n_reads} reads")

    print(f"  round {round_num}: done, {n_reads} reads written to {out_path}")