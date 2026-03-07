#-----------------------------------------------------------------------------
# 2_extract_translate.py
#-----------------------------------------------------------------------------
# Author:  reza.mozafari
# Date:    2025-03-06
#
# Summary: Reads DNA-per-line files from Script 1, finds motif ATGTGC + 30 nt + TGC,
#          extracts 39 nt, translates to 13 aa (M C X10 C). Writes per-round
#          translated tables for Script 3.
#-----------------------------------------------------------------------------
# Input:   data/output/round{N}_DNA.txt (from 1_extract_fastq.py)
# Output:  data/output/round{N}_translated.txt (dna_seq, peptide; one line per read)
#-----------------------------------------------------------------------------

import os
import json
from Bio.Seq import Seq

# Load config
_script_dir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_script_dir, "config.json")) as f:
    cfg = json.load(f)

OUTPUT_DIR = cfg["output_dir"]
ROUNDS = cfg["rounds"]
MOTIF_LEFT = cfg["motif_left"]
MOTIF_RIGHT = cfg["motif_right"]
VARIABLE_NT = cfg["variable_nt"]
MOTIF_TOTAL = cfg["motif_total"]
PROGRESS_EVERY = cfg.get("progress_every", 0)


def find_and_translate(line, try_rc=True):
    """
    Find ATGTGC + 30 nt + TGC in line (or its reverse complement).
    Returns (dna_39nt, peptide) or (None, None) if not found.
    Peptide may contain '*' (stop) from NNK design.
    """
    line = line.strip().upper()
    if not line:
        return None, None

    for seq in [line, str(Seq(line).reverse_complement())] if try_rc else [line]:
        pos = seq.find(MOTIF_LEFT)
        while pos != -1:
            start = pos
            end = start + MOTIF_TOTAL
            if end <= len(seq) and seq[end - len(MOTIF_RIGHT):end] == MOTIF_RIGHT:
                dna = seq[start:end]
                if len(dna) == MOTIF_TOTAL:
                    pep = str(Seq(dna).translate())
                    return dna, pep
            pos = seq.find(MOTIF_LEFT, pos + 1)
    return None, None


for round_num in ROUNDS:
    dna_path = os.path.join(OUTPUT_DIR, f"round{round_num}_DNA.txt")
    out_path = os.path.join(OUTPUT_DIR, f"round{round_num}_translated.txt")

    if not os.path.isfile(dna_path):
        print(f"Skip round {round_num}: not found {dna_path}")
        continue

    print(f"Processing round {round_num}: {dna_path} -> {out_path}")

    n_reads = 0
    n_valid = 0
    with open(dna_path) as infile, open(out_path, "w") as outfile:
        outfile.write("dna_seq\tpeptide\n")
        for line in infile:
            n_reads += 1
            if PROGRESS_EVERY and n_reads % PROGRESS_EVERY == 0:
                print(f"  round {round_num}: {n_reads} reads, {n_valid} valid")
            dna, pep = find_and_translate(line)
            if dna is not None:
                if "*" in pep:
                    continue   # skip stop-codon (truncated) peptides
                outfile.write(f"{dna}\t{pep}\n")
                n_valid += 1

    print(f"  round {round_num}: done, {n_valid} / {n_reads} valid -> {out_path}")