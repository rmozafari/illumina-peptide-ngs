#-----------------------------------------------------------------------------
# 4_merge_rounds.py
#-----------------------------------------------------------------------------
# Author:  reza.mozafari
# Date:    2025-03-06
#
# Summary: Reads per-round count CSVs from Script 3, merges by dna_seq/peptide,
#          adds count and frequency per round and enrichment (e.g. freq_R8/freq_R1).
#          Single round = stream; multiple rounds = system sort then merge
#          (no full load into memory). Progress messages for long steps.
#-----------------------------------------------------------------------------
# Input:   data/output/round{N}_counts.csv (from 3_count_peptides.py)
# Output:  data/output/merged_counts.csv, data/output/QC_summary.txt
#-----------------------------------------------------------------------------

import os
import json
import csv
import subprocess
import tempfile
import heapq

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

# --- Single round: stream row-by-row ---
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
            out_row.append("")
            writer.writerow(out_row)
            n_rows += 1
    print(f"Merged {n_rows} sequences (streamed from round {r1}) -> {out_path}")
    with open(qc_path, "w") as qc:
        qc.write("QC summary\n==========\n")
        qc.write(f"Round {r1}: {total_reads} total reads, {n_rows} unique sequences\n")
        qc.write(f"\nMerged output: {out_path}\n")
    print(f"QC summary -> {qc_path}")
    exit(0)

# --- Multiple rounds: stream QC, then system sort + merge ---
# 1) QC: stream each file (no full load)
round_stats = {}
for r in rounds_present:
    path = os.path.join(OUTPUT_DIR, f"round{r}_counts.csv")
    total_reads = 0
    n_unique = 0
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                total_reads += int(row.get("count", 0))
            except (ValueError, TypeError):
                pass
            n_unique += 1
    round_stats[r] = {"total_reads": total_reads, "n_unique": n_unique}
    print(f"Round {r}: {total_reads} total reads, {n_unique} unique sequences")

# 2) Sort each round by dna_seq using system sort (external sort, low memory)
sorted_paths = []
for r in rounds_present:
    path = os.path.join(OUTPUT_DIR, f"round{r}_counts.csv")
    print(f"Sorting round {r} by dna_seq...")
    fd, sorted_path = tempfile.mkstemp(suffix=".csv", prefix=f"round{r}_", dir=OUTPUT_DIR)
    os.close(fd)
    with open(path) as inf:
        header_line = inf.readline()
        body_path = tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False, dir=OUTPUT_DIR, newline=""
        )
        for line in inf:
            body_path.write(line)
        body_path.close()
    subprocess.run(
        ["sort", "-t,", "-k1,1", "-o", sorted_path, body_path.name],
        check=True,
    )
    os.unlink(body_path.name)
    with open(sorted_path) as f:
        sorted_content = f.read()
    with open(sorted_path, "w") as f:
        f.write(header_line)
        f.write(sorted_content)
    sorted_paths.append((r, sorted_path))
    print(f"  Sorted round {r}.")

# 3) Merge sorted streams (heap)
def row_iter(round_path):
    with open(round_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            dna = row.get("dna_seq", "").strip()
            pep = row.get("peptide", "").strip()
            try:
                c = int(row.get("count", 0))
                freq = float(row.get("frequency", 0))
            except (ValueError, TypeError):
                c, freq = 0, 0.0
            yield (dna, pep, c, freq)

print("Merging sorted rounds...")
current = []
iters = []
for idx, (r, sp) in enumerate(sorted_paths):
    it = row_iter(sp)
    try:
        row = next(it)
        current.append((row[0], idx, row[1], row[2], row[3]))
        iters.append((it, r))
    except StopIteration:
        iters.append((iter([]), r))

heapq.heapify(current)
n_merged = 0
PROGRESS_EVERY = 5_000_000
data_by_round = {r: {"peptide": "", "count": 0, "frequency": 0.0} for r in ROUNDS}
prev_dna = None

with open(out_path, "w", newline="") as outfile:
    writer = csv.writer(outfile)
    writer.writerow(header)
    while current:
        dna, idx, pep, c, freq = heapq.heappop(current)
        r = sorted_paths[idx][0]
        data_by_round[r] = {"peptide": pep, "count": c, "frequency": freq}
        try:
            next_row = next(iters[idx][0])
            heapq.heappush(current, (next_row[0], idx, next_row[1], next_row[2], next_row[3]))
        except StopIteration:
            pass
        if prev_dna is not None and dna != prev_dna:
            peptide = ""
            for rr in ROUNDS:
                peptide = data_by_round[rr]["peptide"]
                if peptide:
                    break
            row = [prev_dna, peptide]
            for rr in ROUNDS:
                row.append(data_by_round[rr]["count"])
                row.append(f"{data_by_round[rr]['frequency']:.6e}" if data_by_round[rr]["frequency"] else "0")
            freq_r1 = data_by_round.get(ROUNDS[0], {}).get("frequency", 0) or 0
            freq_r8 = data_by_round.get(ROUNDS[-1], {}).get("frequency", 0) or 0
            row.append(safe_enrichment(freq_r8, freq_r1))
            writer.writerow(row)
            n_merged += 1
            if PROGRESS_EVERY and n_merged % PROGRESS_EVERY == 0:
                print(f"  Merged {n_merged} sequences...")
        prev_dna = dna
    if prev_dna is not None:
        peptide = ""
        for rr in ROUNDS:
            peptide = data_by_round[rr]["peptide"]
            if peptide:
                break
        row = [prev_dna, peptide]
        for rr in ROUNDS:
            row.append(data_by_round[rr]["count"])
            row.append(f"{data_by_round[rr]['frequency']:.6e}" if data_by_round[rr]["frequency"] else "0")
        freq_r1 = data_by_round.get(ROUNDS[0], {}).get("frequency", 0) or 0
        freq_r8 = data_by_round.get(ROUNDS[-1], {}).get("frequency", 0) or 0
        row.append(safe_enrichment(freq_r8, freq_r1))
        writer.writerow(row)
        n_merged += 1

for r, sp in sorted_paths:
    try:
        os.remove(sp)
    except OSError:
        pass

with open(qc_path, "w") as qc:
    qc.write("QC summary\n==========\n")
    for r in ROUNDS:
        if r in round_stats:
            s = round_stats[r]
            qc.write(f"Round {r}: {s['total_reads']} total reads, {s['n_unique']} unique sequences\n")
    qc.write(f"\nMerged unique sequences: {n_merged}\n")
    qc.write(f"Output: {out_path}\n")

print(f"Merged {n_merged} sequences -> {out_path}")
print(f"QC summary -> {qc_path}")