# Illumina NGS Pipeline — mRNA-Display Peptide Selection

Pipeline for analyzing next-generation sequencing (NGS) data from an mRNA-display peptide selection (rounds 1, 5, 7, 8). Extracts sequences from FASTQ, finds the library motif, translates to peptide, counts and merges across rounds, and produces figures and a showcase report.

## Library design

- **Motif:** `ATGTGC` + 30 nt (10× NNK) + `TGC` → 13 aa peptide (M C X₁₀ C).
- **Rounds:** 1, 5, 7, 8 (configurable in `config.json`).
- **Read:** R1 only.

## Requirements

- Python 3.8+
- Dependencies: `pip install -r requirements.txt`

  (biopython, matplotlib, reportlab; pymupdf optional for building the showcase report from an existing PDF.)

## Configuration

Edit **`config.json`** in the project root:

| Key | Description |
|-----|-------------|
| `fastq_dir` | Directory containing gzipped FASTQ files |
| `output_dir` | Where to write all outputs (default: `data/output`) |
| `rounds` | List of round numbers, e.g. `[1, 5, 7, 8]` |
| `read` | Read index, e.g. `"R1"` |
| `lane` | Lane, e.g. `"001"` |
| `motif_left` | Left motif (e.g. `"ATGTGC"`) |
| `motif_right` | Right motif (e.g. `"TGC"`) |
| `variable_nt` | Variable region length in nt (30) |
| `motif_total` | Total motif length in nt (39) |
| `progress_every` | Print progress every N reads (0 = off) |

**Expected FASTQ names:** `{round}_R1_001.fastq.gz` (e.g. `1_R1_001.fastq.gz`) inside `fastq_dir`.

## Directory layout

```
Illumina_Scripts/
├── config.json
├── requirements.txt
├── README.md
├── 1_extract_fastq.py
├── 2_extract_translate.py
├── 3_count_peptides.py
├── 4_merge_rounds.py
├── 5_report.py
├── build_showcase_report.py
├── data/
│   ├── FASTQ/          # input: {round}_R1_001.fastq.gz
│   └── output/         # all pipeline outputs
└── archive/            # older/alternate scripts (optional)
```

## Pipeline (run in order)

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1 | `1_extract_fastq.py` | `data/FASTQ/{round}_R1_001.fastq.gz` | `round{N}_DNA.txt` |
| 2 | `2_extract_translate.py` | `round{N}_DNA.txt` | `round{N}_translated.txt` |
| 3 | `3_count_peptides.py` | `round{N}_translated.txt` | `round{N}_counts.csv` |
| 4 | `4_merge_rounds.py` | `round{N}_counts.csv` | `merged_counts.csv`, `QC_summary.txt` |
| 5 | `5_report.py` | QC_summary, round1_counts, merged_counts | `report_figures.pdf`, `report_figures.png`, `report_top20.csv` |
| — | `build_showcase_report.py` | report_figures.png, report_top20.csv | `Showcase_Report.pdf` |

**Example (from project root):**

```bash
python 1_extract_fastq.py
python 2_extract_translate.py
python 3_count_peptides.py
python 4_merge_rounds.py
python 5_report.py
python build_showcase_report.py
```

## Main outputs

- **`QC_summary.txt`** — Total reads and unique sequences per round.
- **`merged_counts.csv`** — One row per unique sequence; columns for count/freq per round and enrichment (e.g. R8 vs R1).
- **`report_figures.pdf` / `report_figures.png`** — Three panels: reads per round, top 20 by count (R1), top 20 by enrichment (R8 vs R1).
- **`report_top20.csv`** — Full peptide sequences and metrics for the two top-20 lists (count R1 and enrichment R8 vs R1).
- **`Showcase_Report.pdf`** — Single report for sharing: figures, short interpretation, and top-20 table (full sequences). Use this for collaborators.

## Notation

In peptide sequences, **`*`** denotes a stop codon (or ambiguous base call) from the NNK codon.

## Author / reference

- Author: reza.mozafari  
- Date: 2025-03-06  
- For mRNA-display peptide selection NGS analysis.
