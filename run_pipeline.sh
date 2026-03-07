#!/usr/bin/env bash
set -e
cd "$(dirname "$0")"

echo "=============================================="
echo "Step 1/7: Extracting sequences from FASTQ"
echo "=============================================="
python 1_extract_fastq.py

echo ""
echo "=============================================="
echo "Step 2/7: Finding motif and translating to peptide"
echo "=============================================="
python 2_extract_translate.py

echo ""
echo "=============================================="
echo "Step 3/7: Counting peptides per round"
echo "=============================================="
python 3_count_peptides.py

echo ""
echo "=============================================="
echo "Step 4/7: Merging rounds and computing enrichment"
echo "=============================================="
python 4_merge_rounds.py

echo ""
echo "=============================================="
echo "Step 5/7: Exporting top N peptides (CSV)"
echo "=============================================="
python 6_export_top_peptides.py

echo ""
echo "=============================================="
echo "Step 6/7: Generating report figures and top-20 table"
echo "=============================================="
python 5_report.py

echo ""
echo "=============================================="
echo "Step 7/7: Building showcase report PDF"
echo "=============================================="
python build_showcase_report.py

echo ""
echo "=============================================="
echo "Pipeline complete. Check data/output/ for results."
echo "=============================================="
