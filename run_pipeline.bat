@echo off
cd /d "%~dp0"

echo ==============================================
echo Step 1/6: Extracting sequences from FASTQ
echo ==============================================
python 1_extract_fastq.py
if errorlevel 1 exit /b 1

echo.
echo ==============================================
echo Step 2/6: Finding motif and translating to peptide
echo ==============================================
python 2_extract_translate.py
if errorlevel 1 exit /b 1

echo.
echo ==============================================
echo Step 3/6: Counting peptides per round
echo ==============================================
python 3_count_peptides.py
if errorlevel 1 exit /b 1

echo.
echo ==============================================
echo Step 4/6: Merging rounds and computing enrichment
echo ==============================================
python 4_merge_rounds.py
if errorlevel 1 exit /b 1

echo.
echo ==============================================
echo Step 5/6: Generating report figures and top-20 table
echo ==============================================
python 5_report.py
if errorlevel 1 exit /b 1

echo.
echo ==============================================
echo Step 6/6: Building showcase report PDF
echo ==============================================
python build_showcase_report.py
if errorlevel 1 exit /b 1

echo.
echo ==============================================
echo Pipeline complete. Check data\output\ for results.
echo ==============================================
pause
