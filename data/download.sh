#!/bin/bash
set -euo pipefail

BASE_URL="https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed"
OUTDIR="data/raw_fastq"
mkdir -p "$OUTDIR"

echo "🔍 Crawling site for paired-end FASTQ files from: $BASE_URL"

# Step 1: Crawl the server without downloading
wget --spider -r -np -nd "$BASE_URL" 2>&1 | \
grep '^--' | awk '{print $3}' > all_urls.txt

# Step 2: Keep only paired-end FASTQs (_R1_ or _R2_)
grep -E '_R[12]_.*\.fastq\.gz$' all_urls.txt > paired_fastq_urls.txt

# Step 3: Exclude CHAD samples
grep -Ev 'CHAD' paired_fastq_urls.txt > wanted_fastq_urls.txt

echo "✅ Found $(wc -l < paired_fastq_urls.txt) total paired FASTQs"
echo "✅ Keeping $(wc -l < wanted_fastq_urls.txt) non-CHAD FASTQs"

# Step 4: Download them
echo "⬇️  Downloading selected FASTQ files to $OUTDIR ..."
wget -c -P "$OUTDIR" -i wanted_fastq_urls.txt

echo "🎉 Done! Paired-end FASTQs are in $OUTDIR/"
