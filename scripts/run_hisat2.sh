#!/bin/bash
set -euo pipefail

# ================================
# HISAT2 alignment pipeline for Anopheles funestus
# ================================
# Requirements:
#   - hisat2
#   - samtools
#   - gzip
#
# Input:
#   data/genome/Anopheles_funestus.AfunF3.dna.toplevel.fa.gz
#   data/raw_fastq/*_R1_*.fastq.gz + *_R2_*.fastq.gz
#
# Output:
#   results/hisat2_alignments/*.sorted.bam
# ================================

# -------- CONFIG --------
GENOME_FA_GZ="data/genome/Anopheles_funestus.AfunF3.dna.toplevel.fa.gz"
GENOME_FA="data/genome/Anopheles_funestus.AfunF3.dna.toplevel.fa"
INDEX_PREFIX="data/genome/AfunF3_index"
READ_DIR="data/raw_fastq"
OUT_DIR="results/hisat2_alignments"
THREADS=8
# -------------------------

mkdir -p "$OUT_DIR"

echo "=== Aligning paired-end reads ==="

for R1 in "$READ_DIR"/*_R1_*.fastq.gz; do
    [ -e "$R1" ] || { echo "No R1 files found."; exit 1; }

    SAMPLE=$(basename "$R1" | sed 's/_R1_.*//')
    R2=$(echo "$R1" | sed 's/_R1_/_R2_/')
    OUT_BAM="${OUT_DIR}/${SAMPLE}.sorted.bam"

    if [ ! -f "$R2" ]; then
        echo "⚠️  Missing R2 for sample $SAMPLE — skipping"
        continue
    fi

    echo "🔹 Aligning $SAMPLE ..."
    {
        hisat2 -p "$THREADS" \
            -x "$INDEX_PREFIX" \
            -1 "$R1" -2 "$R2" \
            --rna-strandness RF \
            --dta \
        | samtools view -@ "$THREADS" -bS - \
        | samtools sort -@ "$THREADS" -o "$OUT_BAM"

        # Ensure sorting is complete before indexing
        wait

        echo "🔹 Indexing $OUT_BAM ..."
        samtools index "$OUT_BAM"

        echo "✅ Finished $SAMPLE → $(basename "$OUT_BAM")"
    } || {
        echo "❌ ERROR: Alignment failed for $SAMPLE — skipping this sample."
    }

    echo
done

echo "🎯 All alignments complete! BAMs stored in: $OUT_DIR"
