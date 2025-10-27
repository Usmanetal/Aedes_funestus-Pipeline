mkdir -p results/alignment_qc

for bam in results/hisat2_alignments/*.sorted.bam; do
  sample=$(basename "$bam" .sorted.bam)
  echo "Running flagstat for $sample ..."
  samtools flagstat "$bam" > results/alignment_qc/${sample}_flagstat.txt
done
