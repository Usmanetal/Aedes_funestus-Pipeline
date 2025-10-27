# mkdir -p results/alignment_qc

# for bam in results/hisat2_alignments/*.sorted.bam; do
#   sample=$(basename $bam .sorted.bam)
#   samtools flagstat $bam > results/alignment_qc/${sample}_flagstat.txt
# done

for bam in results/hisat2_alignments/*.sorted.bam; do
  sample=$(basename $bam .sorted.bam)
  samtools stats $bam > results/alignment_qc/${sample}_stats.txt
  plot-bamstats -p results/alignment_qc/${sample}_plots/ results/alignment_qc/${sample}_stats.txt
done
