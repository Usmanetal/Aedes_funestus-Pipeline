mkdir -p resources/featurecounts/

featureCounts -T 8 \
-p \
-s 2 \
-a data/reference/Anopheles_funestus.AfunF3.augustus.gtf \
-o resources/featurecounts/gene_counts.txt \
-g gene_id \
-t exon \
results/hisat2_alignments/*.bam
