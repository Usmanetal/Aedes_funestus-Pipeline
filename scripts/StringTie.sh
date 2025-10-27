# 1️⃣ Decompress the annotation file
gunzip -c data/reference/Anopheles_funestus.AfunF3.augustus.gtf.gz > data/reference/Anopheles_funestus.AfunF3.augustus.gtf

# 2️⃣ Create the output directory for StringTie results
mkdir -p results/stringtie

# 3️⃣ Run StringTie for each sorted BAM file
for bam in results/hisat2_alignments/*.sorted.bam; do
  sample=$(basename "$bam" .sorted.bam)
  mkdir -p "results/stringtie/${sample}"
  
  stringtie "$bam" -p 8 \
    -G data/reference/Anopheles_funestus.AfunF3.augustus.gtf \
    -e -B \
    -o "results/stringtie/${sample}/${sample}.gtf"
done
