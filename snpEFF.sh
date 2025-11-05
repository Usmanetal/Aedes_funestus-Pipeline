mkdir -p snpEff/data/Anopheles_funestus

cp data/genome/Anopheles_funestus.AfunF3.dna.toplevel.fa \
   snpEff/data/Anopheles_funestus/sequences.fa

cp data/reference/Anopheles_funestus.AfunF3.augustus.gtf \
   snpEff/data/Anopheles_funestus/genes.gtf


mkdir -p snpEff
cd snpEff

unzip snpEff_latest_core.zip

cd snpEff
java -jar snpEff/snpEff.jar download -v Anopheles_funestus


for f in bcftools_scikit/results/vcf/*.vcf.gz; do
    base=$(basename $f .vcf.gz)
    java -jar snpEff/snpEff.jar -v Anopheles_funestus $f > bcftools_scikit/results/vcf/${base}.ann.vcf
done

# summary stats for one sample
java -jar snpEff/snpEff.jar -v -stats GGAFUN_DELTA_2.html Anopheles_funestus \
  bcftools_scikit/results/vcf/GGAFUN_DELTA_2.ann.vcf > GGAFUN_DELTA_2.snpEff.vcf

# Allele | Effect | Impact | Gene_Name | Gene_ID | Feature_Type | Transcript_ID | Biotype | Rank/Total | cDNA_change | Protein_change | CDS_pos/Length | AA_pos/Length | Distance | ... 

java -jar snpEff/snpEff.jar -v -stats FANG.ann.html Anopheles_funestus \
  bcftools_scikit/results/vcf/FANG.ann.vcf > FANG.ann.snpEff.vcf