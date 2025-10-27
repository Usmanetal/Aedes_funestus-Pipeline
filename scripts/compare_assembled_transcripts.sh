
mkdir -p results/stringtie_merged

stringtie --merge \
  -p 8 \
  -G data/reference/Anopheles_funestus.AfunF3.augustus.gtf \
  -o results/stringtie_merged/stringtie_merged.gtf \
  results/stringtie/*/*.gtf

mkdir -p results/gffcompare
gffcompare -r data/reference/Anopheles_funestus.AfunF3.augustus.gtf \
           -o results/gffcompare/merged \
           results/stringtie_merged/stringtie_merged.gtf

