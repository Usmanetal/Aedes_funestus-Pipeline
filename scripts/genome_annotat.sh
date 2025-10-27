# Create a directory
mkdir -p data/anopheles_funestus && cd data/anopheles_funestus

# Genome FASTA
wget ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/release-59/metazoa/fasta/anopheles_funestus/dna/Anopheles_funestus.AfunF3.dna.toplevel.fa.gz

# Gene annotation GTF
wget ftp://ftp.ensemblgenomes.ebi.ac.uk/pub/release-59/metazoa/gtf/anopheles_funestus/Anopheles_funestus.AfunF3.augustus.gtf.gz
