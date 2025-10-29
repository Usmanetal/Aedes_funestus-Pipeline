############################################################
# RNA-seq Pipeline: Alignment and Feature Counting with Rsubread
# Author: Ibrahim's Lab
# Date: August 31, 2025
############################################################

# 1. Install and load Rsubread
# if (!requireNamespace("BiocManager", quietly = TRUE))
# 		install.packages("BiocManager")
# if (!requireNamespace("Rsubread", quietly = TRUE))
# 		BiocManager::install("Rsubread")
library(Rsubread)

# 2. Create output directories
dir.create("rsubread", showWarnings = FALSE)
dir.create("rsubread/bam", showWarnings = FALSE)
dir.create("rsubread/sam", showWarnings = FALSE)
dir.create("rsubread/counts", showWarnings = FALSE)
dir.create("rsubread/index", showWarnings = FALSE)

# 3. Build genome index (only needs to be done once)
# Replace 'genome.fasta' with your actual fasta filename in 'ref/'
buildindex(
	basename = "rsubread/index/genome_index",
	reference = "data/genome/Anopheles_funestus.AfunF3.dna.toplevel.fa"
)


# 4. List paired-end FASTQ files in HADRES/
fq_dir <- "/mnt/d/fastq_katsina"

fq1 <- list.files(fq_dir, pattern = "_R1_001\\.fastq\\.gz$", full.names = TRUE)
fq2 <- list.files(fq_dir, pattern = "_R2_001\\.fastq\\.gz$", full.names = TRUE)

# Ensure matching pairs by sample name
get_sample_name <- function(x) sub("_R1\\.fastq.gz$|_R2\\.fastq.gz$", "", basename(x))
samples <- intersect(get_sample_name(fq1), get_sample_name(fq2))
fq1 <- fq1[get_sample_name(fq1) %in% samples]
fq2 <- fq2[get_sample_name(fq2) %in% samples]
fq1 <- fq1[order(get_sample_name(fq1))]
fq2 <- fq2[order(get_sample_name(fq2))]

bam_files <- file.path("rsubread/bam", paste0(samples, ".BAM"))
sam_files <- file.path("rsubread/sam", paste0(samples, ".SAM"))

# 5. Align paired-end reads (output both BAM and SAM)
for (i in seq_along(samples)) {
	align(
		index = "rsubread/index/genome_index",
		readfile1 = fq1[i],
		readfile2 = fq2[i],
		output_file = bam_files[i],
		output_format = "BAM",
		nthreads = 8
	)
	align(
		index = "rsubread/index/genome_index",
		readfile1 = fq1[i],
		readfile2 = fq2[i],
		output_file = sam_files[i],
		output_format = "SAM",
		nthreads = 8
	)
}


# 6. Feature counting (paired-end)
# Replace 'genes.gtf' with your actual GTF filename in 'ref/'
fc <- featureCounts(
	files = bam_files,
	annot.ext = "data/reference/Anopheles_funestus.AfunF3.augustus.gtf",
	isGTFAnnotationFile = TRUE,
	GTF.featureType = "exon",
	GTF.attrType = "gene_id",
	useMetaFeatures = TRUE,
	isPairedEnd = TRUE,
	nthreads = 8
)

# 7. Save count matrix
write.csv(fc$counts, file = "rsubread/counts/gene_counts.csv")

# ---
# For paired-end data, use:
# align(..., readfile1 = <R1>, readfile2 = <R2>, ...)
# featureCounts(..., isPairedEnd = TRUE, ...)
 
# ---
# This script will:
# - Build an index from your reference genome (in rsubread/index/)
# - Align all FASTQ files in HADRES/ to the reference, outputting BAM (rsubread/bam/) and SAM (rsubread/sam/)
# - Count reads per gene using your GTF annotation
# - Save the count matrix as rsubread/counts/gene_counts.csv
#
# ---
# For communication and reporting, you can copy this code into an Rmarkdown (.Rmd) file.
# Rmarkdown allows you to combine code, results, and narrative text for reproducible reports.
# ---
# Make sure to update the filenames if your fasta/gtf are named differently.

counts_df <- read.table("rsubread/counts/gene_counts.txt", header = TRUE, sep = "\t", row.names = 1)