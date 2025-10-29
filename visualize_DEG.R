#Differential expression

#load the libraries




library(ggplot2)
library(ballgown)
library(genefilter)
library(RSkittleBrewer)
library(devtools)
library(ggrepel)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(gplots)
library(GenomicRanges)
library(viridis)
library(stringr)


#lets load the sample information
pheno_data <- read.csv("resources/sample_metadata.csv")

#let's show information for first 6 samples
head(pheno_data[pheno_data$description != "control", ])
pheno_data<- pheno_data[pheno_data$description != "control", ]

# Get a list of all StringTie sample directories
sample_dirs <- list.dirs("results/stringtie", full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[!grepl("GGAFUN_UNX", sample_dirs)]
sample_dirs

#Load the expression data using ballgown
bg_rs <- ballgown(
  samples = sample_dirs,
  pData = pheno_data,
  meas = "FPKM"
)
methods(class="ballgown")
#Lets filter out transcripts with low variance
#This is done to remove some genes that have few counts. Filtering improves the statistical power of differential expression analysis. 
#We use variance filter to remove transcripts with low variance( 1 or less)

bg_rs_filt<- subset(bg_rs,"rowVars(texpr(bg_rs))>1",genomesubset=TRUE)

# pData(bg_R-S_filt)= data.frame(
#   id=sampleNames(bg_chrX_filt),
#   group= rep(c("control","Resistant"),each=3))
# -----------------------------------------------------------
# Differential expression: transcript level (3 groups initially)
# -----------------------------------------------------------
de_transcripts <- stattest(
  bg_rs_filt,
  feature   = "transcript",
  covariate = "description",
  getFC     = TRUE,
  meas      = "FPKM"
)

#add identifiers
de_transcripts = data.frame(geneNames=ballgown::geneNames(bg_rs_filt), geneIDs=ballgown::geneIDs(bg_rs_filt), de_transcripts)

# Let's test on genes
de_genes <- stattest(bg_rs_filt,
          feature="gene",
          covariate="description",
          getFC=TRUE, meas="FPKM")


#lets get the gene names
bg_filt_table=texpr(bg_rs_filt,'all')
gene_names=unique(bg_filt_table[,9:10])
features=de_genes$id
mapped_gene_names=vector()
for (i in features) 
{  query=gene_names%>%filter(gene_id==i & gene_name != '.') ; n_hit=dim(query)[1]; if (n_hit==1) {mapped_gene_names=append(mapped_gene_names,query$gene_name[[1]]) } else
{mapped_gene_names=append(mapped_gene_names,'.') }    
}
#add the mapped gene names to the de genes table
de_genes$gene_name <- mapped_gene_names
de_genes <- de_genes[, c('feature','gene_name','id','fc','pval','qval')]



de_genes[,"log2fc"] <- log2(de_genes[,"fc"])
de_transcripts[,"log2fc"] <- log2(de_transcripts[,"fc"])


# Create results directory (and parent folders if missing)
dir.create("ballgown/de_results", recursive = TRUE, showWarnings = FALSE)

# Arrange results by p-value
de_transcripts <- arrange(de_transcripts, pval)
de_genes <- arrange(de_genes, pval)

# Write to CSV files inside the new folder
write.csv(de_transcripts,
          file = "ballgown/de_results/DE_transcripts_sorted.csv",
          row.names = FALSE)

write.csv(de_genes,
          file = "ballgown/de_results/DE_genes_sorted.csv",
          row.names = FALSE)


# Create results folder if it doesn't already exist
dir.create("results/DE_subsets", recursive = TRUE, showWarnings = FALSE)

# Subset significant transcripts and genes
subset_transcripts <- subset(de_transcripts, qval < 0.05)
subset_genes <- subset(de_genes, qval < 0.05)

# Save to CSV
write.csv(subset_transcripts,
          file = "results/DE_subsets/significant_transcripts_qval_0.05.csv",
          row.names = FALSE)

write.csv(subset_genes,
          file = "results/DE_subsets/significant_genes_qval_0.05.csv",
          row.names = FALSE)

cat("✅ Significant DE results saved to 'results/DE_subsets/'\n")
cat(paste("  - Transcripts:", nrow(subset_transcripts), "\n"))
cat(paste("  - Genes:", nrow(subset_genes), "\n"))


dir.create('ballgown/plots')

print('generating plots')


#gene expression for a isoforms of gene XIST
png('ballgown/plots/CYP6P5.png')
# Open PNG output
png("ballgown/plots/CYP6P5.png", width = 1800, height = 1000, res = 150)

# Draw transcript structure directly
plotTranscripts(
  gene = ballgown::geneIDs(bg_rs)[ballgown::geneNames(bg_rs) == "CYP6P5"],
  gown = bg_rs,
  main = "Gene CYP6P5 in sample 10-GGAFUN_DELTA_1",
  sample = "10-GGAFUN_DELTA_1_RNA_CACGCAAT-ATGGAAGG_L001"
)

# Close PNG device (very important!)
dev.off()

cat("✅ Plot saved to ballgown/plots/CYP6P5.png\n")
print(myplot)
dev.off()
#DONE



#volcano plot


#https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html

#mycolors <- c("red", "green", "black")
#names(mycolors) <- c("DOWN", "UP", "NO")
#p3 <- p2 + scale_colour_manual(values = mycolors)

#de_genes$diffexpressed[de_genes$log2fc < -0.6 & de_genes$qval < 0.05] <- "DOWN"
# Assign differential expression categories
de_genes$diffexpressed <- "NO"
de_genes$diffexpressed[de_genes$log2fc > 1 & de_genes$pval < 0.05] <- "UP"
de_genes$diffexpressed[de_genes$log2fc < -1 & de_genes$pval < 0.05] <- "DOWN"

# Add gene name labels for significant genes
de_genes$delabel <- NA
de_genes$delabel[de_genes$diffexpressed != "NO"] <- de_genes$gene_name[de_genes$diffexpressed != "NO"]

# Allow all labels to be plotted
options(ggrepel.max.overlaps = Inf)

# Create output folder if not exists
dir.create("ballgown/plots", showWarnings = FALSE)

# Volcano plot
png('ballgown/plots/volcano.png', width = 1800, height = 1000)

volcano <- ggplot(data = de_genes, aes(x = log2fc, y = -log10(pval),
                                       color = diffexpressed, label = delabel)) +
  geom_point(alpha = 0.8) +
  # geom_text_repel() +
  theme_minimal() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-1, 1), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dotted") +
  labs(title = "Volcano Plot of Differentially Expressed Genes",
       x = "log2 Fold Change",
       y = "-log10(p-value)",
       color = "Expression") +
  theme(text = element_text(size = 18))

print(volcano)
dev.off()

cat("✅ Volcano plot saved to: plots/volcano.png\n")




# Ensure results directory exists
dir.create("ballgown/plot2", showWarnings = FALSE, recursive = TRUE)

# Define thresholds
logFC_cutoff <- 1
pval_cutoff <- 0.05

# Categorize genes
de_genes <- de_genes %>%
  mutate(
    diffexpressed = case_when(
      log2fc > logFC_cutoff & pval < pval_cutoff  ~ "Upregulated",
      log2fc < -logFC_cutoff & pval < pval_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Select top genes for labeling
top_genes <- de_genes %>%
  filter(diffexpressed != "Not significant") %>%
  arrange(pval) %>%
  slice_head(n = 10)  # label top 10 most significant genes

# Volcano plot
volcano <- ggplot(de_genes, aes(x = log2fc, y = -log10(pval))) +
  geom_point(aes(color = diffexpressed), alpha = 0.8, size = 2.2) +
  scale_color_manual(
    values = c("Downregulated" = "#0072B2",   # blue
               "Not significant" = "gray70",
               "Upregulated" = "#D55E00")     # orange/red
  ) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff),
             linetype = "dashed", color = "black", linewidth = 0.7) +
  geom_hline(yintercept = -log10(pval_cutoff),
             linetype = "dotted", color = "black", linewidth = 0.7) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "gray50",
    segment.size = 0.3,
    min.segment.length = 0,
    fontface = "italic"
  ) +
  labs(
    title = "Differential Expression Volcano Plot",
    subtitle = paste0("Up/Down genes at |log2FC| >", logFC_cutoff,
                      " & p < ", pval_cutoff),
    x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~pvalue),
    color = "Expression"
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30"),
    axis.title = element_text(face = "bold"),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    legend.background = element_blank()
  )

# Save high-resolution publication-quality plot
ggsave("ballgown/plot2/Volcano_DE_genes_publication.png",
       plot = volcano, width = 10, height = 8, dpi = 400)

cat("✅ Volcano plot saved: plots/Volcano_DE_genes_publication.png\n")


# Example thresholds
logFC_cutoff <- 1
pval_cutoff <- 0.05

# Categorize genes
de_genes <- de_genes %>%
  mutate(
    diffexpressed = case_when(
      log2fc > logFC_cutoff & pval < pval_cutoff  ~ "Upregulated",
      log2fc < -logFC_cutoff & pval < pval_cutoff ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Select top genes to label (most significant ones)
top_genes <- de_genes %>%
  filter(pval < 0.05 & abs(log2fc) > 1) %>%
  arrange(pval) %>%
  slice_head(n = 10)  # label top 10

# Volcano plot
volcano <- ggplot(de_genes, aes(x = log2fc, y = -log10(pval), color = diffexpressed)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Downregulated" = "blue",
                                "Not significant" = "black",
                                "Upregulated" = "red")) +
  geom_vline(xintercept = c(-1, 1), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dotted") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),       # Use gene_name instead of delabel if available
    size = 4,
    box.padding = 0.3,
    point.padding = 0.1,
    segment.size = 0.2,
    segment.color = "gray60",
    min.segment.length = 0.2,
    max.overlaps = Inf,
    force = 1,                    # reduce repelling strength
    nudge_y = 0.5
  ) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "log2 Fold Change",
    y = "-log10(p-value)",
    color = "Expression"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

# Save
ggsave("ballgown/plots/Volcano_DE_genes_clean.png", plot = volcano, width = 10, height = 8, dpi = 300)

#MAPLOT
#https://davetang.org/muse/2017/10/25/getting-started-hisat-stringtie-ballgown/

png('plots/maplot.png',width = 1800, height = 1000)
de_transcripts$mean <- rowMeans(texpr(bg_chrX_filt))
maplot=ggplot(de_transcripts, aes(log2(mean), log2(fc), colour = qval<0.05)) +
  scale_color_manual(values=c("#999999", "#FF0000")) +
  geom_point() +
  theme(legend.text=element_text(size=20),legend.title=element_text(size=20)) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20)) +
  geom_hline(yintercept=0)
  



print(maplot)
dev.off()
#DONE


#heatmap
#Visualise the differentially expressed genes as a heatmap

#https://pydupont.github.io/BallGownTuto/Protocol.html



#get fpkm values 

png('ballgown/plots/heatmap_clustered.png')

fpkm = gexpr(bg_rs_filt)
fpkm = log2(fpkm+1)
g_ids=subset_genes$id
hits= which (de_genes$id %in% g_ids)
hit_frame=fpkm[hits,]
row.names(hit_frame) <- g_ids
heatmap_image=pheatmap(hit_frame)
print(heatmap_image)
dev.off()

png('ballgown/plots/heatmap_unclustered.png')
heatmap_image=pheatmap(hit_frame,cluster_rows = F,cluster_cols=F)
print(heatmap_image)
dev.off()


print("plots have been generated and saved in 'plots directory' ")



# Extract and log-transform FPKM
fpkm <- log2(gexpr(bg_rs_filt) + 1)

# Use only significant genes
g_ids <- subset_genes$id
hits <- which(de_genes$id %in% g_ids)
hit_frame <- fpkm[hits, ]
row.names(hit_frame) <- g_ids

# --- 📘 Shorten Sample Names ---
# Example pattern removal: remove long Illumina run suffixes
short_names <- colnames(hit_frame) %>%
  str_replace("_RNA.*", "") %>%          # removes everything after "_RNA"
  str_replace("-L001.*", "") %>%         # removes lane details
  str_replace("_.*", "")                 # optional: keep first section
colnames(hit_frame) <- short_names

# --- 🎨 Publication-Grade Aesthetics ---
png("ballgown/plots/heatmap_clustered.png", width = 2200, height = 1600, res = 300)

pheatmap(
  hit_frame,
  scale = "row",                        # Z-score normalization per gene
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  color = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(200), # blue-white-red gradient
  show_rownames = FALSE,                # hide gene IDs for cleaner look
  show_colnames = TRUE,
  fontsize_col = 12,
  fontsize = 10,
  main = "Clustered Heatmap of Differentially Expressed Genes",
  border_color = NA
)

dev.off()

cat("✅ High-resolution clustered heatmap saved to 'ballgown/plots/heatmap_clustered.png'\n")
