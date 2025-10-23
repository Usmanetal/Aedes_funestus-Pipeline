#!/bin/bash

# List of sample URLs separated by semicolon
URLS="https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_8-FANG_RNA_6/;https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_9-FG7_RNA/;https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_10-GGAFUN_DELTA_1_RNA/;https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_11-GGAFUN_DELTA_2_RNA/;https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_12-GGAFUN_DELTA_3_RNA/;https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_13-GGAFUN_UNX1_RNA/;https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_14-GGAFUN_UNX2_RNA/;https://cgr.liv.ac.uk/illum/LIMS18211_0f587c04e76299c6/Trimmed/Sample_15-GGAFUN_UNX3_RNA/"

# Create download directory if it doesn't exist
DOWNLOAD_DIR="downloads"
mkdir -p "$DOWNLOAD_DIR"

# Split URLs by semicolon and loop through each
IFS=';' read -ra ADDR <<< "$URLS"
for url in "${ADDR[@]}"; do
    echo "Processing $url"
    
    # Fetch the listing page HTML for this sample folder
    files=$(curl -s "$url" | grep -oP 'href="[^"]+\.fastq\.gz"' | sed 's/href="//;s/"$//')
    
    # Download each fastq.gz file found
    for file in $files; do
        file_url="${url}${file}"
        echo "Downloading $file_url"
        wget -c -P "$DOWNLOAD_DIR" "$file_url"
    done
done

echo "All downloads completed."
