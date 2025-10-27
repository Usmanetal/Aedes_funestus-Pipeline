#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepDE.py3 - Generate gene/transcript count matrices from StringTie output
"""

import os
import sys
import glob
import pandas as pd

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Generate gene/transcript count matrices from StringTie outputs')
    parser.add_argument('-i', '--input', required=True, help='Text file with paths to sample GTFs, one per line')
    parser.add_argument('-o', '--output', default='.', help='Output directory (default=current folder)')
    args = parser.parse_args()

    with open(args.input) as f:
        samples = [line.strip() for line in f if line.strip()]

    gene_count = pd.DataFrame()
    transcript_count = pd.DataFrame()

    for gtf_file in samples:
        sample_name = os.path.basename(os.path.dirname(gtf_file))
        df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=[
            'chrom','source','feature','start','end','score','strand','frame','attribute'
        ])

        # Keep only transcripts or exons for counting
        df_transcripts = df[df['feature']=='transcript']
        df_exons = df[df['feature']=='exon']

        # Extract gene_id and transcript_id
        df_exons['gene_id'] = df_exons['attribute'].str.extract('gene_id "([^"]+)"')
        df_exons['transcript_id'] = df_exons['attribute'].str.extract('transcript_id "([^"]+)"')

        # Count exons per gene and transcript
        gene_counts = df_exons.groupby('gene_id').size().rename(sample_name)
        transcript_counts = df_exons.groupby('transcript_id').size().rename(sample_name)

        gene_count = pd.concat([gene_count, gene_counts], axis=1)
        transcript_count = pd.concat([transcript_count, transcript_counts], axis=1)

    # Fill NaNs with 0
    gene_count = gene_count.fillna(0).astype(int)
    transcript_count = transcript_count.fillna(0).astype(int)

    os.makedirs(args.output, exist_ok=True)
    gene_count.to_csv(os.path.join(args.output,'gene_count_matrix.csv'))
    transcript_count.to_csv(os.path.join(args.output,'transcript_count_matrix.csv'))

    print(f"✅ Gene and transcript count matrices generated in {args.output}")

if __name__ == '__main__':
    main()
