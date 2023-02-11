#!/usr/bin/env python
import pandas as pd
import sys

unfilteredDataFile =  sys.argv[1]
filteredDataFile =  sys.argv[2]


unfilteredData = pd.read_table(unfilteredDataFile)

for i in unfilteredData.columns[11:]:
    current_sample = (unfilteredData.loc[:, ['gene_ID', 'transcript_ID', 'annot_gene_id', 'annot_transcript_id', 'annot_gene_name',
                       'annot_transcript_name', 'n_exons', 'length', 'gene_novelty', 'transcript_novelty', 'ISM_subtype',
                       i]])
    current_sample.to_csv(i+"_unfiltered_sample.tsv")


filteredData = pd.read_table(filteredDataFile)
for i in filteredData.columns[11:]:
    current_sample = (filteredData.loc[:, ['gene_ID', 'transcript_ID', 'annot_gene_id', 'annot_transcript_id', 'annot_gene_name',
                       'annot_transcript_name', 'n_exons', 'length', 'gene_novelty', 'transcript_novelty', 'ISM_subtype',
                       i]])
    current_sample.to_csv(i+"_filtered_sample.tsv")
