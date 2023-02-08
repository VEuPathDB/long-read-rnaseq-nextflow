#!/usr/bin/env python
import pandas as pd

data = pd.read_table("/Users/saikouybah/Documents/Long_Read_Workflow/Results/counts/results_talon_abundance_filtered.tsv")

for i in data.columns[11:]:
    current_sample = (data.loc[:, ['gene_ID', 'transcript_ID', 'annot_gene_id', 'annot_transcript_id', 'annot_gene_name',
                       'annot_transcript_name', 'n_exons', 'length', 'gene_novelty', 'transcript_novelty', 'ISM_subtype',
                       i]])
    current_sample.to_csv(i+"sample.tsv")
