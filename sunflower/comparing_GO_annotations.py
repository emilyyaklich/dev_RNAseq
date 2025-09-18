# Name: comparing_go_annotations
# Author: EY
# Date: September 2025
# Version: Python 3.9
# Description: Will compare annotations of entap vs InterProScan

import pandas as pd
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools import obo_parser
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np


go_obo = '/home/emilyyaklich/Documents/dev/go-basic_09_18_2025.obo'
go = obo_parser.GODag(go_obo)

# intialize dictionary
gene_go_dict = defaultdict(set)

with open('/home/emilyyaklich/Documents/dev/Ha412GO_terms_interpro_09_2025.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue  # skip empty lines
        transcript, go_str = line.split('\t')
        gene = transcript.split('.')[0]  # remove .t1/.t2/etc.
        go_terms = go_str.split(',')
        gene_go_dict[gene].update(go_terms)

gene_go_dict = {gene: sorted(terms) for gene, terms in gene_go_dict.items()}

gene_go_dict['g51546']


genes_with_no_go_txt = [
    gene for gene, go_terms in gene_go_dict.items()
    if not go_terms  # empty list evaluates to False
]



# Nested dict: gene -> { 'eggnog': set(), 'interpro': set() }
gene_go_dict_eb = defaultdict(lambda: {'eggnog': set(), 'interpro': set()})

# Define GO-related columns to combine
egg_cols = [
    'Database_EggNOG_GO_Biological',
    'Database_EggNOG_GO_Cellular',
    'Database_EggNOG_GO_Molecular'
]

interpro_cols = [
    'Database_InterProScan_GO_Biological',
    'Database_InterProScan_GO_Cellular',
    'Database_InterProScan_GO_Molecular'
]

# Open and parse the file
with open('/home/emilyyaklich/Documents/dev/annotated_entap_EB.tsv', 'r') as f:
    # Parse header to get column index mapping
    header = f.readline().strip().split('\t')
    col_idx = {col: i for i, col in enumerate(header)}

    for line in f:
        fields = line.strip().split('\t')
        if len(fields) < len(header):
            continue  # skip incomplete or malformed lines

        transcript = fields[col_idx['Query_Sequence']]
        gene = transcript.split('.')[0]  # strip .t1/.t2

        # Add GO terms from EggNOG
        for col in egg_cols:
            idx = col_idx.get(col)
            if idx is not None:
                raw = fields[idx].strip()
                if raw and raw != 'NaN':
                    terms = [term.strip() for term in raw.split(',') if term.strip().startswith('GO:')]
                    gene_go_dict_eb[gene]['eggnog'].update(terms)

        # Add GO terms from InterProScan
        for col in interpro_cols:
            idx = col_idx.get(col)
            if idx is not None:
                raw = fields[idx].strip()
                if raw and raw != 'NaN':
                    terms = [term.strip() for term in raw.split(',') if term.strip().startswith('GO:')]
                    gene_go_dict_eb[gene]['interpro'].update(terms)

# Convert sets to sorted lists
gene_go_eb = {
    gene: {
        'eggnog': sorted(go_terms['eggnog']),
        'interpro': sorted(go_terms['interpro'])
    }
    for gene, go_terms in gene_go_dict_eb.items()
}

len(gene_go_dict_eb['g51546']['eggnog'])
len(gene_go_dict_eb['g51546']['interpro'])

# plot GO terms per gene
go_counts_txt = [len(go_terms) for go_terms in gene_go_dict.values()]

plt.figure(figsize=(10,6))
plt.hist(go_counts_txt, bins=30, color='skyblue', edgecolor='black')
plt.title('Distribution of GO Terms per Gene InterProScan')
plt.xlabel('Number of GO Terms per Gene')
plt.grid(axis='y', alpha=0.75)
plt.show()


eggnog_counts = [len(terms['eggnog']) for terms in gene_go_dict_eb.values()]
interpro_counts = [len(terms['interpro']) for terms in gene_go_dict_eb.values()]
total_counts = [e + i for e, i in zip(eggnog_counts, interpro_counts)]


plt.figure(figsize=(12,7))
plt.hist([eggnog_counts, interpro_counts], bins=30, stacked=True,
         label=['EggNOG GO Terms', 'InterPro GO Terms'],
         color=['orange', 'green'], edgecolor='black')
plt.title('Stacked GO Terms per Gene enTap')
plt.xlabel('Number of GO Terms per Gene')
plt.legend()
plt.grid(axis='y', alpha=0.75)
plt.show()


genes_txt = set(gene_go_dict.keys())
genes_tsv = set(gene_go_dict_eb.keys())
only_in_txt = genes_txt - genes_tsv       # Genes in txt but NOT in tsv
only_in_tsv = genes_tsv - genes_txt       # Genes in tsv but NOT in txt
in_both = genes_txt & genes_tsv


print(f"Genes only in .txt dict: {len(only_in_txt)}")
print(f"Genes only in .tsv dict: {len(only_in_tsv)}")
print(f"Genes in both: {len(in_both)}")