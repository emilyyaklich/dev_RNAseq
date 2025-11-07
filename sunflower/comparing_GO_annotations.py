# Name: comparing_go_annotations
# Author: EY
# Date: October 2025
# Version: Python 3.9
# Description: Will compare annotations of entap vs InterProScan vs gogetter

import pandas as pd
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools import obo_parser
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from goatools.go_enrichment import GOEnrichmentStudy



def parse_go_terms_from_file(file_path):
    """
    parses an entap TSV file and processes the GO terms for each gene, grouping them by gene ID.

    Args:
    file_path (str): Path to the TSV file to parse.

    Returns:
    dict: A dictionary where keys are gene IDs and values are dictionaries of GO terms for each category.
    """

    df = pd.read_csv(file_path, sep='\t', low_memory=False)


    gene_dict = {}

    columns_of_interest = [
        'Database_EggNOG_GO_Biological',
        'Database_EggNOG_GO_Cellular',
        'Database_EggNOG_GO_Molecular',
        'Database_UniProt_GO_Biological',
        'Database_UniProt_GO_Cellular',
        'Database_UniProt_GO_Molecular',
        'Database_InterProScan_GO_Biological',
        'Database_InterProScan_GO_Cellular',
        'Database_InterProScan_GO_Molecular'
    ]

    # Iterate through each row
    for index, row in df.iterrows():
        # extract transcript ID
        query_sequence = row['Query_Sequence']

        # get the gene ID by removing the transcript suffix (.t1, .t2)
        gene_id = query_sequence.rsplit('.', 1)[0]

        for col in columns_of_interest:
            go_terms = row[col]
            if pd.notna(go_terms):
                # initialize a dictionary to store GO terms for this gene if not already present
                if gene_id not in gene_dict:
                    gene_dict[gene_id] = {
                        'Database_EggNOG_GO_Biological': set(),
                        'Database_EggNOG_GO_Cellular': set(),
                        'Database_EggNOG_GO_Molecular': set(),
                        'Database_UniProt_GO_Biological': set(),
                        'Database_UniProt_GO_Cellular': set(),
                        'Database_UniProt_GO_Molecular': set(),
                        'Database_InterProScan_GO_Biological': set(),
                        'Database_InterProScan_GO_Cellular': set(),
                        'Database_InterProScan_GO_Molecular': set()
                    }
                # split the GO terms string by commas and remove the trailing empty element
                go_terms_list = go_terms.strip(',').split(',')
                # add the GO terms to the set for the appropriate column
                gene_dict[gene_id][col].update(go_terms_list)

    # convert sets back to lists
    for gene_id, go_terms_dict in gene_dict.items():
        for col in go_terms_dict:
            gene_dict[gene_id][col] = list(go_terms_dict[col])

    return gene_dict

# function to extract the number of GO terms per gene
def get_go_terms_per_gene(go_dict):
    return [len(go_terms) for go_terms in go_dict.values()]


# function to calculate the number of unique (non-redundant) GO terms per gene for Entap
def get_entap_go_terms_per_gene_unique(entap_dict):
    go_terms_per_gene = []
    for gene, categories in entap_dict.items():
        unique_go_terms = set()

        for go_terms in categories.values():
            unique_go_terms.update(go_terms)

        go_terms_per_gene.append(len(unique_go_terms))

    return go_terms_per_gene


# get the genes with the highest number of GO terms for Entap v1 and v2
def get_high_go_term_genes(entap_dict, threshold):
    high_go_genes = {}

    for gene, categories in entap_dict.items():
        unique_go_terms = set()
        for go_terms in categories.values():
            unique_go_terms.update(go_terms)

        if len(unique_go_terms) > threshold:
            high_go_genes[gene] = unique_go_terms

    return high_go_genes


def get_go_term_info(go_terms, go):
    go_term_info = defaultdict(dict)
    for go_term in go_terms:
        if go_term in go:
            term = go[go_term]
            go_term_info[go_term]['name'] = term.name  # Guaranteed to be present
            go_term_info[go_term]['definition'] = getattr(term, 'def', None)  # 'def' is used for definitions
        else:
            go_term_info[go_term]['name'] = None
            go_term_info[go_term]['definition'] = None
    return go_term_info



# load interpro annotations
go_obo = '/home/emilyyaklich/Documents/dev/go-basic_09_18_2025.obo'
go = obo_parser.GODag(go_obo)
interpo_dict = defaultdict(set)

# read in the interpro results from my run
with open('/home/emilyyaklich/Documents/dev/Ha412GO_terms_interpro_09_2025.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue  # skip empty lines
        transcript, go_str = line.split('\t')
        gene = transcript.split('.')[0]  # remove .t1/.t2/etc.
        go_terms = go_str.split(',')
        interpo_dict[gene].update(go_terms)

interpo = {gene: sorted(terms) for gene, terms in interpo_dict.items()}


# read in the gogetter data (uses GOSlim terms)
file_path = '/home/emilyyaklich/Documents/dev/ha412_transcript_b3.fasta.blast.besthits.tsv.genes_to_GOSlim_terms_with_TAIR.tsv'
df = pd.read_csv(file_path, sep='\t')

gogetter_dict = {}
for index, row in df.iterrows():
    # extract gene and get rid of transcript ID
    gene_id = row['YourGene'].split('.')[0]

    go_term = row['GOID']

    # if the gene_id is not already in the dictionary add it with an empty set
    if gene_id not in gogetter_dict:
        gogetter_dict[gene_id] = set()

    # add the GO term to the set
    gogetter_dict[gene_id].add(go_term)

# convert the sets back to lists
for gene_id in gogetter_dict:
    gogetter_dict[gene_id] = list(gogetter_dict[gene_id])



# read in the files
entap_v1 = parse_go_terms_from_file('/home/emilyyaklich/Documents/dev/annotated_entap_EB.tsv')
entap_v2 = parse_go_terms_from_file('/home/emilyyaklich/Documents/dev/annotated_entap_full_db_EB.tsv')



genes_in_entap_v1 = set(entap_v1.keys())
genes_in_entap_v2 = set(entap_v2.keys())

# genes that are in entap_v1 but not in entap_v2
genes_in_entap_v1_not_in_entap_v2 = genes_in_entap_v1 - genes_in_entap_v2

# genes that are in entap_v2 but not in entap_v1
genes_in_entap_v2_not_in_entap_v1 = genes_in_entap_v2 - genes_in_entap_v1

# genes that are in both v1 and v2
genes_in_both_entap = genes_in_entap_v1 & genes_in_entap_v2




# see how many genes were annotated by interproscan in v1
interpro_gene_count_v1 = 0
for gene_id, go_terms_dict in entap_v1.items():
    if (go_terms_dict['Database_InterProScan_GO_Biological'] or
        go_terms_dict['Database_InterProScan_GO_Cellular'] or
        go_terms_dict['Database_InterProScan_GO_Molecular']):
        interpro_gene_count_v1 += 1
# see how many genes were annotated by interproscan in v2
interpro_gene_count_v2 = 0
for gene_id, go_terms_dict in entap_v2.items():
    if (go_terms_dict['Database_InterProScan_GO_Biological'] or
        go_terms_dict['Database_InterProScan_GO_Cellular'] or
        go_terms_dict['Database_InterProScan_GO_Molecular']):
        interpro_gene_count_v2 += 1

interpro_gene_count = len(interpo)

# For Gogetter
gogetter_gene_count = len(gogetter_dict)

# For Entap v1
entap_v1_gene_count = len([gene for gene, go_terms in entap_v1.items() if any(go_terms.values())])

# For Entap v2
entap_v2_gene_count = len([gene for gene, go_terms in entap_v2.items() if any(go_terms.values())])

# Create the bar chart
methods = ['InterProScan', 'Gogetter', 'Entap v1', 'Entap v2']
gene_counts = [interpro_gene_count, gogetter_gene_count, entap_v1_gene_count, entap_v2_gene_count]

plt.figure(figsize=(8, 5))
plt.bar(methods, gene_counts, color=['skyblue', 'lightgreen', 'salmon', 'lightcoral'])
plt.xlabel('Method', fontsize=12)
plt.ylabel('Number of Genes Annotated with GO Terms', fontsize=12)
plt.title('Comparison of GO Term Annotations Across Methods', fontsize=14)
plt.tight_layout()
plt.savefig('number_of_genes_w_go_terms.png', dpi=300)
plt.show()




# prepare the data for plotting

# get the number of GO terms per gene for each method
interpro_go_terms_per_gene = get_go_terms_per_gene(interpo)
gogetter_go_terms_per_gene = get_go_terms_per_gene(gogetter_dict)

# get the number of unique GO terms per gene for Entap v1 and Entap v2 (non-redundant)
entap_v1_go_terms_per_gene_unique = get_entap_go_terms_per_gene_unique(entap_v1)
entap_v2_go_terms_per_gene_unique = get_entap_go_terms_per_gene_unique(entap_v2)

# combine the data into a single df for easier plotting
plot_data_unique = pd.DataFrame({
    'Method': ['InterProScan'] * len(interpro_go_terms_per_gene) +
              ['Gogetter'] * len(gogetter_go_terms_per_gene) +
              ['Entap v1'] * len(entap_v1_go_terms_per_gene_unique) +
              ['Entap v2'] * len(entap_v2_go_terms_per_gene_unique),
    'GO_Terms_Per_Gene': interpro_go_terms_per_gene +
                         gogetter_go_terms_per_gene +
                         entap_v1_go_terms_per_gene_unique +
                         entap_v2_go_terms_per_gene_unique
})


plt.figure(figsize=(10, 6))
sns.boxplot(x='Method', y='GO_Terms_Per_Gene', data=plot_data_unique, palette="Set2")
plt.title('Distribution of Unique GO Terms per Gene', fontsize=14)
plt.xlabel('Method', fontsize=12)
plt.ylabel('Number of Unique GO Terms per Gene', fontsize=12)
plt.tight_layout()
plt.savefig('unique_go_terms_per_gene.png', dpi=300)
plt.show()




high_go_genes_v1 = get_high_go_term_genes(entap_v1, threshold=300)
high_go_genes_v2 = get_high_go_term_genes(entap_v2, threshold=100)

print("High GO Term Genes (Entap v1):")
for gene, go_terms in list(high_go_genes_v1.items()):
    print(f"Gene: {gene}, Number of GO Terms: {len(go_terms)}")

print("\nHigh GO Term Genes (Entap v2):")
for gene, go_terms in list(high_go_genes_v2.items())[:5]:
    print(f"Gene: {gene}, Number of GO Terms: {len(go_terms)}")



wus_dict_entap = {'g51546': []}
for item, value in entap_v1['g51546'].items():
    for thing in entap_v1['g51546'][item]:
        wus_dict_entap['g51546'].append(thing)

go_term_info = get_go_term_info(wus_dict_entap['g51546'], go)

for go_term, info in go_term_info.items():
    print(f"{go_term}:")
    print(f"  Name: {info['name']}")
    print(f"  Definition: {info['definition']}\n")



# compare entap v1 with GOgetter

genes_go_getter = set(gogetter_dict.keys())

# genes that are in entap_v1 but not in entap_v2
genes_in_entap_v1_not_in_gogetter = genes_in_entap_v1 - genes_go_getter

# genes that are in entap_v2 but not in entap_v1
genes_in_gogetter_not_in_entap_v1 = genes_go_getter - genes_in_entap_v1

# genes that are in both v1 and v2
genes_in_both_entap_go_getter = genes_in_entap_v1 & genes_go_getter

only_entap_v1 = len(genes_in_entap_v1_not_in_gogetter)
only_gogetter = len(genes_in_gogetter_not_in_entap_v1)
print(genes_in_gogetter_not_in_entap_v1)
both = len(genes_in_both_entap_go_getter)

categories = ['Only Entap_v1', 'Only GoGetter', 'Both']
counts = [only_entap_v1, only_gogetter, both]

fig, ax = plt.subplots(figsize=(8, 6))

bar_width = 0.6

ax.barh(categories, [only_entap_v1, only_gogetter, both], color=['#1f77b4', '#ff7f0e', '#2ca02c'])

ax.set_xlabel('Number of Genes')
ax.set_title('Comparison of Genes in Entap_v1 vs. GoGetter')

for i, v in enumerate(counts):
    ax.text(v + 50, i, str(v), ha='center', va='center', color='black', fontweight='bold')

plt.tight_layout()
plt.savefig('go_getter_vs_entap_v1.png', dpi=300)

plt.show()

filtered_go_getter = {gene: go_terms for gene, go_terms in gogetter_dict.items() if gene in genes_in_gogetter_not_in_entap_v1}

all_go_terms = set()
for go_terms in gogetter_dict.values():
    all_go_terms.update(go_terms)


filtered_go_terms = set()
for go_terms in filtered_go_getter.values():
    filtered_go_terms.update(go_terms)


unique_filtered_terms = filtered_go_terms - all_go_terms

common_terms = set(filtered_go_terms & all_go_terms)

print(f"Unique GO terms in the 12,000 genes: {unique_filtered_terms}")
print(f"Common GO terms: {common_terms}")


gene2gos = {}
for gene, go_terms_list in gogetter_dict.items():
    gene2gos[gene] = set(go_terms_list)

background_genes = list(gogetter_dict.keys())
study_genes = list(genes_in_gogetter_not_in_entap_v1)



study_gene2gos = {gene: set(gogetter_dict[gene]) for gene in study_genes if gene in gogetter_dict}
background_gene2gos = {gene:  set(gogetter_dict[gene]) for gene in background_genes if gene in gogetter_dict}

goeaobj = GOEnrichmentStudy(
    background_genes,             # Background gene set (all genes in your reference)
    background_gene2gos,          # Background gene-to-GO mappings (dictionary)
    go,                        # GO DAG (Ontology)
    methods=['fdr_bh'],  # Multiple comparison correction methods
)

results = goeaobj.run_study(study_genes)

goea_results_sig = [r for r in results if r.p_fdr_bh < 0.05]

print("Enrichment Results:")
for res in results:
    print(f"GO Term: {res.GO} - p-value: {res.p_fdr_bh:.5f} - FDR Corrected p-value: {res.p_fdr_bh:.5f}")
    print(f"Description: {res.name}")
    print("----")
print("\n" + "="*50 + "\n")  # Divider between results


enriched_terms_with_defs = []
for res in results:
    if res.p_fdr_bh < 0.05 and res.enrichment == 'e':  # Example: assuming positive score means increasing
        enriched_terms_with_defs.append({
            'GO Term': res.GO,
            'Description': res.name,
            'Enrichment Score': res.enrichment  # Include enrichment score in the output
        })



for term in enriched_terms_with_defs:
    print(f"GO Term: {term['GO Term']}")
    print(f"Description: {term['Description']}")
    print("-" * 50)


