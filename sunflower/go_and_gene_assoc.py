# Name: go_and_gene_assoc.py
# Author: EY
# Date: December 2024
# Version: Python 3.9
# Description: 

import pandas as pd

import functools as ft
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


def go_gene_assoc(go_result, go_gene_mapping, deseq_results):
    result_dict = {}
    for category in go_result['category']:
        # Find genes associated with each GO term
        associated_genes = []
        for gene in go_gene_mapping.loc[
            go_gene_mapping['go_terms'].apply(lambda terms: category in terms), 'gene']:
            # Check if the gene is in deseq_results and if the pvalue is <= 0.05
            if gene in deseq_results['Unnamed: 0'].values:
                gene_row = deseq_results.loc[deseq_results['Unnamed: 0'] == gene]
                if gene_row['padj'].values[0] <= 0.05:
                    associated_genes.append(gene)
        # Initialize the category in the results_dict
        result_dict[category] = {'genes': associated_genes, 'up_down': []}
        for gene in associated_genes:
            if gene in deseq_results['Unnamed: 0'].values:
                # Find the log2FoldChange value
                log2fc = deseq_results.loc[deseq_results['Unnamed: 0'] == gene, 'log2FoldChange'].values[0]
                # Append "up" or "down" to the 'up_down' list
                result_dict[category]['up_down'].append('up' if log2fc > 0 else 'down')
    return(result_dict)


def calculate_up_fraction(result_dict):
    # Initialize the new dictionary
    fraction_dict = {}

    # Loop through the result_dict
    for category, data in result_dict.items():
        # Count occurrences of "up" and "down"
        up_count = data['up_down'].count('up')
        down_count = data['up_down'].count('down')
        total_count = up_count + down_count

        # Calculate the fraction of "up"
        if total_count == 0:
            fraction = 0.0  # Avoid division by zero if no "up" or "down" values exist
        else:
            fraction = up_count / total_count

        # Add the fraction to the new dictionary
        fraction_dict[category] = fraction

    return fraction_dict


wd = '/Users/emilyyaklich/Documents/'
# read in GO_result dfs
go_10v20 = pd.read_csv(wd + 'deg/10_v20.csv')
go_20v30 = pd.read_csv(wd + 'deg/20_v30.csv')
go_30v35 = pd.read_csv(wd + 'deg/30_v35.csv')

# read in go-gene annotations
go_gene = pd.read_csv(wd + 'Ha412GO_terms_interpro_b3.txt', sep='\t', header=None, names=['gene', 'go_terms'])

# read in deseq results
deseq_10v20 = pd.read_csv(wd + 'pairwise/result_10D_v_20D.csv')
deseq_20v30 = pd.read_csv(wd + 'pairwise/result_20D_v_30D.csv')
deseq_30v35 = pd.read_csv(wd + 'pairwise/result_30D_v_35D.csv')

result_10v20 = go_gene_assoc(go_10v20, go_gene, deseq_10v20)
result_20v30 = go_gene_assoc(go_20v30, go_gene, deseq_20v30)
result_30v35 = go_gene_assoc(go_30v35, go_gene, deseq_30v35)

ratio_10v20 = calculate_up_fraction(result_10v20)
ratio_20v30 = calculate_up_fraction(result_20v30)
ratio_30v35 = calculate_up_fraction(result_30v35)



def prepare_go_heatmap(treatment1, treatment2, treatment3, name1, name2, name3, ratio_dict1, ratio_dict2, ratio_dict3):

    # Creating the ratio column based on the dictionaries
    treatment1["ratio"] = treatment1["category"].map(ratio_dict1)
    treatment2["ratio"] = treatment2["category"].map(ratio_dict2)
    treatment3["ratio"] = treatment3["category"].map(ratio_dict3)
    
    # Creating the 'category' column by concatenating 'category' and 'term'
    treatment1["category"] = treatment1['category'].astype(str) + ":" + treatment1["term"]
    treatment2["category"] = treatment2['category'].astype(str) + ":" + treatment2["term"]
    treatment3["category"] = treatment3['category'].astype(str) + ":" + treatment3["term"]


    # Subsetting the data to include only 'category' and 'ratio' columns
    treatment1_subset = treatment1[["category", "ratio"]]
    treatment2_subset = treatment2[["category", "ratio"]]
    treatment3_subset = treatment3[["category", "ratio"]]

    # Renaming columns for clarity
    treatment1_subset.columns = ['category', name1]
    treatment2_subset.columns = ['category', name2]
    treatment3_subset.columns = ['category', name3]

    # Merging the dataframes
    dfs = [treatment1_subset, treatment2_subset, treatment3_subset]
    df_final = ft.reduce(lambda left, right: pd.merge(left, right, on='category', how='outer'), dfs)
    df_final = df_final.set_index(list(df_final)[0])

    # Finding intersections and differences for each treatment category
    treatment1_treatment2 = set(treatment1_subset['category'].values).intersection(set(treatment2_subset['category'].values))
    treatment1_treatment2_treatment3 = treatment1_treatment2.intersection(set(treatment3_subset['category'].values))

    treatment1_treatment3 = set(treatment1_subset['category'].values).intersection(set(treatment3_subset['category'].values))
    treatment1_treatment3_no_treatment2 = treatment1_treatment3.difference(set(treatment2_subset['category'].values))

    treatment2_treatment3 = set(treatment2_subset['category'].values).intersection(set(treatment3_subset['category'].values))
    treatment2_treatment3_no_treatment1 = treatment2_treatment3.difference(set(treatment1_subset['category'].values))

    treatment1_treatment3_dif = set(treatment1_subset['category'].values).difference(set(treatment3_subset['category'].values))
    treatment1_only = treatment1_treatment3_dif.difference(set(treatment2_subset['category'].values))

    treatment2_treatment3_dif = set(treatment2_subset['category'].values).difference(set(treatment3_subset['category'].values))
    treatment2_only = treatment2_treatment3_dif.difference(set(treatment1_subset['category'].values))

    treatment3_treatment1_dif = set(treatment3_subset['category'].values).difference(set(treatment1_subset['category'].values))
    treatment3_only = treatment3_treatment1_dif.difference(set(treatment2_subset['category'].values))

    # Reordering the categories to maintain specific subsets
    index_order = [
        list(treatment1_treatment2_treatment3), list(treatment1_treatment3_no_treatment2),
        list(treatment2_treatment3_no_treatment1), list(treatment1_only),
        list(treatment2_only), list(treatment3_only)
    ]

    # Flattening the list of indices
    index_order_flat = [element for innerList in index_order for element in innerList]
    df_final = df_final.reindex(index_order_flat)

    # Selecting only the columns corresponding to the treatment ratios
    df_final = df_final[[name1, name2, name3]]

    return df_final

all_go = prepare_go_heatmap(go_10v20,go_20v30,go_30v35, "10v20", "20v30", "30v35", ratio_10v20, ratio_20v30,ratio_30v35)
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))
labels = list(all_go.index.values)
plt.figure(figsize=(12, 8))
go_terms = sns.heatmap(all_go, cbar_kws={'format': formatter}, yticklabels=labels)
go_terms.figure.savefig(wd + 'go_heatmap_ratio.png', dpi='figure', format="png", bbox_inches='tight')
plt.close()

