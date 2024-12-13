# Name: go_heatmap_w_ratio.py
# Author: EY
# Date: December 2024
# Version: Python 3.9
# Description: Will plot heatmap of go terms with the ratio of up/down

import pandas as pd
import functools as ft
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


def go_gene_assoc(go_result, go_gene_mapping, deseq_results):
    """Will read in the .csv GO result from GOseq, the .txt file of gene ID and go term mapping, and the deseq
    result output from DEseq2. The significant GO term results will then be mapped to the gene IDs and the
    direction of expression will be determined by the log2foldchange in the deseq result output. The output is a nested
    dictionary - result_dict = {'go_term1':{'genes': ['geneID','geneID', 'geneID'], 'up_down':['up','down','up']},
    'go_term2':{'genes': ['geneID','geneID', 'geneID'], 'up_down':['up','down','down']}}"""
    result_dict = {}
    for category in go_result['category']:
        # find genes associated with each GO term
        associated_genes = []
        for gene in go_gene_mapping.loc[go_gene_mapping['go_terms'].apply(lambda terms: category in terms), 'gene']:
            # ensure gene is in deseq_results and if the pvalue is <= 0.05
            if gene in deseq_results['Unnamed: 0'].values:
                gene_row = deseq_results.loc[deseq_results['Unnamed: 0'] == gene]
                if gene_row['padj'].values[0] <= 0.05:
                    associated_genes.append(gene)
        # create result dictionary
        result_dict[category] = {'genes': associated_genes, 'up_down': []}
        for gene in associated_genes:
            if gene in deseq_results['Unnamed: 0'].values:
                # find the log2foldchange value
                log2fc = deseq_results.loc[deseq_results['Unnamed: 0'] == gene, 'log2FoldChange'].values[0]
                # append "up" or "down" to the 'up_down' list depending on the direction of expression
                result_dict[category]['up_down'].append('up' if log2fc > 0 else 'down')
    return(result_dict)


def calculate_up_fraction(result_dict):
    """will take the resulting dictionary from the function go_gene_assoc and from the up_down will calculate the ratio
    of genes with increasing expression (genes up/total genes) for each GO term. The result is a dictionary -
    result = {'go_term': ratio, 'go_term2': ratio}"""
    fraction_dict = {}
    # Loop through the result_dict
    for category, data in result_dict.items():
        # Count occurrences of "up" and "down"
        up_count = data['up_down'].count('up')
        down_count = data['up_down'].count('down')
        total_count = up_count + down_count
        # Calculate the fraction of "up"
        if total_count == 0:
            fraction = 0.0  # avoid division by zero if no "up" or "down" values exist
        else:
            fraction = up_count / total_count

        fraction_dict[category] = fraction

    return fraction_dict


def prepare_go_heatmap_ratio(treatment1, treatment2, treatment3, name1, name2, name3, ratio_dict1, ratio_dict2,
                             ratio_dict3):
    """Input is the result .csv from GOseq for the treatments of interest, the names of the treatments (will appear on
    the plot) and the ratio dictionaries from calculate_up_fraction. For all treatments, the intersection and differences
    between overlapping GO terms will be found and the result is a dataframe with one column for each treatment, the index
    is the GO terms, and the values within are the ratio of expression (this will be the value color-scaled in the heat
    map).

    similar to prepare_go_heatmap in go_heatmap.py, but it adds the direction of expression into the output"""

    # create a ratio column of the dataframe (combining input dictionary and df)
    treatment1["ratio"] = treatment1["category"].map(ratio_dict1)
    treatment2["ratio"] = treatment2["category"].map(ratio_dict2)
    treatment3["ratio"] = treatment3["category"].map(ratio_dict3)

    # rename the category column to contain the GO description as well
    treatment1["category"] = treatment1['category'].astype(str) + ":" + treatment1["term"]
    treatment2["category"] = treatment2['category'].astype(str) + ":" + treatment2["term"]
    treatment3["category"] = treatment3['category'].astype(str) + ":" + treatment3["term"]

    # subset the data to include only 'category' and 'ratio' columns
    treatment1_subset = treatment1[["category", "ratio"]]
    treatment2_subset = treatment2[["category", "ratio"]]
    treatment3_subset = treatment3[["category", "ratio"]]

    # rename the columns (to match what we want in the plot)
    treatment1_subset.columns = ['category', name1]
    treatment2_subset.columns = ['category', name2]
    treatment3_subset.columns = ['category', name3]

    # merge the dataframes
    dfs = [treatment1_subset, treatment2_subset, treatment3_subset]
    df_final = ft.reduce(lambda left, right: pd.merge(left, right, on='category', how='outer'), dfs)
    df_final = df_final.set_index(list(df_final)[0])

    # Finding intersections and differences for each treatment category (we will use this to reindex/reorder our dataframe

    # in all 3 treatments
    treatment1_treatment2 = set(treatment1_subset['category'].values).intersection(
        set(treatment2_subset['category'].values))
    treatment1_treatment2_treatment3 = treatment1_treatment2.intersection(set(treatment3_subset['category'].values))

    # treatment 1 and 3 only
    treatment1_treatment3 = set(treatment1_subset['category'].values).intersection(
        set(treatment3_subset['category'].values))
    treatment1_treatment3_no_treatment2 = treatment1_treatment3.difference(set(treatment2_subset['category'].values))

    # treatment 2 and 3 only
    treatment2_treatment3 = set(treatment2_subset['category'].values).intersection(
        set(treatment3_subset['category'].values))
    treatment2_treatment3_no_treatment1 = treatment2_treatment3.difference(set(treatment1_subset['category'].values))

    # treatment 1 only
    treatment1_treatment3_dif = set(treatment1_subset['category'].values).difference(
        set(treatment3_subset['category'].values))
    treatment1_only = treatment1_treatment3_dif.difference(set(treatment2_subset['category'].values))

    # treatment 2 only
    treatment2_treatment3_dif = set(treatment2_subset['category'].values).difference(
        set(treatment3_subset['category'].values))
    treatment2_only = treatment2_treatment3_dif.difference(set(treatment1_subset['category'].values))

    # treatment 3 only
    treatment3_treatment1_dif = set(treatment3_subset['category'].values).difference(
        set(treatment1_subset['category'].values))
    treatment3_only = treatment3_treatment1_dif.difference(set(treatment2_subset['category'].values))

    # reordering the categories to maintain specific subsets
    index_order = [
        list(treatment1_treatment2_treatment3), list(treatment1_treatment3_no_treatment2),
        list(treatment2_treatment3_no_treatment1), list(treatment1_only),
        list(treatment2_only), list(treatment3_only)
    ]

    # flattening the list of indices (needs to be flat to reindex)
    index_order_flat = [element for innerList in index_order for element in innerList]
    df_final = df_final.reindex(index_order_flat)
    # set order of final df
    df_final = df_final[[name1, name2, name3]]

    return df_final


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

# get the gene IDs and whether increasing/decreasing expression for each GO term in each comparison of interest
result_10v20 = go_gene_assoc(go_10v20, go_gene, deseq_10v20)
result_20v30 = go_gene_assoc(go_20v30, go_gene, deseq_20v30)
result_30v35 = go_gene_assoc(go_30v35, go_gene, deseq_30v35)

# get the ratio of increasing expression for each GO term in each comparison of interest
ratio_10v20 = calculate_up_fraction(result_10v20)
ratio_20v30 = calculate_up_fraction(result_20v30)
ratio_30v35 = calculate_up_fraction(result_30v35)

# create a dataframe that contains info at the intersection of all the GO terms between treatments
all_go = prepare_go_heatmap_ratio(go_10v20,go_20v30,go_30v35, "10v20", "20v30", "30v35", ratio_10v20, ratio_20v30,ratio_30v35)
# plot the heat map
formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))
labels = list(all_go.index.values)
plt.figure(figsize=(12, 8))
go_terms = sns.heatmap(all_go, cbar_kws={'format': formatter}, yticklabels=labels)
go_terms.figure.savefig(wd + 'go_heatmap_ratio.png', dpi='figure', format="png", bbox_inches='tight')
plt.close()

