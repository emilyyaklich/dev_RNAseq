# Name: go_heatmap.py
# Author: EY
# Date: December 2024
# Version: Python 3.9
# Description: Will plot heatmap of GO Enrichment

import pandas as pd
import functools as ft
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


def prepare_go_heatmap(treatment1, treatment2, treatment3, name1, name2, name3):
    """Input is the result .csv from GOseq for the treatments of interest and the names of the treatments (will appear on
    the plot) For all treatments, the intersection and differences
    between overlapping GO terms will be found and the result is a dataframe with one column for each treatment, the index
    is the GO terms, and the values within are the adjusted p-values (this will be the value color-scaled in the heat
    map)."""
    treatment1["category"] = treatment1['category'].astype(str) + ":" + treatment1["term"]
    treatment2["category"] = treatment2['category'].astype(str) + ":" + treatment2["term"]
    treatment3["category"] = treatment3['category'].astype(str) + ":" + treatment3["term"]


    treatment1_subset = treatment1[["category", "over_represented_pvalue"]]
    treatment2_subset = treatment2[["category", "over_represented_pvalue"]]
    treatment3_subset = treatment3[["category", "over_represented_pvalue"]]

    treatment1_subset.columns = ['category', '10v20']
    treatment2_subset.columns = ['category', '20v30']
    treatment3_subset.columns = ['category', '30v35']

    dfs = [treatment1_subset, treatment2_subset, treatment3_subset]
    df_final = ft.reduce(lambda left, right: pd.merge(left, right, on='category', how='outer'), dfs)
    df_final = df_final.set_index(list(df_final)[0])

    # all 3
    treatment1_treatment2 = set(treatment1_subset['category'].values).intersection(set(treatment2_subset['category'].values))
    treatment1_treatment2_treatment3 = treatment1_treatment2.intersection(set(treatment3_subset['category'].values))

    # treatment 1 and 3
    treatment1_treatment3 = set(treatment1_subset['category'].values).intersection(set(treatment3_subset['category'].values))
    treatment1_treatment3_no_treatment2 = treatment1_treatment3.difference(set(treatment2_subset['category'].values))

    # treatment 2 and 3
    treatment2_treatment3 = set(treatment2_subset['category'].values).intersection(set(treatment3_subset['category'].values))
    treatment2_treatment3_no_treatment1 = treatment2_treatment3.difference(set(treatment1_subset['category'].values))

    # treatment 1 only

    treatment1_treatment3_dif = set(treatment1_subset['category'].values).difference(set(treatment3_subset['category'].values))
    treatment1_only = treatment1_treatment3_dif.difference(set(treatment2_subset['category'].values))

    # treatment 2 only
    treatment2_treatment3_dif = set(treatment2_subset['category'].values).difference(set(treatment3_subset['category'].values))
    treatment2_only = treatment2_treatment3_dif.difference(set(treatment1_subset['category'].values))
    # treatment 3 only
    treatment3_treatment1_dif = set(treatment3_subset['category'].values).difference(set(treatment1_subset['category'].values))
    treatment3_only = treatment3_treatment1_dif.difference(set(treatment2_subset['category'].values))

    index_order = [list(treatment1_treatment2_treatment3), list(treatment1_treatment3_no_treatment2), list(treatment2_treatment3_no_treatment1),
                   list(treatment1_only)
        , list(treatment2_only), list(treatment3_only)]

    index_order_flat = [element for innerList in index_order for element in innerList]
    df_final = df_final.reindex(index_order_flat)
    df_final = df_final[[name1, name2, name3]]


    return df_final

if __name__ == '__main__':
    wd = '/Users/emilyyaklich/Documents/deg/'
    # high salt
    _10v20 = pd.read_csv(wd + '10_v20.csv')
    # low nutrient
    _20v30 = pd.read_csv(wd + '20_v30.csv')
    # combo
    _30v35 = pd.read_csv(wd + '30_v35.csv')
    all_go = prepare_go_heatmap(_10v20,_20v30,_30v35, "10v20", "20v30", "30v35")
    formatter = tkr.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    labels = list(all_go.index.values)
    plt.figure(figsize=(12, 8))
    go_terms = sns.heatmap(all_go, cbar_kws={'format': formatter}, yticklabels=labels)
    go_terms.figure.savefig(wd + 'go_heatmap.png', dpi='figure', format="png", bbox_inches='tight')
    plt.close()

