# Name: go_heatmap.py
# Author: EY
# Date: July 2023
# Version: Python 3.9
# Description: Will plot heatmap of GO Enrichment

import pandas as pd
import functools as ft
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr


def prepare_go_heatmap(treatment1, lownut, combo):
    treatment1["category"] = treatment1['category'].astype(str) + ":" + treatment1["term"]
    lownut["category"] = lownut['category'].astype(str) + ":" + lownut["term"]
    combo["category"] = combo['category'].astype(str) + ":" + combo["term"]


    highsalt_subset = treatment1[["category", "over_represented_pvalue"]]
    lownut_subset = lownut[["category", "over_represented_pvalue"]]
    combo_subset = combo[["category", "over_represented_pvalue"]]

    highsalt_subset.columns = ['category', 'highsalt']
    lownut_subset.columns = ['category', 'lownut']
    combo_subset.columns = ['category', 'combo']

    dfs = [highsalt_subset, lownut_subset, combo_subset]
    df_final = ft.reduce(lambda left, right: pd.merge(left, right, on='category', how='outer'), dfs)
    df_final = df_final.set_index(list(df_final)[0])

    # all 3
    highsalt_lownut = set(highsalt_subset['category'].values).intersection(set(lownut_subset['category'].values))
    highsalt_lownut_combo = highsalt_lownut.intersection(set(combo_subset['category'].values))

    # combo and high salt
    highsalt_combo = set(highsalt_subset['category'].values).intersection(set(combo_subset['category'].values))
    highsalt_combo_no_lownut = highsalt_combo.difference(set(lownut_subset['category'].values))

    # combo and low nut
    lownut_combo = set(lownut_subset['category'].values).intersection(set(combo_subset['category'].values))
    lownut_combo_no_highsalt = lownut_combo.difference(set(highsalt_subset['category'].values))

    # high salt only

    highsalt_combo_dif = set(highsalt_subset['category'].values).difference(set(combo_subset['category'].values))
    highsalt_only = highsalt_combo_dif.difference(set(lownut_subset['category'].values))

    lownut_combo_dif = set(lownut_subset['category'].values).difference(set(combo_subset['category'].values))
    lownut_only = lownut_combo_dif.difference(set(highsalt_subset['category'].values))

    combo_highsalt_dif = set(combo_subset['category'].values).difference(set(highsalt_subset['category'].values))
    combo_only = combo_highsalt_dif.difference(set(lownut_subset['category'].values))

    index_order = [list(highsalt_lownut_combo), list(highsalt_combo_no_lownut), list(lownut_combo_no_highsalt),
                   list(highsalt_only)
        , list(lownut_only), list(combo_only)]

    index_order_flat = [element for innerList in index_order for element in innerList]
    df_final = df_final.reindex(index_order_flat)
    df_final = df_final[['highsalt', 'combo', 'lownut']]


    return df_final

if __name__ == '__main__':
    wd = '/Users/emilyyaklich/Documents/deg/'
    # high salt
    highsalt = pd.read_csv(wd + '10_v20.csv')
    # low nutrient
    lownut = pd.read_csv(wd + '20_v30.csv')
    # combo
    combo = pd.read_csv(wd + '30_v35.csv')
    all_go = prepare_go_heatmap(highsalt,lownut,combo)
    formatter = tkr.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    labels = list(all_go.index.values)
    plt.figure(figsize=(12, 8))
    go_terms = sns.heatmap(all_go, cbar_kws={'format': formatter}, yticklabels=labels)
    go_terms.figure.savefig(wd + 'go_heatmap.png', dpi='figure', format="png", bbox_inches='tight')
    plt.close()


# BP
    # high salt
    highsalt_bp = pd.read_csv(wd + 'go_highsalt_bp.csv')
    # low nutrient
    lownut_bp = pd.read_csv(wd + 'go_lownut_bp.csv')
    # combo
    combo_bp = pd.read_csv(wd + 'go_combo_bp.csv')

    bp_go = prepare_go_heatmap(highsalt_bp, lownut_bp, combo_bp)
    labels = list(bp_go.index.values)
    go_terms_bp=sns.heatmap(bp_go,cbar_kws={'format':formatter},cmap="Greens",center=0.00001, yticklabels=labels)
    go_terms_bp.set_title("BP")
    go_terms_bp.figure.savefig(wd+'go_heatmap_bp.png', dpi='figure', format="png",bbox_inches='tight')
    plt.close()

# MF
    # high salt
    highsalt_mf = pd.read_csv(wd + 'go_highsalt_mf.csv')
    # low nutrient
    lownut_mf = pd.read_csv(wd + 'go_lownut_mf.csv')
    # combo
    combo_mf = pd.read_csv(wd + 'go_combo_mf.csv')

    mf_go = prepare_go_heatmap(highsalt_mf, lownut_mf, combo_mf)
    labels = list(mf_go.index.values)
    go_terms_mf=sns.heatmap(mf_go,cbar_kws={'format':formatter},cmap="Blues",center=0.00001, yticklabels=labels)
    go_terms_mf.set_title("MF")
    go_terms_mf.figure.savefig(wd + 'go_heatmap_mf.png', dpi='figure', format="png", bbox_inches='tight')
    plt.close()

# CC
    # high salt
    highsalt_cc = pd.read_csv(wd + 'go_highsalt_cc.csv')
    # low nutrient
    lownut_cc = pd.read_csv(wd + 'go_lownut_cc.csv')
    # combo
    combo_cc = pd.read_csv(wd + 'go_combo_cc.csv')

    cc_go = prepare_go_heatmap(highsalt_cc, lownut_cc, combo_cc)
    labels = list(cc_go.index.values)
    go_terms_cc = sns.heatmap(cc_go, cbar_kws={'format': formatter}, cmap="BuPu", center=0.00001, yticklabels=labels)
    go_terms_cc.set_title("CC")
    go_terms_cc.figure.savefig(wd + 'go_heatmap_cc.png', dpi='figure', format="png", bbox_inches='tight')
    plt.close()




highsalt_subset = highsalt[["category", "over_represented_pvalue"]]
lownut_subset = lownut[["category", "over_represented_pvalue"]]
combo_subset = combo[["category", "over_represented_pvalue"]]

highsalt_subset.columns = ['category','highsalt']
lownut_subset.columns = ['category','lownut']
combo_subset.columns = ['category','combo']

dfs=[highsalt_subset,lownut_subset,combo_subset]
df_final = ft.reduce(lambda left, right: pd.merge(left, right, on='category', how='outer'), dfs)
df_final=df_final.set_index(list(df_final)[0])


# all 3
highsalt_lownut=set(highsalt_subset['category'].values).intersection(set(lownut_subset['category'].values))
highsalt_lownut_combo=highsalt_lownut.intersection(set(combo_subset['category'].values))

# combo and high salt
highsalt_combo=set(highsalt_subset['category'].values).intersection(set(combo_subset['category'].values))
highsalt_combo_no_lownut=highsalt_combo.difference(set(lownut_subset['category'].values))

# combo and low nut
lownut_combo=set(lownut_subset['category'].values).intersection(set(combo_subset['category'].values))
lownut_combo_no_highsalt=lownut_combo.difference(set(highsalt_subset['category'].values))

# high salt only

highsalt_combo_dif=set(highsalt_subset['category'].values).difference(set(combo_subset['category'].values))
highsalt_only=highsalt_combo_dif.difference(set(lownut_subset['category'].values))

lownut_combo_dif=set(lownut_subset['category'].values).difference(set(combo_subset['category'].values))
lownut_only=lownut_combo_dif.difference(set(highsalt_subset['category'].values))

combo_highsalt_dif=set(combo_subset['category'].values).difference(set(highsalt_subset['category'].values))
combo_only=combo_highsalt_dif.difference(set(lownut_subset['category'].values))

index_order=[list(highsalt_lownut_combo), list(highsalt_combo_no_lownut), list(lownut_combo_no_highsalt),list(highsalt_only)
             ,list(lownut_only), list(combo_only)]

index_order_flat = [element for innerList in index_order for element in innerList]
df_final=df_final.reindex(index_order_flat)
df_final = df_final[['highsalt','combo','lownut']]

formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))
labels=list(df_final.index.values)

plt.figure(figsize = (12,8))
go_terms=sns.heatmap(df_final,cbar_kws={'format':formatter}, yticklabels=labels)
go_terms.figure.savefig(wd+'go_heatmap.png', dpi='figure', format="png",bbox_inches='tight')


#go_terms=sns.heatmap(df_final,cbar_kws={'format':formatter},cmap="Greens",center=0.00001, yticklabels=labels)
#go_terms.set_title("BP")
#go_terms.figure.savefig(wd+'go_heatmap_bp.png', dpi='figure', format="png",bbox_inches='tight')


#go_terms=sns.heatmap(df_final,cbar_kws={'format':formatter},cmap="Blues",center=0.00001, yticklabels=labels)
#go_terms.set_title("MF")
#go_terms.figure.savefig(wd+'go_heatmap_mf.png', dpi='figure', format="png",bbox_inches='tight')


#go_terms=sns.heatmap(df_final,cbar_kws={'format':formatter},cmap="BuPu",center=0.00001)
#go_terms.set_title("CC")
#go_terms.figure.savefig(wd+'go_heatmap_cc.png', dpi='figure', format="png",bbox_inches='tight')
plt.close()

