# Name: add orthogroup
# Author: EY
# Date: May 15 2024
# Version: Python 3.10
# Description: add Arabidopsis orthogroup to DESeq output
import os
import pandas as pd

def process_files(directory, output_filename):
    """function will read in all .csv files in a directory, add Arabidopsis orthogroups and Orthogroup, and output as a single excel
    file with different sheets"""
    dataframes = {}
    orthogroups = pd.read_csv("Aug30_Orthogroups_longest_iso.csv", sep=",")
    column_to_search = 'Lactuca'
    values_lists_orthologs = []
    values_lists_orthogroups = []

    for filename in os.listdir(directory):
        if filename.endswith('.csv'):
            filepath = os.path.join(directory, filename)
            df = pd.read_csv(filepath)
            df_subset = df[df['padj'] < 0.05]
            dataframes[filename] = df_subset

            values_list_orthologs = []
            values_list_orthogroups = []
            for index, row in df_subset.iterrows():
                value_to_check = row['Unnamed: 0']
                print(value_to_check)
                value_found = False
                for separate_index, separate_row in orthogroups.iterrows():
                    try:
                        if pd.notna(separate_row[column_to_search]) and value_to_check in separate_row[column_to_search]:
                            # Add both Arabidopsis orthologs and Orthogroup
                            ortholog_value = separate_row['Arabidopsis']
                            orthogroup_value = separate_row['Orthogroup']  # Assuming "Orthogroup" is the column name
                            values_list_orthologs.append(ortholog_value)
                            values_list_orthogroups.append(orthogroup_value)
                            value_found = True
                            break
                    except TypeError:
                        values_list_orthologs.append('nan')
                        values_list_orthogroups.append('nan')
                if not value_found:
                    values_list_orthologs.append('nan')
                    values_list_orthogroups.append('nan')

            values_lists_orthologs.append(values_list_orthologs)
            values_lists_orthogroups.append(values_list_orthogroups)

    new_column_name_orthologs = 'Arabidopsis_orthologs'
    new_column_name_orthogroups = 'Orthogroup'

    for df, orthologs_list, orthogroups_list in zip(dataframes.values(), values_lists_orthologs, values_lists_orthogroups):
        df[new_column_name_orthologs] = orthologs_list
        df[new_column_name_orthogroups] = orthogroups_list

    excel_writer = pd.ExcelWriter(output_filename, engine='xlsxwriter')

    for df_name, df in dataframes.items():
        df.to_excel(excel_writer, sheet_name=df_name, index=False)

    excel_writer.save()
    excel_writer.close()


# process for lettuce
process_files("lettuce/pairwise/", "lettuce/staged_lettuce_w_orthogroup.xlsx")

