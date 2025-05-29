import pandas as pd

def add_gene_info_from_excel(tsv_file, excel_file, output_excel):
    # Step 1: Load the tsv file with genes of interest
    genes_df = pd.read_csv(tsv_file, sep="\t", encoding="ISO-8859-1")

    # Step 2: Load the Excel file with the orthogroups data (all sheets)
    excel_df = pd.read_excel(excel_file, sheet_name=None)  # Load all sheets
    sheet_names = excel_df.keys()  # Get sheet names

    # Mapping for day values to corresponding labels
    day_to_label = {
        '10D': 'VM',
        '20D': 'TM',
        '30D': 'IM',
        '35D': 'IM/FM'
    }

    # Step 3: Create lists to store new columns' values
    gene_in_h_annuus = []
    stage_de_h_annuus = []

    # Step 4: Iterate through each row in the genes_df and check for OrthogroupID
    for _, gene_row in genes_df.iterrows():
        orthogroup_id = gene_row['OrthogroupID']

        # Initialize empty lists for this row's new columns
        gene_matches = []
        stage_matches = []

        # Step 5: Check each sheet in the Excel file for matches
        for sheet_name in sheet_names:
            df = excel_df[sheet_name]  # Get dataframe for the current sheet
            orthogroups_column = df['Orthogroup']  # The Orthogroup column in the sheet

            # Check if the OrthogroupID is in this sheet's Orthogroup column
            if orthogroup_id in orthogroups_column.values:
                # Get the matching Unnamed: 0 value
                matched_gene = df.loc[orthogroups_column == orthogroup_id, 'Unnamed: 0'].values[0]
                gene_matches.append(matched_gene)

                # Clean up the sheet name: remove 'result_' and '.csv'
                cleaned_sheet_name = sheet_name.replace('result_', '').replace('.csv', '')

                # Replace day values (e.g., 10D) with corresponding labels (e.g., VM)
                for day, label in day_to_label.items():
                    if day in cleaned_sheet_name:
                        cleaned_sheet_name = cleaned_sheet_name.replace(day, label)

                # Append the cleaned and mapped sheet name
                stage_matches.append(cleaned_sheet_name)

        # Step 6: Remove duplicates from gene_matches (by converting to a set) and join by commas
        if gene_matches:
            unique_genes = list(set(gene_matches))  # Remove duplicates
            gene_in_h_annuus.append(', '.join(unique_genes))  # Join by comma if multiple matches
            stage_de_h_annuus.append(', '.join(stage_matches))  # Join sheet names by comma if multiple
        else:
            # If no match, append 'nan'
            gene_in_h_annuus.append('.')
            stage_de_h_annuus.append('.')

    # Step 7: Add the new columns to the genes dataframe
    genes_df['Gene_in_H_annuus'] = gene_in_h_annuus
    genes_df['Stage_DE_H_annuus'] = stage_de_h_annuus

    # Step 8: Save the updated dataframe to a new Excel file
    genes_df.to_excel(output_excel, index=False)


# Example of usage:
add_gene_info_from_excel('pairwise/genes_of_interest_transformed.tsv', 'pairwise/staged_sunflower_w_orthogroup.xlsx',
                         'pairwise/genes_of_interest_transformed_sunflower_data.xlsx')
