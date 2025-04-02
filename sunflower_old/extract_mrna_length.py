# Name: extract mrna length
# Author: EY
# Date: December 10 2024
# Version: Python 3.9
# Description: Will extract mRNA length from gff3 annotation file


import re
import pandas as pd
name_list = []
difference_list = []
with open('/scratch/ely67071/sunflower_inflo_dev_data_b3/genomes/Ha412HO_v2_braker.gff3') as topo_file:
    for line in topo_file:
        if "#" not in line:
            split_line = re.split(r'\t+', line)
            if "mRNA" in split_line[2]:
                difference = (int(split_line[4])-int(split_line[3]))+1
                split_equal = split_line[8].split('=')
                split_semi = split_equal[1].split(';')
                mrna_id = split_semi[0]
                print(mrna_id,difference)
                difference_list.append(difference)
                name_list.append(mrna_id)
d = {'ID': name_list, 'Length': difference_list}
df = pd.DataFrame(d)
df.to_csv('/home/ely67071/dev_RNAseq/sunflower/transcript_length.csv', sep='\t')

