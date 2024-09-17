# Name: create bedfile
# Author: EY
# Date: May 01 2023 (edited October 13 2023)
# Version: Python 3.9
# Description: create bedfile for bedtools input 'getfasta' for HAshoot gtf file

import pandas as pd
import re
# create bedfile for bedtools input 'getfasta'

# chromosome start end name
chrom= []
start=[]
end=[]
name=[]
with open('/scratch/ely67071/sunflower_inflo_dev_data_b3/genomes/Ha412HO_v2_braker.gff3') as topo_file:
    for line in topo_file:
        if '#' not in line:

            split_line = re.split(r'\t+', line)
            if "mRNA" in split_line[2]:
                start.append(int(split_line[3]))
                end.append(int(split_line[4]))
                split_quote = split_line[8].split('=')
                split_quote2=split_quote[1].split(';')
                mrna_id = split_quote2[0]
                name.append(mrna_id)
                chrom_name=split_line[0]
                chrom.append(chrom_name)

d = {'Chrom': chrom, 'Start': start, 'End':end,'Name':name}
df = pd.DataFrame(d)
df.to_csv('/scratch/ely67071/sunflower_inflo_dev_data_b3/interpro/getfasta_bed_input_b3.bed', sep='\t',header=False,index=False)
