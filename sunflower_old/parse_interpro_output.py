# Name: parse interpro output
# Author: EY
# Date: May 16 2023
# Version: Python 3.9
# Description: Will extract mRNA ID and GO Terms from interpro gff3 output

import glob
import re
import pandas as pd

# absolute path to search all text files inside a specific folder
path = r'/scratch/ely67071/sunflower_inflo_dev_data_b3/interpro/interpro_scripts/*.fa.gff3'
files = glob.glob(path)

name = []
go_terms=[]
counter=0
for file in files:
    counter2=0
    with open(file) as topo_file:
        for line in topo_file:
            if "GO:" in line:
                terms=[]
                split_line = re.split(r'\t+', line)
                split_name=split_line[0].split("::")
                name.append(split_name[0])
                split_go=re.split(r';',split_line[8])
                split_id=re.split(r'"',split_go[2])
                for item in split_id:
                    if "GO" in item:
                        terms.append(item)
                string_go = ','.join(terms)
                go_terms.append([terms])
                counter2+=1
    counter+=1
    print(counter2)

d = {'ID': name, 'go_terms': go_terms}
df = pd.DataFrame(d)
df_list = df.values.tolist()

go_dict={}
for i in df_list:
    go_dict[i[0]] = []
for j in df_list:
    for i in j[1]:
        for k in i:
            if k not in go_dict[j[0]]:
                go_dict[j[0]].append(k)
for k,v in go_dict.items():
    string_go = ','.join(v)
    go_dict[k] = string_go

df2=pd.DataFrame(go_dict.items())

df2.to_csv('Ha412GO_terms_interpro_b3.txt', sep='\t',header=False,index=False)
