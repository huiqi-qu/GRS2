import argparse
import re
from collections import defaultdict
import pandas as pd
from pandas.io.parsers import read_table
import numpy as np
import math
import read_tables

parser=argparse.ArgumentParser()
parser.add_argument("-p", "--prefix", help="plink file path prefix /dir/file of the file /dir/file.ped")
parser.add_argument("-aa", "--AfricanAmerican", action="store_true")
parser.add_argument("-t", "--table_path", help="the folder where the tables are")
parser.add_argument("-o", "--output", help="output name")
args = parser.parse_args()

dq_dict = read_tables.read_DQ(args.table_path)
dqi_dict = read_tables.read_DQ_interations(args.table_path)
non_HLA_DQ = read_tables.read_non_DQ(args.table_path)
aa_snps = defaultdict(lambda: None)
if args.AfricanAmerican:
    aa_snps = read_tables.read_aa_snps(args.table_path, args.AfricanAmerican)

dq_map_list = []
dq_map_col_indices = []
ndq_map_list = []
ndq_map_col_indices = []
aa_map_list = []
aa_map_col_indices = []

map_f = open(args.prefix+".map", "r")
#records only the 
#dq annotation
counter = 6

map_l = map_f.readline()
while map_l:
    r_map_l = re.sub("\n", "", map_l).split("\t")
    if dq_dict[r_map_l[1]] is not None:
        dq_map_list.append(r_map_l[1])
        dq_map_col_indices.append(counter)
        dq_map_col_indices.append(counter + 1)
    elif non_HLA_DQ[r_map_l[1]] is not None:
        ndq_map_list.append(r_map_l[1])
        ndq_map_col_indices.append(counter)
        ndq_map_col_indices.append(counter + 1)
    elif aa_snps[r_map_l[1]] is not None:
        aa_map_list.append(r_map_l[1])
        aa_map_col_indices.append(counter)
        aa_map_col_indices.append(counter + 1)
    map_l = map_f.readline()
    counter += 2
map_f.close()


ped_df = pd.read_csv(args.prefix+".ped", header=None, sep=" ")
dq_df = ped_df[dq_map_col_indices].to_numpy()
ndq_df = ped_df[ndq_map_col_indices].to_numpy()
aa_df = 0
if args.AfricanAmerican:
    aa_df = ped_df[aa_map_col_indices].to_numpy()

ids = list(ped_df[[1]].to_numpy().transpose())
scores =[0] * len(dq_df)
for i in range(len(dq_df)):
    dq_vals = []
    dq_haps = []
    for j in range(len(dq_df[i])): 
        if dq_df[i][j] == dq_dict[dq_map_list[math.floor(j/2)]][1]:
            dq_haps.append(dq_dict[dq_map_list[math.floor(j/2)]][2])
            dq_vals.append(dq_dict[dq_map_list[math.floor(j/2)]][3])
    if dqi_dict[",".join(dq_haps)] is not None:
        scores[i] += dqi_dict[",".join(dq_haps)]
    else:
        scores[i] += sum(dq_vals)
    for j in range(len(ndq_df[i])):
        if ndq_df[i][j] == non_HLA_DQ[ndq_map_list[math.floor(j/2)]][1]:
            scores[i] += non_HLA_DQ[ndq_map_list[math.floor(j/2)]][2]
    if args.AfricanAmerican:
        for j in range(len(aa_df[i])):
            if aa_df[i][j] == aa_snps[aa_map_list[math.floor(j/2)]][1]:
                scores[i] += aa_snps[aa_map_list[math.floor(j/2)]][2]

ind_df = ped_df[[0,1]]
ind_df = ind_df.rename(columns={0 : "FID", 1: "IID"})
ind_df["scores"] = scores
ind_df.to_csv(args.output, index=False, header=False)


