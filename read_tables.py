from collections import defaultdict
import re

def read_DQ(table_path):
    dq_dict = defaultdict(lambda: None)
    dq_f = open(table_path + "/DQ_haplotypes.txt", "r")
    dq_l = dq_f.readline()
    dq_l = dq_f.readline()
    dq_count = 0
    while dq_l:
        r_dq_l = re.sub("\n", "", dq_l).split("\t")
        dq_dict[r_dq_l[1]] = (r_dq_l[0], r_dq_l[2], r_dq_l[3], float(r_dq_l[4]))
        dq_l = dq_f.readline()
        dq_count += 1
    dq_f.close()
    return dq_dict

def read_DQ_interations(table_path):
    dqi_dict = defaultdict(lambda: None)
    dqi_f = open(table_path + "/DQ_interactions.txt", "r")
    dqi_l = dqi_f.readline()
    dqi_l = dqi_f.readline()
    dqi_count = 0
    while dqi_l:
        r_dqi_l = re.sub("\n", "", dqi_l).split("\t")
        dqi_dict[r_dqi_l[0]] = float(r_dqi_l[1])
        dqi_l = dqi_f.readline()
        dqi_count += 1
    dqi_f.close()
    return dqi_dict

def read_non_DQ(table_path):
    non_HLA_DQ = defaultdict(lambda: None)
    ndq_f = open(table_path + "/non_HLA_DQ.txt", "r")
    ndq_l = ndq_f.readline()
    ndq_l = ndq_f.readline()
    ndq_l = ndq_f.readline()
    ndq_count = 0
    while ndq_l:
        r_ndq_l = re.sub("\n","", ndq_l).split("\t")
        non_HLA_DQ[r_ndq_l[1]] = (r_ndq_l[0], r_ndq_l[2], float(r_ndq_l[3])) #by marker
        ndq_l = ndq_f.readline()
        ndq_count += 1
    ndq_f.close()
    return non_HLA_DQ

def read_aa_snps(table_path, AfricanAmerican):
    aa_snps = defaultdict(lambda: None)
    aa_count = 0
    if AfricanAmerican:
        aa_f = open(table_path + "/AfricanSNPS.txt", "r")
        aa_l = aa_f.readline()
        aa_l = aa_f.readline()
        aa_l = aa_f.readline()
        while aa_l:
            r_aa_l = re.sub("\n", "", aa_l).split()
            aa_snps[r_aa_l[1]] = (r_aa_l[0], r_aa_l[2],float(r_aa_l[3]))
            aa_l = aa_f.readline()
            aa_count += 1
    return aa_snps
