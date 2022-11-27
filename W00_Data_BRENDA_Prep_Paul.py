#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# The following code ensures the code work properly in 
# MS VS, MS VS CODE and jupyter notebook on both Linux and Windows.
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
if __name__ == "__main__":
    print("="*80)
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")
#--------------------------------------------------#
    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")
#--------------------------------------------------#
    if "__file__" in globals().keys():
        print('CurrentDir: ', os.getcwd())
        try:
            os.chdir(os.path.dirname(__file__))
        except:
            print("Problems with navigating to the file dir.")
        print('CurrentDir: ', os.getcwd())
    else:
        print("Running in python jupyter notebook.")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            print("Problems with navigating to the workbook dir.")
#--------------------------------------------------#

###################################################################################################################
###################################################################################################################
# Imports
import pickle
import argparse
import numpy as np
import pandas as pd

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Basic Functions
def beautiful_print(df): # Print the DataFrame obtained.
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows', 20, 
                           'display.min_rows', 20, 
                           'display.max_columns', 10, 
                           #"display.max_colwidth", None,
                           "display.width", None,
                           "expand_frame_repr", True,
                           "max_seq_items", None,):  # more options can be specified
        # Once the display.max_rows is exceeded, 
        # the display.min_rows options determines 
        # how many rows are shown in the truncated repr.
        print(df)
    return 

def get_longest_smiles_from_aggregates(one_smiles_string):
    return max(one_smiles_string.split("."), key = len)

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# Load database link files.
#====================================================================================================#
# Get dict for reading KEGG compound IDs.
KEGG_id_SMILES_by_ChEBI_CAS_df = pd.read_csv("./CMPD_RXN_tools/KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_ChEBI_CAS.csv", header = 0)

KEGGid_list     = KEGG_id_SMILES_by_ChEBI_CAS_df["KEGGid"].tolist()
SMILES_by_ChEBI = KEGG_id_SMILES_by_ChEBI_CAS_df["SMILES_by_ChEBI"].tolist()
SMILES_by_CAS   = KEGG_id_SMILES_by_ChEBI_CAS_df["SMILES_by_CAS"].tolist()

KEGG_id_SMILES_by_ChEBI_CAS_dict = dict([])
for i, id in enumerate(KEGGid_list):
    KEGG_id_SMILES_by_ChEBI_CAS_dict[id] = [SMILES_by_ChEBI[i], SMILES_by_CAS[i]]

#====================================================================================================#
# Get dict for reading KEGG compound IDs.

KEGG_id_SMILES_by_mol_df = pd.read_csv("./CMPD_RXN_tools/KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_mol_file_no_uniq.csv", header = 0)

KEGGid_list     = KEGG_id_SMILES_by_mol_df["KEGGid"].tolist()
SMILES_by_mol   = KEGG_id_SMILES_by_mol_df["SMILES_by_mol"].tolist()


KEGG_id_SMILES_by_mol_dict = dict([])
for i, id in enumerate(KEGGid_list):
    KEGG_id_SMILES_by_mol_dict[id] = SMILES_by_mol[i]

#====================================================================================================#
# Get dict for reading ALL_KEGG_DATA_DICT.

folder = Path("./CMPD_RXN_tools/KEGGScrapSavings")
filename = "CPD02_ALL_KEGG_DATA_DICT_bk.p"

pickle_in1 = open( folder / filename, "rb")
ALL_KEGG_DATA_DICT = pickle.load(pickle_in1)
pickle_in1.close()

#====================================================================================================#
# Get dict for reading ALL_KEGG_DATA_DICT.

folder = Path("./CMPD_RXN_tools/MNXScrapSavings")
filename = "CPD01_nme_smiles_MNXid_dict.p"

pickle_in1 = open( folder / filename, "rb")
nme_smiles_MNXid_dict = pickle.load(pickle_in1)
pickle_in1.close()


#====================================================================================================#
# KEGG_DATA
KEGG_name_id_dict = dict([])
for one_KEGG_id in ALL_KEGG_DATA_DICT:
    if "NAME" in ALL_KEGG_DATA_DICT[one_KEGG_id]:
        for one_name in ALL_KEGG_DATA_DICT[one_KEGG_id]["NAME"]:
            KEGG_name_id_dict[one_name.lower()] = one_KEGG_id
            
#====================================================================================================#
# cmpd_smls_dict from cactus
cmpd_smls_cactus_df = pd.read_csv("X_DataProcessing/X00_enzyme_datasets/Km_3_BRENDA/cmpd_smls_cactus.csv",
                       encoding = "ISO-8859-1",
                       index_col = None, 
                       header = 0, 
                       sep = ';',
                       engine = 'python')

cmpd_name_list_cactus = cmpd_smls_cactus_df["CMPD"].tolist()
cmpd_smls_list_cactus = cmpd_smls_cactus_df["CMPD_SMILES"].tolist()

cactus_cmpd_smls_dict = dict([])
for i, one_cmpd_name in enumerate(cmpd_name_list_cactus):
    cactus_cmpd_smls_dict[str(one_cmpd_name).lower()] = cmpd_smls_list_cactus[i]

#====================================================================================================#
# cmpd_smls_dict from pubchem
cmpd_smls_pubchem_df = pd.read_csv("X_DataProcessing/X00_enzyme_datasets/Km_3_BRENDA/cmpd_smls_pubchem.csv",
                       encoding = "ISO-8859-1",
                       index_col = None, 
                       header = 0, 
                       sep = ';',
                       engine = 'python')

cmpd_name_list_pubchem = cmpd_smls_pubchem_df["CMPD"].tolist()
cmpd_smls_list_pubchem = cmpd_smls_pubchem_df["CMPD_SMILES"].tolist()

pubchem_cmpd_smls_dict = dict([])
for i, one_cmpd_name in enumerate(cmpd_name_list_pubchem):
    pubchem_cmpd_smls_dict[str(one_cmpd_name).lower()] = cmpd_smls_list_pubchem[i]

###################################################################################################################
###################################################################################################################
# Deal with sequences
AC_SQ_df = pd.read_csv("X_DataProcessing/X00_enzyme_datasets/uniprot_AC_SQ_2.csv", 
                       index_col = None, 
                       header = 0, 
                       sep = ",")

entr_list = AC_SQ_df["Entry"].tolist()
seqs_list = AC_SQ_df["Sequence"].tolist()
#unip_list = AC_SQ_df["UNIPROT ID"].tolist()

AC_SQ_dict = dict([])
for i, entr in enumerate(entr_list):
    AC_SQ_dict[entr] = seqs_list[i]

del seqs_list
del entr_list
#del unip_list

###################################################################################################################
###################################################################################################################
# Deal with sequences
Seq_EC_Org_df = pd.read_csv("X_DataProcessing/X00_enzyme_datasets/Seq_EC_Org_KINMOD.csv", 
                       index_col = None, 
                       header = 0, 
                       sep = ",")

ecno_list = Seq_EC_Org_df["EC"].tolist()
seqs_list = Seq_EC_Org_df["SEQ"].tolist()
orgm_list = Seq_EC_Org_df["ORG"].tolist()

Seq_EC_Org_dict = dict([])
for i, (ecno, orgm) in enumerate(zip(ecno_list, orgm_list)):
    if (ecno, str(orgm).lower()) not in Seq_EC_Org_dict:
        Seq_EC_Org_dict[(ecno, str(orgm).lower())] = [seqs_list[i],]
    else:
        Seq_EC_Org_dict[(ecno, str(orgm).lower())].append(seqs_list[i])


del seqs_list
del ecno_list
del orgm_list

###################################################################################################################
###################################################################################################################
# Get Raw KM file.
Km_df = pd.read_csv("X_DataProcessing/X00_enzyme_datasets/Km_3_BRENDA/brenda_km_raw.csv", 
                    index_col = None, 
                    names = ["EC", "Km", "-", "cmpd", "condition", "organism", "uniprot_id", "publication", "other"], 
                    header = None, 
                    sep = "\t")


#Km_df = Km_df.dropna()
Km_df = Km_df.reset_index()
beautiful_print(Km_df)

ecno_list = Km_df["EC"].tolist()
prpt_list = Km_df["Km"].tolist()
orgm_list = Km_df["organism"].tolist()
unip_list = Km_df["uniprot_id"].tolist()
cmpd_list = Km_df["cmpd"].tolist()
rctt_list = cmpd_list
seqs_list = unip_list


#====================================================================================================#
# 
print("len(rctt_list): ", len(rctt_list))
print("len(prpt_list): ", len(prpt_list))
print("len(orgm_list): ", len(orgm_list))

###################################################################################################################
###################################################################################################################
# Put all cmpd-smls dict together
smls_list = []
cmpd_unknown_list = []
for one_name in rctt_list:
    append_bool = False
    one_name = str(one_name)
    one_name = one_name.replace(";", "")
    one_name = one_name.lower()

    if not append_bool:
        if one_name in pubchem_cmpd_smls_dict:
            cmpd_smls = pubchem_cmpd_smls_dict[one_name]
            if cmpd_smls not in ["NA", "", "nan", "None"]:
                smls_list.append(cmpd_smls)
                append_bool = True

    if not append_bool:
        if one_name in cactus_cmpd_smls_dict:
            cmpd_smls = cactus_cmpd_smls_dict[one_name]
            if cmpd_smls not in ["NA", "", "nan", "None"]:
                smls_list.append(cmpd_smls)
                append_bool = True

    if not append_bool:
        if one_name in nme_smiles_MNXid_dict:
            cmpd_smls = nme_smiles_MNXid_dict[one_name][0]
            if cmpd_smls not in ["NA", "", "nan", "None"]:
                smls_list.append(cmpd_smls)
                append_bool = True

    if not append_bool:
        if one_name in KEGG_name_id_dict:
            cmpd_id = KEGG_name_id_dict[one_name]
            if cmpd_id not in ["NA", "", "nan", "None"]:
                one_KEGG_id = cmpd_id
                if one_KEGG_id in KEGG_id_SMILES_by_mol_dict:
                    if KEGG_id_SMILES_by_mol_dict[one_KEGG_id] != "None":
                        smls_list.append(KEGG_id_SMILES_by_mol_dict[one_KEGG_id])
                        append_bool = True
                        
    if not append_bool:
        smls_list.append("None")
        cmpd_unknown_list.append(one_name)
        #print("cmpd:"+one_name+", smls:None")

'''
with open('X_DataProcessing/X00_enzyme_datasets/Km_3_BRENDA/cmpd_smls_all.csv', 'w') as f:
    f.write(";CMPD;CMPD_SMILES"+"\n")
    for i, (cmpd, smls) in enumerate(zip(cmpd_list, smls_list)):
        print(cmpd)
        if type(cmpd) != str:
            cmpd = str(cmpd)

        f.write(str(i)+";"+cmpd+";"+str(smls)+"\n")
        print(i, cmpd, smls)
        '''

print("\n", smls_list.count("None"), "number of names cannot find smiles!")


print("\n" + str(len(set(rctt_list))), "number of compound names w/o duplicates")
print(len(set(cmpd_unknown_list)), "number of names w/o duplicates cannot find smiles!")




print("len(smls_list): ", len(smls_list))
print("len(rctt_list): ", len(rctt_list))
print("len(unip_list): ", len(unip_list))
print("len(prpt_list): ", len(prpt_list))
print("len(orgm_list): ", len(orgm_list))

all_list = []
new_smls_list = []
new_rctt_list = []
new_seqs_list = []
new_prpt_list = []
new_orgm_list = []
new_ecno_list = []

for i, (cmpd, smls, unip, prpt, orgm, ecno) in enumerate( zip(rctt_list, smls_list, unip_list, prpt_list, orgm_list, ecno_list) ):
    if smls != "None" and isfloat(prpt):
        if unip != "-": # known uniprot accession no. 
            if unip.find(", ") == -1 and unip.find(" AND ") == -1:
                if unip in AC_SQ_dict:
                    seqs = AC_SQ_dict[unip]
                    all_list.append(tuple([smls, seqs, prpt, orgm]))
                    new_smls_list.append(smls)
                    new_rctt_list.append(cmpd)
                    new_seqs_list.append(seqs)
                    new_prpt_list.append(prpt)
                    new_orgm_list.append(orgm)
                    new_ecno_list.append(ecno)
                    #print(unip)

            if unip.find(", ") != -1:
                for unip_x in unip.split(", "):
                    if unip_x in AC_SQ_dict:
                        seqs = AC_SQ_dict[unip_x]
                        all_list.append(tuple([smls, seqs, prpt, orgm]))
                        new_smls_list.append(smls)
                        new_rctt_list.append(cmpd)
                        new_seqs_list.append(seqs)
                        new_prpt_list.append(prpt)
                        new_orgm_list.append(orgm)
                        new_ecno_list.append(ecno)
                        #print(unip_x)

            if unip.find(" AND ") != -1:
                for unip_x in unip.split(" AND "):
                    if unip_x in AC_SQ_dict:
                        seqs = AC_SQ_dict[unip_x]
                        all_list.append(tuple([smls, seqs, prpt, orgm]))
                        new_smls_list.append(smls)
                        new_rctt_list.append(cmpd)
                        new_seqs_list.append(seqs)
                        new_prpt_list.append(prpt)
                        new_orgm_list.append(orgm)
                        new_ecno_list.append(ecno)
                        #print(unip_x)
        
        if unip == "-": # unknown uniprot accession no.
            ec_og = (ecno, orgm.lower())
            if ec_og in Seq_EC_Org_dict:
                ec_og_seqs_list = Seq_EC_Org_dict[ec_og]
                for seqs in ec_og_seqs_list:
                    all_list.append(tuple([smls, seqs, prpt, orgm]))
                    new_seqs_list.append(seqs)
                    new_smls_list.append(smls)
                    new_rctt_list.append(cmpd)
                    new_prpt_list.append(prpt)
                    new_orgm_list.append(orgm)
                    new_ecno_list.append(ecno)
                    

print("len(set(all_list)): ", len(set(all_list)))






#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Generate processed data file.

pair_list = [(smls, seqs) for smls, seqs in zip(new_smls_list, new_seqs_list)]
print("len(set(pair_list)): ", len(set(pair_list)))


orgm_count_dict = {}
orgm_set_list = list(set(new_orgm_list))
for one_orgm in orgm_set_list:
    orgm_count = new_orgm_list.count(one_orgm)
    orgm_count_dict[one_orgm] = orgm_count

orgm_count_dict_sorted = {k: v for k, v in sorted(orgm_count_dict.items(),reverse = True, key=lambda item: item[1])}
orgm_count_list_sorted = list(orgm_count_dict_sorted.keys())

for one_orgm in orgm_count_list_sorted:
    pass


pair_set_list = list(set(pair_list))
pair_orgm_prpt_dict = dict([])


for i, (seqs, cmpd, prpt, orgm, pair) in enumerate( zip(new_seqs_list, new_smls_list, new_prpt_list, new_orgm_list, pair_list) ):
    if pair not in pair_orgm_prpt_dict.keys():
        pair_orgm_prpt_dict[pair] = [tuple([orgm, prpt]), ]
    else:
        pair_orgm_prpt_dict[pair].append(tuple([orgm, prpt])) 

print("len(pair_orgm_prpt_dict): ", len(pair_orgm_prpt_dict))


pair_orgm_prpt_dict_select = dict([])
for pair in list(pair_orgm_prpt_dict.keys()):
    if len(pair_orgm_prpt_dict[pair]) != 1:
        all_orgm_of_pair_list = [orgm_prpt[0] for orgm_prpt in pair_orgm_prpt_dict[pair]]
        orgm = "Cant be me!"
        for one_orgm in orgm_count_list_sorted:
            if one_orgm in all_orgm_of_pair_list:
                orgm = one_orgm
                break   
        pair_one_orgm_prpt_list = []
        for orgm_prpt in pair_orgm_prpt_dict[pair]:
            if orgm_prpt[0] == orgm:
                pair_one_orgm_prpt_list.append(orgm_prpt[1])
        pair_orgm_prpt_dict_select[pair] = max(pair_one_orgm_prpt_list)
    else:
        pair_orgm_prpt_dict_select[pair] = pair_orgm_prpt_dict[pair][0][1]

print("len(pair_orgm_prpt_dict_select): ", len(pair_orgm_prpt_dict_select))
#print(pair_orgm_prpt_dict_select)


with open('X_DataProcessing/X00_enzyme_datasets_processed/KM_BRENDA.csv' , 'w') as f:
    count_x=0
    f.write(",SEQ,CMPD_SMILES,Log_Km_Value"+"\n")
    for one_datapoint in list(pair_orgm_prpt_dict_select.keys()):
        # write a row to the csv file
        f.write(str(count_x)+","+str(one_datapoint[1])+","+str(one_datapoint[0])+","+str(pair_orgm_prpt_dict_select[one_datapoint])+"\n")
        count_x+=1

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


pair_list = [(smls, seqs) for smls, seqs in zip(new_smls_list, new_seqs_list)]

orgm_count_dict = {}
orgm_set_list = list(set(new_orgm_list))
for one_orgm in orgm_set_list:
    orgm_count = new_orgm_list.count(one_orgm)
    orgm_count_dict[one_orgm] = orgm_count

orgm_count_dict_sorted = {k: v for k, v in sorted(orgm_count_dict.items(),reverse = True, key=lambda item: item[1])}
orgm_count_list_sorted = list(orgm_count_dict_sorted.keys())

for one_orgm in orgm_count_list_sorted:
    pass


pair_set_list = list(set(pair_list))
pair_orgm_prpt_dict = dict([])


for i, (seqs, cmpd, prpt, orgm, pair) in enumerate( zip(new_seqs_list, new_smls_list, new_prpt_list, new_orgm_list, pair_list) ):
    if pair not in pair_orgm_prpt_dict.keys():
        pair_orgm_prpt_dict[pair] = [tuple([orgm, prpt]), ]
    else:
        pair_orgm_prpt_dict[pair].append(tuple([orgm, prpt])) 

print("len(pair_orgm_prpt_dict): ", len(pair_orgm_prpt_dict))


pair_orgm_prpt_dict_select = dict([])
for pair in list(pair_orgm_prpt_dict.keys()):
    if len(pair_orgm_prpt_dict[pair]) != 1:
        all_orgm_of_pair_list = [orgm_prpt[0] for orgm_prpt in pair_orgm_prpt_dict[pair]]
        orgm = "Cant be me!"
        for one_orgm in orgm_count_list_sorted:
            if one_orgm in all_orgm_of_pair_list:
                orgm = one_orgm
                break   
        pair_one_orgm_prpt_list = []
        for orgm_prpt in pair_orgm_prpt_dict[pair]:
            if orgm_prpt[0] == orgm:
                pair_one_orgm_prpt_list.append(orgm_prpt[1])
        pair_orgm_prpt_dict_select[pair] = max(pair_one_orgm_prpt_list)
    else:
        pair_orgm_prpt_dict_select[pair] = pair_orgm_prpt_dict[pair][0][1]

print("len(pair_orgm_prpt_dict_select): ", len(pair_orgm_prpt_dict_select))
#print(pair_orgm_prpt_dict_select)


with open('X_DataProcessing/X00_enzyme_datasets_processed/KM_BRENDA.csv' , 'w') as f:
    count_x=0
    f.write(",SEQ,CMPD_SMILES,Log_Km_Value"+"\n")
    for one_datapoint in list(pair_orgm_prpt_dict_select.keys()):
        # write a row to the csv file
        f.write(str(count_x)+","+str(one_datapoint[1])+","+str(one_datapoint[0])+","+str(pair_orgm_prpt_dict_select[one_datapoint])+"\n")
        count_x+=1


print("Done!")

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#




#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #
#   `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'       `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'  #
#     VAMAV          VAMAV          VAMAV          VAMAV          VAMAV           VAMAV          VAMAV          VAMAV          VAMAV          VAMAV    #
#      VVV            VVV            VVV            VVV            VVV             VVV            VVV            VVV            VVV            VVV     #
#       V              V              V              V              V               V              V              V              V              V      #

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
###################################################################################################################
###################################################################################################################
#====================================================================================================#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#--------------------------------------------------#
#------------------------------

#                                                                                                                                                          
#      `MM.              `MM.             `MM.             `MM.             `MM.             `MM.             `MM.             `MM.             `MM.       
#        `Mb.              `Mb.             `Mb.             `Mb.             `Mb.             `Mb.             `Mb.             `Mb.             `Mb.     
# MMMMMMMMMMMMD     MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD   
#         ,M'               ,M'              ,M'              ,M'              ,M'              ,M'              ,M'              ,M'              ,M'     
#       .M'               .M'              .M'              .M'              .M'              .M'              .M'              .M'              .M'       
#                                                                                                                                                          

#------------------------------
#--------------------------------------------------#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#====================================================================================================#
###################################################################################################################
###################################################################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#       A              A              A              A              A               A              A              A              A              A      #
#      MMM            MMM            MMM            MMM            MMM             MMM            MMM            MMM            MMM            MMM     #
#     MMMMM          MMMMM          MMMMM          MMMMM          MMMMM           MMMMM          MMMMM          MMMMM          MMMMM          MMMMM    #
#   ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.       ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.  #
#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #




