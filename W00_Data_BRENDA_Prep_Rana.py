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
#--------------------------------------------------#
import re
import time
import copy
import pickle
import argparse
import numpy as np
import pandas as pd
#--------------------------------------------------#
import requests
import xmltodict
#--------------------------------------------------#
from timeit import timeit
#--------------------------------------------------#
import urllib
import xml.etree.ElementTree as ET
from urllib.request import urlopen

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Basic Functions
# Print the DataFrame obtained.
def beautiful_print(df): # Print the DataFrame obtained.
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'       , 20   , 
                           'display.min_rows'       , 20   , 
                           'display.max_columns'    , 5    , 
                           #"display.max_colwidth"  , None ,
                           "display.width"          , None ,
                           "expand_frame_repr"      , True ,
                           "max_seq_items"          , None , ):  # more options can be specified
        # Once the display.max_rows is exceeded, 
        # the display.min_rows options determines 
        # how many rows are shown in the truncated repr.
        print(df)
    return 

#--------------------------------------------------#
def get_longest_smiles_from_aggregates(one_smiles_string):
    return max(one_smiles_string.split("."), key = len)

#--------------------------------------------------#
def isfloat(num): 
    # Return True if the input can be converted to a float or is a float already.
    try:
        float(num)
        return True
    except ValueError:
        return False



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                    `7MMF'`7MN.   `7MF'`7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM  .M"""bgd                                                                #
#                      MM    MMN.    M    MM   `MM.  MM       M  P'   MM   `7 ,MI    "Y                                                                #
#   ,pP""Yq.           MM    M YMb   M    MM   ,M9   MM       M       MM      `MMb.                                                                    #
#  6W'    `Wb          MM    M  `MN. M    MMmmdM9    MM       M       MM        `YMMNq.                                                                #
#  8M      M8          MM    M   `MM.M    MM         MM       M       MM      .     `MM                                                                #
#  YA.    ,A9 ,,       MM    M     YMM    MM         YM.     ,M       MM      Mb     dM                                                                #
#   `Ybmmd9'  db     .JMML..JML.    YM  .JMML.        `bmmmmd"'     .JMML.    P"Ybmmd"                                                                 #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
# Input Arguments
data_folder_name    = ["Ki_BRENDA"        , "KM_3_BRENDA"      , "kcat_BRENDA"        , "kcat_KM_BRENDA"        ][1]
data_file_name      = ["brenda_Ki_raw.csv", "brenda_KM_raw.csv", "brenda_kcat_raw.csv", "brenda_kcat_KM_raw.csv"][1]
output_file_name    = ["Ki_BRENDA.csv"    , "KM_BRENDA.csv"    , "kcat_BRENDA.csv"    , "kcat_KM_BRENDA.csv"    ][1]
data_name           = ["Ki"               , "Km"               , "kcat"               , "kcat_KM"               ][1]
data_name_0         = ["max_Ki"           , "max_Km"           , "max_kcat"           , "kcat_KM_0"             ][1]

data_folder      = Path("X_DataProcessing/X00_enzyme_datasets/" + data_folder_name)
data_file        = data_file_name

output_folder    = Path("X_DataProcessing/X00_enzyme_datasets_processed/")
output_file      = output_file_name

val_SD_threshold = 1


###################################################################################################################
###################################################################################################################

path = str(data_folder) + "/"

file_path = path + "brenda_km_raw.csv"

raw = pd.read_csv(file_path, on_bad_lines = 'skip', header = None, sep='	', encoding='cp1252')


filtered_raw = raw[raw[1].str.contains("additional information")==False] # skipped
filtered_raw = filtered_raw[filtered_raw[2].str.contains("additional information")==False] # skipped

max_removed = filtered_raw[filtered_raw[2].str.contains('-')==True]
max_included = filtered_raw[filtered_raw[2].str.contains('-')==False] 

df = max_removed.drop(columns=[2, 4, 7, 8])
df = df.rename(columns={0: "EC Number", 1: "Km [mM]", 3:"Substrate", 5: "Organism", 6: "Uniprot ID"})
df = df.iloc[:, [-2, 0, -1, 2, 1]]

df['Uniprot ID'] = df['Uniprot ID'].str.replace('and', ',')
df['Uniprot ID'] = df['Uniprot ID'].str.replace('AND', ',')
df['Uniprot ID'] = df['Uniprot ID'].str.replace(' ', '')

uni_comp = df['Substrate'].unique()

#%% Looking at dataset with UniProt Available Only

data_wi_unip_df = df[df['Uniprot ID'].str.contains('-')==False]
data_wi_unip_df = data_wi_unip_df[data_wi_unip_df["Km [mM]"].str.contains("additional information")==False] 
data_wi_unip_df['cmpd_seqs_pair'] = data_wi_unip_df['Substrate'] + ',' + data_wi_unip_df['Uniprot ID']
data_wi_unip_df['comp-prot-org'] = data_wi_unip_df['Substrate'] + ',' + data_wi_unip_df['Uniprot ID'] + ',' + data_wi_unip_df['Organism']

unique_compound_protein = data_wi_unip_df['cmpd_seqs_pair'].unique()
unique_compound_protein_organism = data_wi_unip_df['comp-prot-org'].unique()

#%% Saving Uniprot IDs

uniprot_ids = filtered_raw.iloc[:,6]
uniprot_ids = [x for x in uniprot_ids if '-' not in x]

def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


uniprot_ids = list(flatten(uniprot_ids))
print("len(uniprot_ids): ", len(uniprot_ids))

np.savetxt(path+'km_uniprot_ids.csv', uniprot_ids, delimiter=",", fmt='%s')






#%% Processing UniProt ID to Sequence Data from UniProt

uniprot_info_df = pd.read_csv(path+'/final_unip_seqs_dict.csv')
uniprot_info_df.rename(columns={'Entry': 'Uniprot ID'}, inplace=True)
#uniprot_to_seq = uniprot_info_df.drop(columns=['Entry', 'Reviewed', 'Entry Name', 'Protein names', 'Gene Names', 'Organism', 'EC number', 'Length'])
uniprot_to_seq_dict_0 = dict([])
unip_list = uniprot_info_df["Uniprot ID"].tolist()
seqs_list = uniprot_info_df["Sequence"].tolist()
for i, (unip, seqs) in enumerate(zip(unip_list, seqs_list)):
    uniprot_to_seq_dict_0[unip] = seqs


df_split_uniprotids = df['Uniprot ID'].str.split(',', expand=True)

data_wi_unip_df = df[df['Uniprot ID'].str.contains('-')==False]
data_wo_unip_df = df[df['Uniprot ID'].str.contains('-')==True]
data_wi_unip_df['Sequence'] = np.nan
data_wo_unip_df['Sequence'] = np.nan

# Exploding data points with multiple Uniprot IDs listed
data_wi_unip_df['Uniprot ID'] = data_wi_unip_df['Uniprot ID'].str.split(',')
data_wi_unip_df = data_wi_unip_df.explode('Uniprot ID')

# Separating data with vs without Uniprot IDs
data_wo_unip_df = df[df['Uniprot ID'].str.contains('-')==True]
data_wi_unip_df['Sequence'] = np.nan
data_wo_unip_df['Sequence'] = np.nan
    
# Quering Uniprot data for sequences (##### ?????: choose to use the first uniprot ID and left out all the rest.)
'''
for i in range(len(data_wi_unip_df)):
    print(i, "out of", len(data_wi_unip_df)) if i % 100 == 0 else 0
    ID = data_wi_unip_df.iloc[i,2].split(',')[0]
    seq_row = uniprot_to_seq[uniprot_to_seq['Uniprot ID'].str.contains(ID)]['Sequence'].tolist()
    if len(seq_row)>0:
        sequence = seq_row[0]
        data_wi_unip_df.iloc[i, -1] = sequence
        '''
data_wi_unip_unip_list = data_wi_unip_df['Uniprot ID'].tolist()

data_wi_unip_seqs_list = []
for unip in data_wi_unip_unip_list:
    if unip in uniprot_to_seq_dict_0:
        data_wi_unip_seqs_list.append(uniprot_to_seq_dict_0[unip])
    else:
        data_wi_unip_seqs_list.append(np.nan)

data_wi_unip_df["Sequence"] = data_wi_unip_seqs_list
        
#%% Separating data with found sequences to data w/out sequences mapped YET

data_wi_unip_wo_seqs_df = data_wi_unip_df[data_wi_unip_df['Sequence'].isnull()]
print("len(data_wi_unip_wo_seqs_df): ", len(data_wi_unip_wo_seqs_df)) # len(data_wi_unip_wo_seqs_df):  9580

data_wo_seqs_df = pd.concat([data_wo_unip_df, data_wi_unip_wo_seqs_df], axis=0)
data_wi_seqs_df = data_wi_unip_df[~data_wi_unip_df['Sequence'].isnull()]
print("len(data_wo_seqs_df): ", len(data_wo_seqs_df)) # len(data_wo_seqs_df):  115017
print("len(data_wi_seqs_df): ", len(data_wi_seqs_df)) # len(data_wi_seqs_df):  58566


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#    pd""b.            .g8""8q.   `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM         ,M' dP         
#   (O)  `8b         .dP'    `YM.   MM       M  P'   MM   `7   MM   `MM.  MM       M  P'   MM   `7         dP .M'    __,  
#        ,89         dM'      `MM   MM       M       MM        MM   ,M9   MM       M       MM           mmmMmmMmm   `7MM  
#      ""Yb.         MM        MM   MM       M       MM        MMmmdM9    MM       M       MM             MP dP       MM  
#         88         MM.      ,MP   MM       M       MM        MM         MM       M       MM          mmdMmmMmmm     MM  
#   (O)  .M'  ,,     `Mb.    ,dP'   YM.     ,M       MM        MM         YM.     ,M       MM           ,M' dP        MM  
#    bmmmd'   db       `"bmmd"'      `bmmmmd"'     .JMML.    .JMML.        `bmmmmd"'     .JMML.         dP ,M'      .JMML.
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 

#%% Getting compound name to SMILES dictionary

cmpd_smls_dict_df = pd.read_csv(path+'../cmpd_smls_all.csv', on_bad_lines='skip', delimiter=';')
cmpd_smls_dict_df = cmpd_smls_dict_df.drop_duplicates(['CMPD'], keep='first')
cmpd_smls_dict_df.rename(columns={'CMPD': 'Substrate', 'CMPD_SMILES': 'SMILES'}, inplace=True)
cmpd_smls_dict_df.drop(columns='Unnamed: 0', inplace=True)


#%% Getting unambigious dataset with Uniprot ID only sequences

data_wi_seqs_df_0 = data_wi_seqs_df 

data_wi_seqs_df_0['cmpd_seqs_pair'] = data_wi_seqs_df_0['Substrate'] + ',' + data_wi_seqs_df_0['Sequence']
data_wi_seqs_df_0['Km [mM]'] = pd.to_numeric(data_wi_seqs_df_0['Km [mM]'], downcast="float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD. 
data_avg_value_df = data_wi_seqs_df_0.groupby('cmpd_seqs_pair').agg(val_avg=pd.NamedAgg(column="Km [mM]",aggfunc='mean'), val_SD=pd.NamedAgg(column="Km [mM]", aggfunc='std'))
data_avg_aggre_df = data_wi_seqs_df_0.groupby('cmpd_seqs_pair').agg(lambda column: "/".join(column))
print("\n\nTaking average and standard deviation of the value, data_avg_value_df: ")
beautiful_print(data_avg_value_df)





data_wi_seqs_avg_df = pd.concat([data_avg_aggre_df, data_avg_value_df], axis=1)                                     
data_wi_seqs_avg_df['Organism'] = data_wi_seqs_avg_df['Organism'].str.split('/').str[0]
data_wi_seqs_avg_df['EC Number'] = data_wi_seqs_avg_df['EC Number'].str.split('/').str[0]
data_wi_seqs_avg_df['Substrate'] = data_wi_seqs_avg_df['Substrate'].str.split('/').str[0]
data_wi_seqs_avg_df['Sequence'] = data_wi_seqs_avg_df['Sequence'].str.split('/').str[0]
data_wi_seqs_avg_df['Uniprot ID'] = data_wi_seqs_avg_df['Uniprot ID'].str.split('/').str[0]

# Removing any datapoints where Km SD is more than or equal to 1 
data_wi_seqs_avg_df_small_sd = data_wi_seqs_avg_df.loc[(data_wi_seqs_avg_df['val_SD']<val_SD_threshold)]
data_wi_seqs_avg_df_one_value = data_wi_seqs_avg_df.loc[data_wi_seqs_avg_df['val_SD'].isnull()]

data_wi_seqs_avg_df_screened = pd.concat([data_wi_seqs_avg_df_small_sd, data_wi_seqs_avg_df_one_value], axis=0)

left_df = data_wi_seqs_avg_df_screened
right_df = cmpd_smls_dict_df

# Getting SMILES for this dataset 
#data_wi_seqs_avg_df_screened = left_df.merge(right_df, on=['Substrate'], how='left')
cmpd_list = cmpd_smls_dict_df["Substrate"].tolist()
smls_list = cmpd_smls_dict_df["SMILES"].tolist()
cmpd_smls_dict_0 = dict([])
for i, (cmpd, smls) in enumerate(zip(cmpd_list, smls_list)):
    cmpd_smls_dict_0[cmpd] = smls

data_cmpd_list = data_wi_seqs_avg_df_screened['Substrate'].tolist()

data_smls_list = []
for cmpd in data_cmpd_list:
    if cmpd in cmpd_smls_dict_0:
        data_smls_list.append(cmpd_smls_dict_0[cmpd])
    else:
        data_smls_list.append(np.nan)

data_wi_seqs_avg_df_screened["SMILES"] = data_smls_list




# Dataframe clean up and saving csv 
data_wi_seqs_avg_df_screened = data_wi_seqs_avg_df_screened.loc[~data_wi_seqs_avg_df_screened['SMILES'].isnull()]
data_wi_seqs_avg_df_screened = data_wi_seqs_avg_df_screened[data_wi_seqs_avg_df_screened['SMILES'].str.contains('None')==False]
data_wi_seqs_avg_df_screened.drop(columns='val_SD', inplace=True)
data_wi_seqs_avg_df_screened.rename(columns={'val_avg': 'Km [mM]'}, inplace=True)
data_wi_seqs_avg_df_screened = data_wi_seqs_avg_df_screened.iloc[:, [0,1, 2, 4, 3, -1, 5]]
data_wi_seqs_avg_df_screened.to_csv(path+'brenda_km_uniprot_only.csv')


seqs_list = data_wi_seqs_avg_df_screened["Sequence"].tolist()
smls_list = data_wi_seqs_avg_df_screened["SMILES"  ].tolist()
valu_list = data_wi_seqs_avg_df_screened["Km [mM]" ].tolist()

pair_list = [(smls, seqs) for smls, seqs in zip(smls_list, seqs_list)]
uniq_pair_list = list(set(pair_list))

print("\n\n_wi_seqs_avg_val_screened, len(uniq_pair_list): ", len(uniq_pair_list))





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#         ,AM             .g8""8q.   `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM         ,M' dP   
#        AVMM           .dP'    `YM.   MM       M  P'   MM   `7   MM   `MM.  MM       M  P'   MM   `7         dP .M'   
#      ,W' MM           dM'      `MM   MM       M       MM        MM   ,M9   MM       M       MM           mmmMmmMmm    pd*"*b.
#    ,W'   MM           MM        MM   MM       M       MM        MMmmdM9    MM       M       MM             MP dP     (O)   j8
#    AmmmmmMMmm         MM.      ,MP   MM       M       MM        MM         MM       M       MM          mmdMmmMmmm       ,;j9
#          MM    ,,     `Mb.    ,dP'   YM.     ,M       MM        MM         YM.     ,M       MM           ,M' dP       ,-='   
#          MM    db       `"bmmd"'      `bmmmmd"'     .JMML.    .JMML.        `bmmmmd"'     .JMML.         dP ,M'      Ammmmmmm
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 


#%% Processing EC Number to Sequence Data from KinMod

Seq_EC_Org_df = pd.read_csv(path + '../Seq_EC_Org_KINMOD.csv')
Seq_EC_Org_df.columns.values[0] = 'EC Number'
Seq_EC_Org_df.columns.values[1] = 'Sequence'
Seq_EC_Org_df.columns.values[2] = 'Organism'

orgm_list = []


for i in range(len(Seq_EC_Org_df)):
    orgn = Seq_EC_Org_df.iloc[i, -1]
    orgn = (str(orgn))
    org = orgn.split(' ')
    org = org[0:2]
    org = ' '.join(org)
    orgm_list.append(org)

Seq_EC_Org_df['Organism_Name'] = orgm_list
Seq_EC_Org_df.drop(columns=['Organism'], inplace=True)
Seq_EC_Org_df.rename(columns={"Organism_Name": "Organism"}, inplace=True)


Seq_EC_Org_grouped_df = Seq_EC_Org_df.groupby(['Organism','EC Number'])['Sequence'].apply(list)
Seq_EC_Org_grouped_df = Seq_EC_Org_grouped_df.to_frame()
Seq_EC_Org_grouped_df.reset_index(inplace=True)

# Create a dictionary that reads org and ec# and find sequence.
Seq_EC_Org_df = pd.read_csv(data_folder / ".." /"Seq_EC_Org_KINMOD.csv", 
                            index_col = None, 
                            header = 0, 
                            sep = ",")

ecno_list = Seq_EC_Org_df["EC" ].tolist()
seqs_list = Seq_EC_Org_df["SEQ"].tolist()
orgm_list = Seq_EC_Org_df["ORG"].tolist()

Seq_EC_Org_dict = dict([])
for i, (ecno, orgm) in enumerate(zip(ecno_list, orgm_list)):
    orgm = " ".join(str(orgm).split(" ")[0:2])
    if (ecno, str(orgm).lower()) not in Seq_EC_Org_dict:
        Seq_EC_Org_dict[(ecno, str(orgm).lower())] = [seqs_list[i],]
    else:
        Seq_EC_Org_dict[(ecno, str(orgm).lower())].append(seqs_list[i])

del seqs_list
del ecno_list
del orgm_list





#%% Merging the Km df with the sequence df to get a combined df

left_df = data_wo_seqs_df.drop(columns='Sequence')
right_df = Seq_EC_Org_grouped_df

data_assgn_seqs_df = left_df.merge(right_df, on=['Organism','EC Number'], how='left')

data_assgn_seqs_df = data_assgn_seqs_df.explode('Sequence') # Exploding the sequence data to ensure 1 substrate is matched to 1 sequence



'''
# Splitting the Ki values
for i in range(len(data_assgn_seqs_df)):
    print(i, "out of", len(data_assgn_seqs_df)) if i%100 == 0 else 0
    data_values = data_assgn_seqs_df['Km [mM]'].iloc[i]
    split = data_values.split(',')
    split = [float(x) for x in split]
    data_assgn_seqs_df['Km [mM]'].iloc[i] = split
    '''


data_values_list = [[float(x) for x in one_value.split(',')] for one_value in data_assgn_seqs_df["Km [mM]"].tolist()]
data_assgn_seqs_df["Km [mM]"] = data_values_list


data_assgn_seqs_df = data_assgn_seqs_df.explode('Km [mM]') # Exploding the Km values to ensure each substrate is matched to 1 Km value and 1 sequence
print("\n\n"+"="*90+"\nDataframe without sequences got values exploded, data_assgn_seqs_df: ")
beautiful_print(data_assgn_seqs_df)



no_zeros_km = data_assgn_seqs_df.sort_values(by=['Km [mM]'])
no_zeros_km.reset_index(inplace=True)
no_zeros_km.drop(columns='index', inplace=True)
no_zeros_km = no_zeros_km[no_zeros_km['Km [mM]']>0].dropna()

final_df = no_zeros_km[~no_zeros_km['Sequence'].isnull()]
data_wi_seqs_df.drop(columns='cmpd_seqs_pair', inplace=True)
final_df = pd.concat([final_df, data_wi_seqs_df], axis=0)

# Removing any duplicate rows
bool_series = final_df.duplicated(keep='first')
final_df = final_df[~bool_series]





#%% Taking Average of Km values
df = final_df
df['cmpd_seqs_pair'] = df['Substrate'] + ',' + df['Sequence']
df['Km [mM]'] = pd.to_numeric(df['Km [mM]'], downcast="float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD.
data_avg_value_df = df.groupby('cmpd_seqs_pair').agg(val_avg=pd.NamedAgg(column="Km [mM]",aggfunc='mean'), val_SD=pd.NamedAgg(column="Km [mM]", aggfunc='std'))
data_avg_aggre_df = df.groupby('cmpd_seqs_pair').agg(lambda column: "/".join(column))
data_wi_seqs_avg_df = pd.concat([data_avg_aggre_df, data_avg_value_df], axis=1)   


                                  
data_wi_seqs_avg_df['Organism'] = data_wi_seqs_avg_df['Organism'].str.split('/').str[0]
data_wi_seqs_avg_df['EC Number'] = data_wi_seqs_avg_df['EC Number'].str.split('/').str[0]
data_wi_seqs_avg_df['Substrate'] = data_wi_seqs_avg_df['Substrate'].str.split('/').str[0]
data_wi_seqs_avg_df['Sequence'] = data_wi_seqs_avg_df['Sequence'].str.split('/').str[0]
data_wi_seqs_avg_df['Uniprot ID'] = data_wi_seqs_avg_df['Uniprot ID'].str.split('/').str[0]

# Removing any entries with SD >= 1 
data_wi_seqs_avg_df_small_sd = data_wi_seqs_avg_df.loc[(data_wi_seqs_avg_df['val_SD']<val_SD_threshold)]
data_wi_seqs_avg_df_one_value = data_wi_seqs_avg_df.loc[data_wi_seqs_avg_df['val_SD'].isnull()]

data_wi_seqs_avg_df = pd.concat([data_wi_seqs_avg_df_small_sd, data_wi_seqs_avg_df_one_value], axis=0)


data_wi_seqs_avg_df.reset_index(inplace=True)
print(data_wi_seqs_avg_df['cmpd_seqs_pair'].nunique())

#%% Obtaining SMILES for Substrates from Dictionary

data_wi_seqs_avg_df_all = data_wi_seqs_avg_df.merge(cmpd_smls_dict_df, on=['Substrate'], how='left')
data_wi_seqs_avg_df_all = data_wi_seqs_avg_df_all[data_wi_seqs_avg_df_all['SMILES'].str.contains('None')==False]

data_wi_seqs_avg_df_all.drop(columns=['cmpd_seqs_pair', 'val_SD'], inplace=True)
data_wi_seqs_avg_df_all.rename(columns={'val_avg' : 'Km [mM]'}, inplace=True)

data_wi_seqs_avg_df_all = data_wi_seqs_avg_df_all.iloc[:, [0, 1, 2, 4, 3, -1, 5]]
data_wi_seqs_avg_df_all.to_csv(path+'brenda_km_pangenomic_chkpt.csv')

print("brenda_km_pangenomic_chkpt.csv generated.")



#%% Filtering the data based on Km variation across different sequences.

file_path = path + 'brenda_km_pangenomic_chkpt.csv'

df = pd.read_csv(file_path)

# Getting dataframes that search for organism, ec number, substrate entries that are the same and search the Km, Km SD, and number of Uniprot IDs
data_assgn_seqs_df = df.groupby(['Organism', 'EC Number', 'Substrate', 'Uniprot ID']).agg(val_avg=pd.NamedAgg(column="Km [mM]",aggfunc='mean'), val_SD=pd.NamedAgg(column="Km [mM]", aggfunc='std'))               
data_assgn_seqs_df_2 = df.groupby(['Organism', 'EC Number', 'Substrate']).agg(val_avg=pd.NamedAgg(column="Km [mM]",aggfunc='mean'), val_SD=pd.NamedAgg(column="Km [mM]", aggfunc='std'), UniProt_IDs=pd.NamedAgg(column='Uniprot ID', aggfunc='count'))
data_assgn_seqs_df_2.reset_index(inplace=True)

# Creating a dictionary df that has Organism, EC number, and substrate triplicate entries that don't have large variation between Km values and different sequences
data_assgn_seqs_df_2 = data_assgn_seqs_df_2[data_assgn_seqs_df_2['val_SD']<val_SD_threshold]


# Merging datasets and filtering the pangenomic dataset for the Km SD condition above
df2 = pd.merge(df, data_assgn_seqs_df_2, on=['Organism', 'EC Number', 'Substrate'], how='left')    

df = df2.dropna(subset=['val_avg'])
df = df.drop(columns=['UniProt_IDs', 'val_avg', 'val_SD'])

df.to_csv(path+'brenda_km_pangenomic.csv')

print("brenda_km_pangenomic.csv generated.")



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




