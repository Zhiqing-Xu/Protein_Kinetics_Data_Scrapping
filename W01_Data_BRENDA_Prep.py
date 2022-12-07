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

data_folder_name    = ["Ki_BRENDA"        , "KM_3_BRENDA"      , "kcat_BRENDA"        , "kcat_KM_BRENDA"        ][2]
data_file_name      = ["brenda_Ki_raw.csv", "brenda_KM_raw.csv", "brenda_kcat_raw.csv", "brenda_kcat_KM_raw.csv"][2]
output_file_name    = ["Ki_BRENDA.csv"    , "KM_BRENDA.csv"    , "kcat_BRENDA.csv"    , "kcat_KM_BRENDA.csv"    ][2]
data_name           = ["Ki"               , "Km"               , "kcat"               , "kcat_KM"               ][2]
data_name_0         = ["max_Ki"           , "max_Km"           , "max_kcat"           , "kcat_KM_0"             ][2]

data_folder      = Path("X_DataProcessing/X00_enzyme_datasets/" + data_folder_name)
data_file        = data_file_name

output_folder    = Path("X_DataProcessing/X00_enzyme_datasets_processed/")
output_file      = output_file_name

output_1_tail    = '_core_unip_only.csv'
output_2_tail    = '_scrnd_all.csv'
output_3_tail    = '_fine_scrnd_all.csv'

if data_name in ["Ki", "Km", "kcat", ]:
    val_SD_threshold = 1
    
if data_name in ["kcat_KM", ]:
    val_SD_threshold = 10

################################################################################################################### KM: 62530, 28466, 19597
#   `7MMM.     ,MMF'      db      `7MMF'`7MN.   `7MF'     `7MM"""YMM `7MMF'`7MMF'      `7MM"""YMM  
#     MMMb    dPMM       ;MM:       MM    MMN.    M         MM    `7   MM    MM          MM    `7  
#     M YM   ,M MM      ,V^MM.      MM    M YMb   M         MM   d     MM    MM          MM   d    
#     M  Mb  M' MM     ,M  `MM      MM    M  `MN. M         MM""MM     MM    MM          MMmmMM    
#     M  YM.P'  MM     AbmmmqMA     MM    M   `MM.M         MM   Y     MM    MM      ,   MM   Y  , 
#     M  `YM'   MM    A'     VML    MM    M     YMM         MM         MM    MM     ,M   MM     ,M 
#   .JML. `'  .JMML..AMA.   .AMMA..JMML..JML.    YM       .JMML.     .JMML..JMMmmmmMMM .JMMmmmmMMM 
###################################################################################################################

# Read the main data file.
raw_df_0 = pd.read_csv(filepath_or_buffer   =   data_folder / data_file, 
                       on_bad_lines         =   'skip', 
                       index_col            =   None, 
                       #names                =   ["EC", data_name, data_name_0, "cmpd", "condition", "organism", "uniprot_id", "publication", "other"], 
                       header               =   None, 
                       sep                  =   '	', 
                       encoding             =   'cp1252')

raw_df_0.rename(columns = {0 : 'EC'          , 
                           1 : data_name     ,
                           2 : data_name_0   ,
                           3 : 'cmpd'        ,
                           4 : 'condition'   ,
                           5 : 'organism'    ,
                           6 : 'uniprot_id'  ,
                           7 : 'publication' ,
                           8 : 'other'       ,
                           }, 
                inplace = True)

print("\n\n"+"="*90+"\n#0.0 Raw data read from the csv file, raw_df_0: ")
beautiful_print(raw_df_0)
print("len(raw_df_0): ", len(raw_df_0))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Filter Something (optional).
raw_df_0_1 = copy.deepcopy(raw_df_0)

raw_df_0_1 = raw_df_0_1[ raw_df_0_1[data_name].str.contains("additional information") == False ] #1
raw_df_0_1 = raw_df_0_1[ raw_df_0_1[data_name].str.strip() != "-" ] #1
#print(len(raw_df_0_1)) # [out]: 162216.





# Use log values for later SD<1.0 filtering.
#raw_df_0_1[data_name] = pd.to_numeric(raw_df_0_1[data_name], downcast="float")
#raw_df_0_1[data_name] = raw_df_0_1[data_name].apply(lambda x: np.log10(x))


raw_df_0_2 = raw_df_0_1[ raw_df_0_1["cmpd"].str.contains("additional information") == False] #2
#print(len(raw_df_0_2)) # [out]: 162208.

"""
There is not a row that contains 'additional information', yet some rows are still removed.
This is because there are rows with too few fields in the dataset.
The function 'read_csv()' read those rows and the value of those fields as NaN.
The #2 filtering remove those rows as the value of those fields are requested.
# [in]: print(raw_df_0[["cmpd",]].iloc[38737]) [out]: NaN
"""

# Keep using raw_df_0, remove those rows with a max_[data_name] value (optional).
raw_df_0_3  = raw_df_0_2[ raw_df_0_2[data_name_0].str.contains('-') == True ] #3
#print("len(raw_df_0): ", len(raw_df_0_3))

raw_df_0 = raw_df_0_2

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Remove useless columnms.
try:
    raw_df_0 = raw_df_0.drop(columns=[data_name_0, "condition", "publication", "other"])
except:
    raw_df_0 = raw_df_0.drop(columns=[data_name_0, "condition", "publication", ])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Change the order of the columns.
raw_df_0 = raw_df_0.iloc[:, [-2, 0, -1, 2, 1]]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Preprocessing.
raw_df_0['uniprot_id'] = raw_df_0['uniprot_id'].str.replace('and', ',')
raw_df_0['uniprot_id'] = raw_df_0['uniprot_id'].str.replace('AND', ',')
raw_df_0['uniprot_id'] = raw_df_0['uniprot_id'].str.replace(' ', '')

print("\n\n"+"="*90+"\n#0.1 Prescreen and simple clean, raw_df_0: ")
beautiful_print(raw_df_0)
print("len(raw_df_0): ", len(raw_df_0))



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#                `7MMF'   `7MF'`7MN.   `7MF'`7MMF'`7MM"""Mq. `7MM"""Mq.    .g8""8q.   MMP""MM""YMM     `7MMF'`7MM"""Yb.           
#    __,           MM       M    MMN.    M    MM    MM   `MM.  MM   `MM. .dP'    `YM. P'   MM   `7       MM    MM    `Yb.         
#   `7MM           MM       M    M YMb   M    MM    MM   ,M9   MM   ,M9  dM'      `MM      MM            MM    MM     `Mb ,pP"Ybd 
#     MM           MM       M    M  `MN. M    MM    MMmmdM9    MMmmdM9   MM        MM      MM            MM    MM      MM 8I   `" 
#     MM           MM       M    M   `MM.M    MM    MM         MM  YM.   MM.      ,MP      MM            MM    MM     ,MP `YMMMa. 
#     MM  ,,       YM.     ,M    M     YMM    MM    MM         MM   `Mb. `Mb.    ,dP'      MM            MM    MM    ,dP' L.   I8 
#   .JMML.db        `bmmmmd"'  .JML.    YM  .JMML..JMML.     .JMML. .JMM.  `"bmmd"'      .JMML.        .JMML..JMMmmmdP'   M9mmmP' 
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 

#====================================================================================================#
#%% Looking at dataset with UniProt Available Only.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Get a dataframe that contains data with uniprot ids only. (Rewritten below)
'''
data_wi_unip_df = copy.deepcopy(raw_df_0)
data_wi_unip_df = data_wi_unip_df[ data_wi_unip_df['uniprot_id'].str.contains('-')                      == False ]
data_wi_unip_df = data_wi_unip_df[ data_wi_unip_df[data_name]   .str.contains("additional information") == False ] 
'''

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Extend the dataset for some test. (Skipped)
'''
# Get combinations of compounds and uniprot ids.
data_wi_unip_df['cmpd_unip_pair'] = data_wi_unip_df['cmpd'] + ',' + data_wi_unip_df['uniprot_id']
data_wi_unip_df['cmpd_unip_orgm'] = data_wi_unip_df['cmpd'] + ',' + data_wi_unip_df['uniprot_id'] + ',' + str(data_wi_unip_df['organism'])

# Get arrays of unique combinations of compounds and uniprot ids.
uniq_cmpd_unip_pair_array = data_wi_unip_df['cmpd_unip_pair'].unique()
uniq_cmpd_unip_orgm_array = data_wi_unip_df['cmpd_unip_orgm'].unique()
'''

#====================================================================================================#
#%%  Rana's extra work of getting a list of UniProt IDs.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Save the dataset that contains data with uniprot ids only. (Skipped)
'''
# Get a list of UniProt IDs.
uniprot_ids_col = raw_df_0['uniprot_id']
uniprot_ids_col = uniprot_ids_col.fillna("")
print("Number of rows of raw data, len(uniprot_ids_col): ", len(uniprot_ids_col))

uniprot_ids_list = [one_uniprot_id for one_uniprot_id in uniprot_ids_col if '-' not in one_uniprot_id]
print("Number of rows of data w/ Uniprot IDs, len(uniprot_ids_list): ", len(uniprot_ids_list))


# Recursively flatten nested list/tuple.
def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


# ? Why did Rana flatten uniprot_ids_list? 
# ? It doesnt seem to contain nested tuples/lists.
uniprot_ids_list_2 = list(flatten(uniprot_ids_list)) 
print("len(uniprot_ids_list): ", len(uniprot_ids_list_2))
print(uniprot_ids_list_2 == uniprot_ids_list)


# Saving Uniprot IDs
np.savetxt(data_folder / (data_name + '_data_wi_unip.csv'), uniprot_ids_list, delimiter = ",", fmt = '%s')
'''

#====================================================================================================#
#%% Separating data with vs without Uniprot IDs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Check UniProt IDs for each record in the dataset and split the data accordingly.
data_wi_unip_df = data_wo_unip_df = copy.deepcopy(raw_df_0)
data_wi_unip_df = data_wi_unip_df[ data_wi_unip_df['uniprot_id'].str.contains('-') == False ]
data_wo_unip_df = data_wo_unip_df[ data_wo_unip_df['uniprot_id'].str.contains('-') == True  ]


'''
# Process in a different way if there are "additional information" in the uniprot_id column.
# It turns out that there is NO "additional information" in the uniprot_id column.
data_wi_unip_df = data_wo_unip_df = copy.deepcopy(raw_df_0)

# Data with UniProt IDs. (First, remove "-". Then get "addi info" for data_wo_unip_df. Finally remove "addi info".)
data_wi_unip_df   = data_wi_unip_df[ data_wi_unip_df['uniprot_id'].str.contains('-')                      == False ]

data_wo_unip_df_1 = data_wi_unip_df[ data_wi_unip_df[data_name]   .str.contains("additional information") == True  ]
data_wi_unip_df   = data_wi_unip_df[ data_wi_unip_df[data_name]   .str.contains("additional information") == False ] 

# Data without UniProt IDs.
data_wo_unip_df_2 = data_wo_unip_df[ data_wo_unip_df['uniprot_id'].str.contains('-') == True ]
print("len(data_wo_unip_df_1): ", len(data_wo_unip_df_1), ", those with `additional information`. ")
print("len(data_wo_unip_df_2): ", len(data_wo_unip_df_2), ", those with `-`. ")

data_wo_unip_df   = pd.concat([data_wo_unip_df_1, data_wo_unip_df_2], axis=0)
print("len(data_wo_unip_df): ", len(data_wo_unip_df), ", combination of two. ")
'''

# Exploding data points with multiple Uniprot IDs listed
print("len(data_wi_unip_df) before/after explode: ", len(data_wi_unip_df), end = ", ")
data_wi_unip_df['uniprot_id'] = data_wi_unip_df['uniprot_id'].str.split(',')
data_wi_unip_df = data_wi_unip_df.explode('uniprot_id')
print(len(data_wi_unip_df))

# Create the Sequence column.
data_wi_unip_df['Sequence'] = np.nan
data_wo_unip_df['Sequence'] = np.nan

# Checkpoint. 
print("len(data_wo_unip_df): ", len(data_wo_unip_df))
print("len(raw_df_0): ", len(raw_df_0))
print("Confirmed that: *** len(raw_df_0) == len(data_wo_unip_df) + len(data_wi_unip_df) *** (before explosion)")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#                   .M"""bgd `7MM"""YMM    .g8""8q.   `7MMF'   `7MF'`7MM"""YMM  `7MN.   `7MF'  .g8"""bgd `7MM"""YMM  
#                  ,MI    "Y   MM    `7  .dP'    `YM.   MM       M    MM    `7    MMN.    M  .dP'     `M   MM    `7  
#   pd*"*b.        `MMb.       MM   d    dM'      `MM   MM       M    MM   d      M YMb   M  dM'       `   MM   d    
#  (O)   j8          `YMMNq.   MMmmMM    MM        MM   MM       M    MMmmMM      M  `MN. M  MM            MMmmMM    
#      ,;j9        .     `MM   MM   Y  , MM.      ,MP   MM       M    MM   Y  ,   M   `MM.M  MM.           MM   Y  , 
#   ,-='    ,,     Mb     dM   MM     ,M `Mb.    ,dP'   YM.     ,M    MM     ,M   M     YMM  `Mb.     ,'   MM     ,M 
#  Ammmmmmm db     P"Ybmmd"  .JMMmmmmMMM   `"bmmd"'      `bmmmmd"'  .JMMmmmmMMM .JML.    YM    `"bmmmd'  .JMMmmmmMMM 
#                                           MMb                                                                   
#                                            `bood'                                                               
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
# Processing UniProt IDs and Sequences.

###################################################################################################################
#                .g8"""bgd  `7MM"""YMM  MMP""MM""YMM     `7MM"""Yb.   `7MMF'  .g8"""bgd MMP""MM""YMM              #
#              .dP'     `M    MM    `7  P'   MM   `7       MM    `Yb.   MM  .dP'     `M P'   MM   `7              #
#              dM'       `    MM   d         MM            MM     `Mb   MM  dM'       `      MM                   #
#              MM             MMmmMM         MM            MM      MM   MM  MM               MM                   #
#              MM.    `7MMF'  MM   Y  ,      MM            MM     ,MP   MM  MM.              MM                   #
#              `Mb.     MM    MM     ,M      MM            MM    ,dP'   MM  `Mb.     ,'      MM                   #
#                `"bmmmdPY  .JMMmmmmMMM    .JMML.        .JMMmmmdP'   .JMML.  `"bmmmd'     .JMML.                 #
###################################################################################################################
# Get a dictionary in order to convert UniProt IDs to sequences.
#====================================================================================================#
# Read a table contains a large number of AA_seqs info.
# Return a dictionary of UniProt IDs -> sequences.
def look_up_seqs_dict(dictionary         =  "uniprot_annotation_3_201_300.tsv" ,
                      search_unip_list   =  []                         ,
                      ):

    #--------------------------------------------------#
    # Initialize dicts.
    # Dict to be output, contains the uniprot_ids for those in the dataset only.
    unip_seqs_dict_x = dict([]) 
    # List to be output, contains the uniprot_ids that are NOT found in the dict file.
    missing_unip_list = []
    # Dict that used to store UniProt IDs -> sequences through reading the dict file.
    all_unip_seqs_dict = dict([])
    #--------------------------------------------------#
    # Read the dict file.
    begin_time = time.time()
    
    if dictionary.find(".tsv") != -1:
        uniprot_info_all_df = pd.read_csv(filepath_or_buffer   =   data_folder / ".." /  dictionary,
                                          header               =   0       , 
                                          sep                  =   '\t'    , 
                                          encoding             =   'cp1252', 
                                          low_memory           =   False   , )
    elif dictionary.find(".csv") != -1:
        uniprot_info_all_df = pd.read_csv(filepath_or_buffer   =   data_folder / ".." /  dictionary,
                                          header               =   0        , 
                                          sep                  =   ','      , 
                                          encoding             =   'cp1252' , 
                                          low_memory           =   True    , )
    else:
        print("File Name Error. Expect a csv or tsv file.")

    seqs_list = uniprot_info_all_df["Sequence"].tolist()
    unip_list = uniprot_info_all_df["Entry"].tolist()

    for i, (seqs, unip) in enumerate(zip(seqs_list, unip_list)):
        all_unip_seqs_dict[unip] = seqs
    del seqs_list, unip_list
    #--------------------------------------------------#
    print("len(all_unip_seqs_dict): ", len(all_unip_seqs_dict), ". # of mapping: UniProt IDs -> sequences. ")
    print("read_csv takes ", time.time() - begin_time, "seconds.")

    #--------------------------------------------------#
    # Use all_unip_seqs_dict to look up the list of uniprot_ids to be searched.
    # Get a dictionary for the list of uniprot_ids to be searched.
    for unip in search_unip_list:
        if unip in all_unip_seqs_dict:
            unip_seqs_dict_x[unip] = all_unip_seqs_dict[unip]
        else:
            missing_unip_list.append(unip)

    return unip_seqs_dict_x, missing_unip_list, len(missing_unip_list)


#====================================================================================================#
# Look up dictionaries. Generate final_unip_seqs_dict.csv if the file doesnt exist.
# After the final_unip_seqs_dict.csv is generated for the first time, use it directly.

# Full unip -> seqs dictionary downloaded from UniProtKB is very large. 
# Once a unip -> seqs dict is generated, do not check those large files.
# This is totally different from processing cmpd_name -> smiles, 
#  where dictionaries are NOT that large.

if not Path.exists(data_folder / "final_unip_seqs_dict.csv"):
    # If there is NOT a final dictionary, get one.
    final_unip_seqs_dict = dict([])
    # Get a list of unique uniprot ids.
    unique_unip_list = list(set(data_wi_unip_df['uniprot_id'].tolist()))
    print("len(unique_unip_list): (Number of uniprot IDs in dataset)", len(unique_unip_list))
    # Start with a small dict file. 
    final_unip_seqs_dict, _, _ = look_up_seqs_dict(dictionary = "uniprot_AC_SQ_2_1.csv"   , search_unip_list = unique_unip_list )
    #--------------------------------------------------#
    # The following dictionaries are more than 25 GB combined. Avoid this step as long as final_unip_seqs_dict.csv is generated.
    try:
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_5.tsv"         , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_4.tsv"         , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_001_200.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_201_300.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_301_400.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_401_600.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_601_800.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_801_inf.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
    except:
        pass
    #--------------------------------------------------#
    # Collect those uniprot ids that are left out by the dictionaries.
    uniprot_IDs_missing_list = []
    for one_unip in unique_unip_list:
        if one_unip not in final_unip_seqs_dict.keys():
            uniprot_IDs_missing_list.append(one_unip)
    print("Number of uniprot IDs not found: ", len(uniprot_IDs_missing_list))

    #--------------------------------------------------#
    # Search for sequences directly in UniProt DB.
    for i, one_unip in enumerate(uniprot_IDs_missing_list):
        print(i, "out of", len(uniprot_IDs_missing_list)) if i % 100 == 0 else 0
        try:
            xml_url = 'https://rest.uniprot.org/uniprotkb/' + str(one_unip) + '.xml'
            xml_response = requests.get(xml_url)
            xml_data = xmltodict.parse(xml_response.content)
            seqs = xml_data["uniprot"]["entry"]["sequence"]['#text']
            final_unip_seqs_dict[one_unip] = seqs
        except:
            print("Unidentified UniProt ID: ", one_unip)

    #--------------------------------------------------#
    # Finally get the final_unip_seqs_dict and save it to a csv file.
    final_unip_seqs_dict_df_data = {'Entry': list(final_unip_seqs_dict.keys()), 'Sequence': list(final_unip_seqs_dict.values())}
    final_unip_seqs_dict_df = pd.DataFrame.from_dict(final_unip_seqs_dict_df_data)
    final_unip_seqs_dict_df.to_csv( data_folder / 'final_unip_seqs_dict.csv' )

#====================================================================================================#
# Use the saved final_unip_seqs_dict.csv to look up seqs directly if there is one.
else: 
    #--------------------------------------------------#
    # Read final_unip_seqs_dict.csv
    final_unip_seqs_dict_df = pd.read_csv(data_folder / 'final_unip_seqs_dict.csv')

    seqs_list = final_unip_seqs_dict_df["Sequence"].tolist()
    unip_list = final_unip_seqs_dict_df["Entry"].tolist()
    final_unip_seqs_dict = dict([])

    for i, (seqs, unip) in enumerate(zip(seqs_list, unip_list)):
        final_unip_seqs_dict[unip] = seqs

    del seqs_list, unip_list, final_unip_seqs_dict_df
    #print(final_unip_seqs_dict)

print("Done getting one final_unip_seqs_dict AND final_unip_seqs_dict_df.")



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#   pd""b.             .g8"""bgd   .g8""8q.   `7MMM.     ,MMF'`7MM"""Mq.   .g8""8q.   `7MMF'   `7MF'`7MN.   `7MF'`7MM"""Yb.   
#  (O)  `8b          .dP'     `M .dP'    `YM.   MMMb    dPMM    MM   `MM..dP'    `YM.   MM       M    MMN.    M    MM    `Yb. 
#       ,89          dM'       ` dM'      `MM   M YM   ,M MM    MM   ,M9 dM'      `MM   MM       M    M YMb   M    MM     `Mb 
#     ""Yb.          MM          MM        MM   M  Mb  M' MM    MMmmdM9  MM        MM   MM       M    M  `MN. M    MM      MM 
#        88          MM.         MM.      ,MP   M  YM.P'  MM    MM       MM.      ,MP   MM       M    M   `MM.M    MM     ,MP 
#  (O)  .M'  ,,      `Mb.     ,' `Mb.    ,dP'   M  `YM'   MM    MM       `Mb.    ,dP'   YM.     ,M    M     YMM    MM    ,dP' 
#   bmmmd'   db        `"bmmmd'    `"bmmd"'   .JML. `'  .JMML..JMML.       `"bmmd"'      `bmmmmd"'  .JML.    YM  .JMMmmmdP'   
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
# Processing Compounds
# This is totally different from UniProt IDs -> sequences.
# Since the compounds -> SMILES dictionaries are relatively small.
# Import all dictionaries and get a comprehensive one.
# Use the combined dictionary to look up SMILES.

###################################################################################################################
#                .g8"""bgd  `7MM"""YMM  MMP""MM""YMM     `7MM"""Yb.   `7MMF'  .g8"""bgd MMP""MM""YMM              #
#              .dP'     `M    MM    `7  P'   MM   `7       MM    `Yb.   MM  .dP'     `M P'   MM   `7              #
#              dM'       `    MM   d         MM            MM     `Mb   MM  dM'       `      MM                   #
#              MM             MMmmMM         MM            MM      MM   MM  MM               MM                   #
#              MM.    `7MMF'  MM   Y  ,      MM            MM     ,MP   MM  MM.              MM                   #
#              `Mb.     MM    MM     ,M      MM            MM    ,dP'   MM  `Mb.     ,'      MM                   #
#                `"bmmmdPY  .JMMmmmmMMM    .JMML.        .JMMmmmdP'   .JMML.  `"bmmmd'     .JMML.                 #
###################################################################################################################
# Read a table contains a large number of AA_seqs info.
def look_up_cmpd_dict(additional_dict_list  =  [
                                                "cmpd_smls_pubchem.csv" , 
                                                "cmpd_smls_cactus.csv"  , 
                                                #"cmpd_smls_all.csv"     , 
                                               ] ,
                      search_cmpd_list      =  [] ,
                      data_folder           =  data_folder,
                      ):

    # Import dictionaries previously extracted in a separate program.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Dictionary #1: KEGG Database (Compound Name --> KEGG ID --> SMILES)
    # Get dict for reading ALL_KEGG_DATA_DICT. (Compound Name --> KEGG ID)
    folder = Path("./CMPD_RXN_tools/KEGGScrapSavings")
    filename = "CPD02_ALL_KEGG_DATA_DICT_bk.p"
    ALL_KEGG_DATA_DICT = pickle.load(open( folder / filename, "rb"))

    KEGG_name_id_dict = dict([])
    for one_KEGG_id in ALL_KEGG_DATA_DICT:
        if "NAME" in ALL_KEGG_DATA_DICT[one_KEGG_id]:
            for one_name in ALL_KEGG_DATA_DICT[one_KEGG_id]["NAME"]:
                KEGG_name_id_dict[one_name.lower()] = one_KEGG_id
    #--------------------------------------------------#
    # Get dict for reading KEGG compound IDs. (KEGG ID --> SMILES)
    KEGG_id_SMILES_by_mol_df = pd.read_csv("./CMPD_RXN_tools/KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_mol_file_no_uniq.csv", header = 0)

    KEGGid_list     = KEGG_id_SMILES_by_mol_df["KEGGid"].tolist()
    SMILES_by_mol   = KEGG_id_SMILES_by_mol_df["SMILES_by_mol"].tolist()

    KEGG_id_SMILES_by_mol_dict = dict([])
    for i, id in enumerate(KEGGid_list):
        KEGG_id_SMILES_by_mol_dict[id] = SMILES_by_mol[i]

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Dictionary #2: MetaNetX Database (Compound Name --> SMILES)
    folder = Path("./CMPD_RXN_tools/MNXScrapSavings")
    filename = "CPD01_nme_smiles_MNXid_dict.p"
    nme_smiles_MNXid_dict = pickle.load(open( folder / filename, "rb"))

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Dictionary #3: Import addtional dictionary.
    additional_cmpd_smls_dict = dict([])
    for one_dict_nme in additional_dict_list:
        additional_cmpd_smls_df = pd.read_csv(data_folder / ".." / one_dict_nme,
                                              encoding = "ISO-8859-1",
                                              index_col = None, 
                                              header = 0, 
                                              sep = ';',
                                              engine = 'python')

        cmpd_name_list_additional = additional_cmpd_smls_df["CMPD"].tolist()
        cmpd_smls_list_additional = additional_cmpd_smls_df["CMPD_SMILES"].tolist()

        for i, one_cmpd_name in enumerate(cmpd_name_list_additional):
            additional_cmpd_smls_dict[str(one_cmpd_name).lower()] = cmpd_smls_list_additional[i]

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Imported:
    # pubchem_cmpd_smls_dict
    # cactus_cmpd_smls_dict
    # nme_smiles_MNXid_dict
    # KEGG_name_id_dict, KEGG_id_SMILES_by_mol_dict
    # additional_cmpd_smls_dict

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Look up dictionaries in a certain order.
    # Output a final dictionary for the list of compound names input.
    cmpd_smls_dict_x = dict([])
    missing_cmpd_list = []

    for i, one_cmpd in enumerate(search_cmpd_list):
        append_bool = False
        one_cmpd = str(one_cmpd)
        one_cmpd = one_cmpd.replace(";", "")
        one_cmpd = one_cmpd.lower()

        if (not append_bool) and one_cmpd in additional_cmpd_smls_dict:
            cmpd_smls = additional_cmpd_smls_dict[one_cmpd]
            if cmpd_smls not in ["NA", "", "nan", "None"]:
                cmpd_smls_dict_x[one_cmpd] = cmpd_smls
                append_bool = True

        
        if (not append_bool) and one_cmpd in nme_smiles_MNXid_dict:
            cmpd_smls = nme_smiles_MNXid_dict[one_cmpd][0]
            if cmpd_smls not in ["NA", "", "nan", "None"]:
                cmpd_smls_dict_x[one_cmpd] = cmpd_smls
                append_bool = True

        if (not append_bool) and one_cmpd in KEGG_name_id_dict:
            if KEGG_name_id_dict[one_cmpd] not in ["NA", "", "nan", "None"] and KEGG_name_id_dict[one_cmpd] in KEGG_id_SMILES_by_mol_dict:
                if KEGG_id_SMILES_by_mol_dict[KEGG_name_id_dict[one_cmpd]] != "None":
                    cmpd_smls_dict_x[one_cmpd] = KEGG_id_SMILES_by_mol_dict[KEGG_name_id_dict[one_cmpd]]
                    append_bool = True
                            
        if (not append_bool):
            cmpd_smls_dict_x[one_cmpd] = "None"
            missing_cmpd_list.append(one_cmpd)
            #print("cmpd:"+one_cmpd+", smls: None")
            
    return cmpd_smls_dict_x, missing_cmpd_list, len(missing_cmpd_list)


###################################################################################################################
#                 `7MMF'        .g8""8q.     .g8""8q.   `7MMF' `YMM'     `7MMF'   `7MF'`7MM"""Mq.                 #
#                   MM        .dP'    `YM. .dP'    `YM.   MM   .M'         MM       M    MM   `MM.                #
#                   MM        dM'      `MM dM'      `MM   MM .d"           MM       M    MM   ,M9                 #
#                   MM        MM        MM MM        MM   MMMMM.           MM       M    MMmmdM9                  #
#                   MM      , MM.      ,MP MM.      ,MP   MM  VMA          MM       M    MM                       #
#                   MM     ,M `Mb.    ,dP' `Mb.    ,dP'   MM   `MM.        YM.     ,M    MM                       #
#                 .JMMmmmmMMM   `"bmmd"'     `"bmmd"'   .JMML.   MMb.       `bmmmmd"'  .JMML.                     #
###################################################################################################################
# Getting cmpd_name -> SMILES dictionary.
#====================================================================================================#
#%% A list of all compounds without duplicates.
uniq_cmpd_list = list(raw_df_0['cmpd'].unique())
print("Number of different compounds in the set: ", len(uniq_cmpd_list))


#%% Getting compound name to SMILES dictionary
cmpd_smls_dict, missing_cmpd_list, missing_cmpd_list_len = \
    look_up_cmpd_dict(additional_dict_list  =  ["cmpd_smls_cactus.csv"  ,
                                                "cmpd_smls_pubchem.csv" , 
                                                #"cmpd_smls_all.csv"     , 
                                               ]                        ,
                        search_cmpd_list    =  uniq_cmpd_list           ,
                        data_folder         =  data_folder              ,
                        )

print("Number of compounds with smiles found: ", len(cmpd_smls_dict) - missing_cmpd_list_len)
print("Number of compounds with no smiles found: ", missing_cmpd_list_len)

#====================================================================================================#
#%% Get a dataframe with compound names and SMILES.
cmpd_smls_dict_df = pd.DataFrame.from_dict({"CMPD": list(cmpd_smls_dict.keys()), "CMPD_SMILES": cmpd_smls_dict.values()})
cmpd_smls_dict_df = cmpd_smls_dict_df.drop_duplicates(['CMPD'], keep='first')
#cmpd_smls_dict_df.rename(columns={'CMPD': 'Substrate', 'CMPD_SMILES': 'SMILES'}, inplace=True)



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#         ,AM              .g8""8q.   `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM         ,M' dP         
#        AVMM            .dP'    `YM.   MM       M  P'   MM   `7   MM   `MM.  MM       M  P'   MM   `7         dP .M'    __,  
#      ,W' MM            dM'      `MM   MM       M       MM        MM   ,M9   MM       M       MM           mmmMmmMmm   `7MM  
#    ,W'   MM            MM        MM   MM       M       MM        MMmmdM9    MM       M       MM             MP dP       MM  
#    AmmmmmMMmm          MM.      ,MP   MM       M       MM        MM         MM       M       MM          mmdMmmMmmm     MM  
#          MM     ,,     `Mb.    ,dP'   YM.     ,M       MM        MM         YM.     ,M       MM           ,M' dP        MM  
#          MM     db       `"bmmd"'      `bmmmmd"'     .JMML.    .JMML.        `bmmmmd"'     .JMML.         dP ,M'      .JMML.
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
# Process Data with UniProt IDs. 
###################################################################################################################
#         `7MM           mm                                    db                              db                 #
#           MM           MM                                                                                       #
#      ,M""bMM   ,6"Yb.mmMMmm  ,6"Yb.      `7M'    ,A    `MF'`7MM      `7MM  `7MM `7MMpMMMb. `7MM `7MMpdMAo.      #
#    ,AP    MM  8)   MM  MM   8)   MM        VA   ,VAA   ,V    MM        MM    MM   MM    MM   MM   MM   `Wb      #
#    8MI    MM   ,pm9MM  MM    ,pm9MM         VA ,V  VA ,V     MM        MM    MM   MM    MM   MM   MM    M8      #
#    `Mb    MM  8M   MM  MM   8M   MM          VVV    VVV      MM        MM    MM   MM    MM   MM   MM   ,AP      #
#     `Wbmd"MML.`Moo9^Yo.`Mbmo`Moo9^Yo.         W      W     .JMML.      `Mbod"YML.JMML  JMML.JMML. MMbmmd'       #
#                                                                                                   MM            #
#                                                                                                 .JMML.          #
###################################################################################################################
#%% Look up seqs, and add sequence column, (for those without uniprot ids, assign a np.nan.)
data_wi_unip_unip_list = data_wi_unip_df['uniprot_id'].tolist()

data_wi_unip_seqs_list = []
for unip in data_wi_unip_unip_list:
    if unip in final_unip_seqs_dict:
        data_wi_unip_seqs_list.append(final_unip_seqs_dict[unip])
    else:
        data_wi_unip_seqs_list.append(np.nan)

data_wi_unip_df["Sequence"] = data_wi_unip_seqs_list

print("\n\n"+"="*90+"\n#4.0 Dataframe after adding sequence column, (only for those rows with a uniprot ID,) data_wi_unip_df: ")
beautiful_print(data_wi_unip_df)

# A relatively slower version.
'''
for i in range(len(data_wi_unip_df)):
    print(i, "out of", len(data_wi_unip_df)) if i % 1000 == 0 else 0
    unip = data_wi_unip_df['uniprot_id'].iloc[i]
    if unip in final_unip_seqs_dict:
        seqs = final_unip_seqs_dict[unip]
        data_wi_unip_df.iloc[i, -1] = seqs
'''

#====================================================================================================#
#%% Separate data with found sequences and data w/out sequences mapped yet.

# Get a dataframe of those data points with uniprot IDs but without a sequence found through searching UniProt DB.
data_wi_unip_wo_seqs_df = data_wi_unip_df[data_wi_unip_df['Sequence'].isnull()]
print("len(data_wi_unip_wo_seqs_df): ", len(data_wi_unip_wo_seqs_df), "# data points with uniprot IDs but without a sequence found") 
# len(data_wi_unip_wo_seqs_df): 217

# Get data_wi_seqs_df & data_wo_seqs_df.
data_wo_seqs_df = pd.concat([data_wo_unip_df, data_wi_unip_wo_seqs_df], axis=0)
data_wi_seqs_df = data_wi_unip_df[~data_wi_unip_df['Sequence'].isnull()]

# After processing uniprot_ids, redefine data_wi_seqs_df and data_wo_seqs_df.
print("After processing uniprot_ids, redefine data_wi_seqs_df and data_wo_seqs_df.")
print("len(data_wi_seqs_df): ", len(data_wi_seqs_df)) # len(data_wi_seqs_df): 68316
print("len(data_wo_seqs_df): ", len(data_wo_seqs_df)) # len(data_wo_seqs_df): 107420


# Some redundant code.
'''
#%% Getting unambigious dataset with Uniprot ID only sequences
data_wi_seqs_df_0 = copy.deepcopy(data_wi_seqs_df)

data_wi_seqs_df_0['cmpd_seqs_pair'] = data_wi_seqs_df_0['cmpd'] + ',' + data_wi_seqs_df_0['Sequence']
data_wi_seqs_df_0[data_name] = pd.to_numeric(data_wi_seqs_df_0[data_name], downcast="float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD. 
data_avg_value_df = data_wi_seqs_df_0.groupby('cmpd_seqs_pair').agg(
    val_avg = pd.NamedAgg(column = data_name, aggfunc = 'mean') , 
    val_SD  = pd.NamedAgg(column = data_name, aggfunc = 'std' ) ,
                                                                   )
#data_avg_aggre_df = data_wi_seqs_df_0.groupby('cmpd_seqs_pair').agg(lambda column: "/".join(column))

print("\n\n"+"="*90+"\nTaking average and standard deviation of the value, data_avg_value_df: ")
beautiful_print(data_avg_value_df)
'''

#====================================================================================================#
# Prepare Output #1
print("\n\n" + "#"*100 + "\n" + "OUTPUT #1\n" + "#"*100 + "\n")

#====================================================================================================#
#%% Add smiles to the dataframe using the dictionary.
data_wi_seqs_df_0 = copy.deepcopy(data_wi_seqs_df)
data_wi_seqs_df_0['smiles'] = data_wi_seqs_df_0["cmpd"].apply(lambda x: cmpd_smls_dict[x.lower()])
data_wi_seqs_df_0 = data_wi_seqs_df_0[data_wi_seqs_df_0['smiles'].str.contains('None') == False]

print("\n\n"+"="*90+"\n#4.1 Dataframe after adding smiles string through looking up compound-smiles dictionary, data_wi_seqs_df_0: ")
beautiful_print(data_wi_seqs_df_0)




###################################################################################################################
#            .g8""8q. `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq.`7MMF'   `7MF'MMP""MM""YMM                             #
#          .dP'    `YM. MM       M  P'   MM   `7   MM   `MM. MM       M  P'   MM   `7     __,         __,         #
#          dM'      `MM MM       M       MM        MM   ,M9  MM       M       MM         `7MM        `7MM         #
#          MM        MM MM       M       MM        MMmmdM9   MM       M       MM           MM          MM         #
#          MM.      ,MP MM       M       MM        MM        MM       M       MM           MM          MM         #
#          `Mb.    ,dP' YM.     ,M       MM        MM        YM.     ,M       MM           MM  ,,      MM         #
#            `"bmmd"'    `bmmmmd"'     .JMML.    .JMML.       `bmmmmd"'     .JMML.       .JMML.db    .JMML.       #
###################################################################################################################
# The following two parts generate  output 1.1 === output 1.2 

#%% Getting unambigious dataset with Uniprot ID only sequences
print("Getting unambigious dataset with Uniprot ID only sequences. First approach.")
# Add SMILES to the dataframe.
data_wi_seqs_df_1 = copy.deepcopy(data_wi_seqs_df_0)
data_wi_seqs_df_1['cmpd_seqs_pair'] = data_wi_seqs_df_1['smiles'] + ',' + data_wi_seqs_df_1['Sequence']
data_wi_seqs_df_1[data_name] = pd.to_numeric(data_wi_seqs_df_1[data_name], downcast="float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD. 
data_avg_value_df = data_wi_seqs_df_1.groupby('cmpd_seqs_pair').agg(
    val_avg = pd.NamedAgg(column = data_name, aggfunc = 'mean') , 
    val_SD  = pd.NamedAgg(column = data_name, aggfunc = 'std' ) ,
                                                                   )


print("\n\n"+"="*90+"\n#4.2 Taking average and standard deviation of the value, data_avg_value_df: ")
beautiful_print(data_avg_value_df)


data_avg_aggre_df = data_wi_seqs_df_1.groupby('cmpd_seqs_pair').agg(lambda column: "~".join(column))
data_wi_seqs_avg_df = pd.concat([data_avg_aggre_df, data_avg_value_df], axis=1)                                   
data_wi_seqs_avg_df['organism'  ] = data_wi_seqs_avg_df['organism'  ].apply(lambda x : str(x).split("~")[0])
data_wi_seqs_avg_df['smiles'    ] = data_wi_seqs_avg_df['smiles'    ].apply(lambda x : str(x).split("~")[0])
data_wi_seqs_avg_df['EC'        ] = data_wi_seqs_avg_df['EC'        ].apply(lambda x : str(x).split("~")[0])
data_wi_seqs_avg_df['cmpd'      ] = data_wi_seqs_avg_df['cmpd'      ].apply(lambda x : str(x).split("~")[0])
data_wi_seqs_avg_df['Sequence'  ] = data_wi_seqs_avg_df['Sequence'  ].apply(lambda x : str(x).split("~")[0])
data_wi_seqs_avg_df['uniprot_id'] = data_wi_seqs_avg_df['uniprot_id'].apply(lambda x : str(x).split("~")[0])

print("\n\n"+"="*90+"\n#4.3 Add avg and std.dev to the aggregated dataframe, data_wi_seqs_avg_df: ")
beautiful_print(data_wi_seqs_avg_df)


# Removing any datapoints where Km SD is more than or equal to 1 
data_wi_seqs_avg_df_small_sd  = data_wi_seqs_avg_df.loc[data_wi_seqs_avg_df['val_SD'] < val_SD_threshold ]
data_wi_seqs_avg_df_one_value = data_wi_seqs_avg_df.loc[data_wi_seqs_avg_df['val_SD'].isnull()           ]
data_wi_seqs_avg_df_screened = pd.concat([data_wi_seqs_avg_df_small_sd, data_wi_seqs_avg_df_one_value], axis=0)
data_wi_seqs_avg_df_screened = data_wi_seqs_avg_df_screened.dropna(subset=['val_avg'])
data_wi_seqs_avg_df_screened.rename(columns={'val_avg': data_name}, inplace=True)
data_wi_seqs_avg_df_screened.reset_index(inplace=True)
data_wi_seqs_avg_df_screened.drop_duplicates(subset=["cmpd_seqs_pair", ], inplace=True)

print("\n\n"+"="*90+"\n#4.4 Removing any datapoints where Km SD is more than or equal to 1, data_wi_seqs_avg_df_screened: ")
beautiful_print(data_wi_seqs_avg_df_screened)


# Final Check, print the number of sequence and smiles pairs.
seqs_list = data_wi_seqs_avg_df_screened["Sequence"].tolist()
smls_list = data_wi_seqs_avg_df_screened["smiles"  ].tolist()
valu_list = data_wi_seqs_avg_df_screened[data_name ].tolist()

pair_list = [(smls, seqs) for smls, seqs in zip(smls_list, seqs_list)]
uniq_pair_list = list(set(pair_list))

print("\n\n_wi_seqs_avg_val_screened, len(uniq_pair_list): ", len(uniq_pair_list))

# Output a processed file in which,
# 1. Only data with clarified uniprot ID's are used.
# 2. Sequences are added to the dataset based on uniprot ID's.
# 3. Compound SMILES string representations are obtained and added.
# 4. Average values are taken if multiple records exist for the same sequence and compound.
# 5. Leave out any records if multiple records for the same sequence and compound exist and variance is big.

data_wi_seqs_avg_df_screened.reset_index(inplace=True)
data_wi_seqs_avg_df_screened = data_wi_seqs_avg_df_screened[["smiles", "Sequence", data_name]]
data_wi_seqs_avg_df_screened.rename(columns={'smiles': "CMPD_SMILES", "Sequence": "SEQ"}, inplace=True)
#data_wi_seqs_avg_df_screened.to_csv(output_folder / ("BRENDA_" + data_name + '_wi_unip_avg_val_screened.csv'))
data_wi_seqs_avg_df_screened.to_csv(output_folder / ("BRENDA_" + data_name + output_1_tail))
del data_wi_seqs_avg_df_screened


###################################################################################################################
#            .g8""8q. `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq.`7MMF'   `7MF'MMP""MM""YMM                             #
#          .dP'    `YM. MM       M  P'   MM   `7   MM   `MM. MM       M  P'   MM   `7     __,                     #
#          dM'      `MM MM       M       MM        MM   ,M9  MM       M       MM         `7MM        pd*"*b.      #
#          MM        MM MM       M       MM        MMmmdM9   MM       M       MM           MM       (O)   j8      #
#          MM.      ,MP MM       M       MM        MM        MM       M       MM           MM           ,;j9      #
#          `Mb.    ,dP' YM.     ,M       MM        MM        YM.     ,M       MM           MM  ,,    ,-='         #
#            `"bmmd"'    `bmmmmd"'     .JMML.    .JMML.       `bmmmmd"'     .JMML.       .JMML.db   Ammmmmmm      #
###################################################################################################################
# A different Approach.
#%% Getting unambigious dataset with Uniprot ID only sequences.
print("Getting unambigious dataset with Uniprot ID only sequences. A different approach.")
data_wi_seqs_df_2 = copy.deepcopy(data_wi_seqs_df_0)

data_wi_seqs_df_2[data_name] = pd.to_numeric(data_wi_seqs_df_2[data_name], downcast = "float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD.
data_avg_value_df = data_wi_seqs_df_2.groupby(['smiles', 'Sequence']).agg(
                                                     val_avg = pd.NamedAgg(column = data_name, aggfunc = 'mean'), 
                                                     val_SD  = pd.NamedAgg(column = data_name, aggfunc = 'std' ), 
                                                    )

# Need to keep those without a val_SD, which means there is only one record for that sequence and smiles.
# Fillna in order to keep those rows when applying the val_SD < val_SD_threshold filter.
data_avg_value_df["val_SD"] = data_avg_value_df["val_SD"].fillna(value=0)                                                         
data_avg_value_df.reset_index(inplace=True)

# Removing any datapoints where Km SD is more than or equal to 1.
data_avg_value_df_screened = data_avg_value_df[data_avg_value_df['val_SD'] < val_SD_threshold]


data_wi_seqs_avg_df_screened = pd.merge(data_wi_seqs_df_2, data_avg_value_df_screened, on=['smiles', 'Sequence'], how = 'left')
data_wi_seqs_avg_df_screened.drop(columns = [data_name, ], inplace = True)


data_wi_seqs_avg_df_screened = data_wi_seqs_avg_df_screened.dropna(subset=['val_avg'])
data_wi_seqs_avg_df_screened.rename(columns={'val_avg': data_name}, inplace=True)
data_wi_seqs_avg_df_screened.drop_duplicates(subset=["smiles", "Sequence"], keep='first', inplace=True)

print("\n\n"+"="*90+"\n#4.? Removing any datapoints where Km SD is more than or equal to 1, data_wi_seqs_avg_df_screened: ")
beautiful_print(data_wi_seqs_avg_df_screened)


# Final Check, print the number of sequence and smiles pairs.
seqs_list = data_wi_seqs_avg_df_screened["Sequence"].tolist()
smls_list = data_wi_seqs_avg_df_screened["smiles"  ].tolist()
valu_list = data_wi_seqs_avg_df_screened[data_name ].values.tolist()
pair_list = [(smls, seqs) for smls, seqs in zip(smls_list, seqs_list)]
uniq_pair_list = list(set(pair_list))
print("\n\n_wi_seqs_avg_val_screened, len(uniq_pair_list): ", len(uniq_pair_list))

data_wi_seqs_avg_df_screened.reset_index(inplace=True)
data_wi_seqs_avg_df_screened = data_wi_seqs_avg_df_screened[["smiles", "Sequence", data_name]]
data_wi_seqs_avg_df_screened.rename(columns={'smiles': "CMPD_SMILES", "Sequence": "SEQ"}, inplace=True)
# No need to output.
#data_wi_seqs_avg_df_screened.to_csv(output_folder / ("BRENDA_" + data_name + '_wi_unip_avg_val_screened_2.csv'))






#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#     M******           .g8""8q.   `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM         ,M' dP   
#    .M               .dP'    `YM.   MM       M  P'   MM   `7   MM   `MM.  MM       M  P'   MM   `7         dP .M'   
#    |bMMAg.          dM'      `MM   MM       M       MM        MM   ,M9   MM       M       MM           mmmMmmMmm    pd*"*b.
#         `Mb         MM        MM   MM       M       MM        MMmmdM9    MM       M       MM             MP dP     (O)   j8
#          jM         MM.      ,MP   MM       M       MM        MM         MM       M       MM          mmdMmmMmmm       ,;j9
#    (O)  ,M9  ,,     `Mb.    ,dP'   YM.     ,M       MM        MM         YM.     ,M       MM           ,M' dP       ,-='   
#     6mmm9    db       `"bmmd"'      `bmmmmd"'     .JMML.    .JMML.        `bmmmmd"'     .JMML.         dP ,M'      Ammmmmmm
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
# Output #2 Dealing with data without uniprot ids.
# We had a list of options here. However, the most reasonable approach should be the last one.
# Output2.1.1 : (1) combine wi_unip and wo_unip -> (2) groupby smiles & sequence -> (3) screen SD < 1.0 [METHOD #1]
# Output2.1.2 : (1) combine wi_unip and wo_unip -> (2) groupby smiles & sequence -> (3) screen SD < 1.0 [METHOD #2] (=== Output2.1.1)
# Output2.1.3 : (1) combine wi_unip and wo_unip -> (2) groupby org, EC & cmpd -> (3) screen SD < 1.0 [METHOD #2]
# Output2.1.4 : (1) groupby org, EC & cmpd for wo_unip_data -> (2) screen SD < 1.0 -> (3) combine wi_unip and wo_unip [METHOD #2] <<< The largest dataset
# Output final: (1) combine wi_unip and wo_unip -> (2) groupby org, EC & cmpd -> (3) screen SD < 1.0 -> (4) Add back all wi_unip_data [METHOD #2]


###################################################################################################################
#         `7MM           mm                                                                       db              #
#           MM           MM                                                                                       #
#      ,M""bMM   ,6"Yb.mmMMmm  ,6"Yb.     `7M'    ,A    `MF' ,pP""Yq.     `7MM  `7MM `7MMpMMMb. `7MM `7MMpdMAo.   #
#    ,AP    MM  8)   MM  MM   8)   MM       VA   ,VAA   ,V  6W'    `Wb      MM    MM   MM    MM   MM   MM   `Wb   #
#    8MI    MM   ,pm9MM  MM    ,pm9MM        VA ,V  VA ,V   8M      M8      MM    MM   MM    MM   MM   MM    M8   #
#    `Mb    MM  8M   MM  MM   8M   MM         VVV    VVV    YA.    ,A9      MM    MM   MM    MM   MM   MM   ,AP   #
#     `Wbmd"MML.`Moo9^Yo.`Mbmo`Moo9^Yo.        W      W      `Ybmmd9'       `Mbod"YML.JMML  JMML.JMML. MMbmmd'    #
#                                                                                                   MM            #
#                                                                                                 .JMML.          #
###################################################################################################################
print("\n\n"+"#"*100 + "\n" + "OUTPUT #2\n" + "#"*100 + "\n")

#====================================================================================================#
#%% Processing EC Number to Sequence Data from KinMod.
Seq_EC_Org_df = pd.read_csv(data_folder / ".." /"Seq_EC_Org_KINMOD.csv")
Seq_EC_Org_df.columns.values[0] = 'EC'
Seq_EC_Org_df.columns.values[1] = 'Sequence'
Seq_EC_Org_df.columns.values[2] = 'organism'

# Modify the organism string.
orgm_list = [" ".join(str(orgm).lower().split(" ")[0:2]) for orgm in Seq_EC_Org_df["organism"]]
Seq_EC_Org_df['organism'] = orgm_list
del orgm_list

# Get a dataframe that works as a organism + EC -> Sequence dictionary.
Seq_EC_Org_grouped_df = Seq_EC_Org_df.groupby(['organism', 'EC'])['Sequence'].apply(list)
Seq_EC_Org_grouped_df = Seq_EC_Org_grouped_df.to_frame().reset_index(inplace = False)


#====================================================================================================#
# Create a dictionary that reads org and ec# and find sequence. (Alternative).
'''
Seq_EC_Org_df = pd.read_csv(data_folder / ".." /"Seq_EC_Org_KINMOD.csv", 
                            index_col = None , 
                            header    = 0    , 
                            sep       = ","  , )

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
'''


#====================================================================================================#
#%% Assign sequence for pairs of organism and ec number (if found in the dict.)
# Modify the organism string.
orgm_list = [" ".join(str(orgm).lower().split(" ")[0:2]) for orgm in data_wo_seqs_df["organism"]]
data_wo_seqs_df['organism'] = orgm_list
del orgm_list

# Merging the Km df with the sequence df to get a combined df.
data_assgn_seqs_df = data_wo_seqs_df.drop(columns = 'Sequence').merge(Seq_EC_Org_grouped_df, on=['organism', 'EC'], how='left')
data_assgn_seqs_df = data_assgn_seqs_df.explode('Sequence') # Exploding the sequence data to ensure 1 substrate is matched to 1 sequence
#data_assgn_seqs_df = data_assgn_seqs_df.loc[data_assgn_seqs_df["Sequence"].isnull() == False]

print("\n\n"+"="*90+"\n#5.0 Dataframe without sequences got sequences assigned, data_assgn_seqs_df: ")
beautiful_print(data_assgn_seqs_df)


#====================================================================================================#
# Splitting the values

data_values_list = [[float(x) for x in one_value.split(',')] for one_value in data_assgn_seqs_df[data_name].tolist()]
data_assgn_seqs_df[data_name] = data_values_list

# for i in range(len(data_assgn_seqs_df)):
#     data_values = data_assgn_seqs_df[data_name].iloc[i]
#     data_assgn_seqs_df[data_name].iloc[i] = [float(x) for x in data_values.split(',')]

# Exploding the Km values to ensure each substrate is matched to 1 Km value and 1 sequence.
data_assgn_seqs_df = data_assgn_seqs_df.explode(data_name) 

print("\n\n"+"="*90+"\n#5.1.1 Dataframe without sequences got values exploded, data_assgn_seqs_df: ")
beautiful_print(data_assgn_seqs_df)


#====================================================================================================#
#%% Prepare outputs for data without uniprot ids.
# Drop data smaller than or equal to zero.
data_assgn_seqs_df.sort_values(by = [data_name], inplace = True)
data_assgn_seqs_df.reset_index(inplace = True)
data_assgn_seqs_df = data_assgn_seqs_df[data_assgn_seqs_df[data_name] > 0].dropna(subset=[data_name])
print("After removing zeros and negative values, len(data_assgn_seqs_df): ", len(data_assgn_seqs_df))

# Check Null values in the dataframe
data_assgn_seqs_df = data_assgn_seqs_df[~data_assgn_seqs_df['Sequence'].isnull()]
print("After removing null in sequence column, len(data_assgn_seqs_df): ", len(data_assgn_seqs_df))

# Get a combined dataframe that contains data with known sequences and data with inferred sequences.
data_all_seqs_combined_df = pd.concat([data_assgn_seqs_df, data_wi_seqs_df], axis=0)
print("After combining with data_wi_seqs_df, len(data_all_seqs_combined_df): ", len(data_all_seqs_combined_df))

# Removing any duplicate rows, === drop_duplicates().
bool_series = data_all_seqs_combined_df.duplicated(keep = 'first')
data_all_seqs_combined_df = data_all_seqs_combined_df[~bool_series]
print("After removing duplicated rows, len(data_all_seqs_combined_df): ", len(data_all_seqs_combined_df))

data_all_seqs_combined_df.drop(columns = ["index", ], inplace = True)
print("\n\n"+"="*90+"\n#5.1.2 Process data both with & without uniprot and sequence, data_all_seqs_combined_df: ")
beautiful_print(data_all_seqs_combined_df)


# Now going to make output files based on:
#--------------------------------------------------------------------------------------------------|
#      DATA WO SMILES       |        DATA WI SMILES       |        DATA AVERAGED & SCREENED        |
#---------------------------|-----------------------------|----------------------------------------|
#      data_wi_seqs_df      |       data_wi_seqs_df_0     |       data_wi_seqs_avg_df_screened     |
#    data_assgn_seqs_df     |    data_assgn_seqs_df_0     |     data_assgn_seqs_avg_df_screened    |
# data_all_seqs_combined_df | data_all_seqs_combined_df_0 | data_all_seqs_combined_avg_df_screened |
#--------------------------------------------------------------------------------------------------|
# Already have data_wi_seqs_df, data_wi_seqs_df_0, data_wi_seqs_avg_df_screened AND data_assgn_seqs_df, data_all_seqs_combined_df.

# Add smiles to the dataframe.
data_assgn_seqs_df_0 = copy.deepcopy(data_assgn_seqs_df)
data_assgn_seqs_df_0['smiles'] = data_assgn_seqs_df_0["cmpd"].apply(lambda x: cmpd_smls_dict[x.lower()])
data_assgn_seqs_df_0 = data_assgn_seqs_df_0[data_assgn_seqs_df_0['smiles'].str.contains('None') == False]

# Add smiles to the dataframe.
data_all_seqs_combined_df_0 = copy.deepcopy(data_all_seqs_combined_df)
data_all_seqs_combined_df_0['smiles'] = data_all_seqs_combined_df_0["cmpd"].apply(lambda x: cmpd_smls_dict[x.lower()])
data_all_seqs_combined_df_0 = data_all_seqs_combined_df_0[data_all_seqs_combined_df_0['smiles'].str.contains('None') == False]


###################################################################################################################
#                                                          __,           __,                                      #
#                                         pd*"*b.         `7MM          `7MM                                      #
#                                        (O)   j8           MM            MM                                      #
#                                            ,;j9           MM            MM                                      #
#                                         ,-='     ,,       MM   ,,       MM                                      #
#                                        Ammmmmmm  db     .JMML. db     .JMML.                                    #
###################################################################################################################
# Getting data_all_seqs_combined_avg_df_screened, output2.1.1 (which is same as output 2.1.2)
# (1) combine wi_unip and wo_unip -> (2) groupby smiles & sequence -> (3) screen SD < 1.0 [METHOD #1]
#%% Taking Average of Km values
# Get pairs.
data_all_seqs_combined_df_1 = copy.deepcopy(data_all_seqs_combined_df_0)
data_all_seqs_combined_df_1['cmpd_seqs_pair'] = data_all_seqs_combined_df_1['smiles'] + ',' + data_all_seqs_combined_df_1['Sequence']
data_all_seqs_combined_df_1[data_name] = pd.to_numeric(data_all_seqs_combined_df_1[data_name], downcast = "float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD.
data_avg_value_df = data_all_seqs_combined_df_1.groupby('cmpd_seqs_pair').agg( 
                            val_avg = pd.NamedAgg(column = data_name, aggfunc = 'mean') , 
                            val_SD  = pd.NamedAgg(column = data_name, aggfunc = 'std' ) ,
                                                                                         )
data_avg_aggre_df = data_all_seqs_combined_df_1.groupby('cmpd_seqs_pair').agg(lambda column: "~".join(column))

print("\n\n"+"="*90+"\n#5.2.1.1 Taking average and standard deviation of the value, data_avg_value_df: ")
beautiful_print(data_avg_value_df)

print("\n\n"+"="*90+"\n#5.2.1.2 Aggregate the dataframe, data_avg_aggre_df: ")
beautiful_print(data_avg_aggre_df)

# Get a dataframe with averaged values.
data_all_seqs_combined_avg_df = pd.concat([data_avg_aggre_df, data_avg_value_df], axis=1)                                     
data_all_seqs_combined_avg_df['organism'  ] = data_all_seqs_combined_avg_df['organism'  ].apply(lambda x : str(x).split("~")[0])
data_all_seqs_combined_avg_df['smiles'    ] = data_all_seqs_combined_avg_df['smiles'    ].apply(lambda x : str(x).split("~")[0])
data_all_seqs_combined_avg_df['EC'        ] = data_all_seqs_combined_avg_df['EC'        ].apply(lambda x : str(x).split("~")[0])
data_all_seqs_combined_avg_df['cmpd'      ] = data_all_seqs_combined_avg_df['cmpd'      ].apply(lambda x : str(x).split("~")[0])
data_all_seqs_combined_avg_df['Sequence'  ] = data_all_seqs_combined_avg_df['Sequence'  ].apply(lambda x : str(x).split("~")[0])
data_all_seqs_combined_avg_df['uniprot_id'] = data_all_seqs_combined_avg_df['uniprot_id'].apply(lambda x : str(x).split("~")[0])

print("\n\n"+"="*90+"\n#5.2.1.3 Add avg and std.dev to the aggregated dataframe, data_all_seqs_combined_avg_df: ")
beautiful_print(data_all_seqs_combined_avg_df)

# Removing any datapoints where Km SD is more than or equal to 1 
data_avg_df_small_sd  = data_all_seqs_combined_avg_df.loc[data_all_seqs_combined_avg_df['val_SD'] < val_SD_threshold ]
data_avg_df_one_value = data_all_seqs_combined_avg_df.loc[data_all_seqs_combined_avg_df['val_SD'].isnull()           ]
data_all_seqs_combined_avg_df_screened = pd.concat([data_avg_df_small_sd, data_avg_df_one_value], axis=0)
data_all_seqs_combined_avg_df_screened = data_all_seqs_combined_avg_df_screened.dropna(subset=['val_avg'])
data_all_seqs_combined_avg_df_screened.rename(columns={'val_avg': data_name}, inplace=True)
data_all_seqs_combined_avg_df_screened.reset_index(inplace=True)
data_all_seqs_combined_avg_df_screened.drop_duplicates(subset=["cmpd_seqs_pair", ], inplace=True)

print("\n\n"+"="*90+"\n#5.2.1.4 Removing any datapoints where Km SD is more than or equal to 1, data_all_seqs_combined_avg_df_screened: ")
beautiful_print(data_all_seqs_combined_avg_df_screened)

# Output.
data_all_seqs_combined_avg_df_screened.reset_index(inplace=True)
data_all_seqs_combined_avg_df_screened = data_all_seqs_combined_avg_df_screened[["smiles", "Sequence", data_name]]
data_all_seqs_combined_avg_df_screened.rename(columns={'smiles': "CMPD_SMILES", "Sequence": "SEQ"}, inplace=True)
#data_all_seqs_combined_avg_df_screened.to_csv(output_folder / ("BRENDA_" + data_name + '_wi_wo_unip_avg_val_screened_2_1_1.csv'))
del data_all_seqs_combined_avg_df_screened

###################################################################################################################
#                                                          __,                                                    #
#                                         pd*"*b.         `7MM           pd*"*b.                                  #
#                                        (O)   j8           MM          (O)   j8                                  #
#                                            ,;j9           MM              ,;j9                                  #
#                                         ,-='     ,,       MM   ,,      ,-='                                     #
#                                        Ammmmmmm  db     .JMML. db     Ammmmmmm                                  #
###################################################################################################################
# Getting data_all_seqs_combined_avg_df_screened, output2.1.2 (which is same as output 2.1.1)
# (1) combine wi_unip and wo_unip -> (2) groupby smiles & sequence -> (3) screen SD < 1.0 [METHOD #2]
#%% Taking Average of Km values
data_all_seqs_combined_df_2 = copy.deepcopy(data_all_seqs_combined_df_0)
data_all_seqs_combined_df_2[data_name] = pd.to_numeric(data_all_seqs_combined_df_2[data_name], downcast = "float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD.
data_avg_value_df = data_all_seqs_combined_df_2.groupby(['smiles', 'Sequence']).agg(
                                                     val_avg = pd.NamedAgg(column = data_name, aggfunc = 'mean'), 
                                                     val_SD  = pd.NamedAgg(column = data_name, aggfunc = 'std' ), 
                                                    )

# Need to keep those without a val_SD, which means there is only one record for that sequence and smiles.
# Fillna in order to keep those rows when applying the val_SD < val_SD_threshold filter.
data_avg_value_df["val_SD"] = data_avg_value_df["val_SD"].fillna(value=0)                                                         
data_avg_value_df.reset_index(inplace=True)

# Removing any datapoints where Km SD is more than or equal to 1.
data_avg_value_df_screened = data_avg_value_df[data_avg_value_df['val_SD'] < val_SD_threshold]


data_all_seqs_combined_avg_df_screened = pd.merge(data_all_seqs_combined_df_2, data_avg_value_df_screened, on=['smiles', 'Sequence'], how = 'left')
data_all_seqs_combined_avg_df_screened.drop(columns = [data_name, ], inplace = True)


data_all_seqs_combined_avg_df_screened = data_all_seqs_combined_avg_df_screened.dropna(subset=['val_avg'])
data_all_seqs_combined_avg_df_screened.rename(columns={'val_avg': data_name}, inplace=True)
data_all_seqs_combined_avg_df_screened.drop_duplicates(subset=["smiles", "Sequence"], keep='first', inplace=True)

print("\n\n"+"="*90+"\n#5.3.? Removing any datapoints where Km SD is more than or equal to 1, data_all_seqs_combined_avg_df_screened: ")
beautiful_print(data_all_seqs_combined_avg_df_screened)

data_all_seqs_combined_avg_df_screened.reset_index(inplace=True)
data_all_seqs_combined_avg_df_screened = data_all_seqs_combined_avg_df_screened[["smiles", "Sequence", data_name]]
data_all_seqs_combined_avg_df_screened.rename(columns={'smiles': "CMPD_SMILES", "Sequence": "SEQ"}, inplace=True)
#data_all_seqs_combined_avg_df_screened.to_csv(output_folder / ("BRENDA_" + data_name + '_wi_wo_unip_avg_val_screened_2_1_2.csv'))

###################################################################################################################
#                                                       __,                                                       #
#                                       pd*"*b.        `7MM              ,AM                                      #
#                                      (O)   j8          MM             AVMM                                      #
#                                          ,;j9          MM           ,W' MM                                      #
#                                       ,-='    ,,       MM  ,,     ,W'   MM                                      #
#                                      Ammmmmmm db     .JMML.db     AmmmmmMMmm                                    #
#                                                                         MM                                      #
#                                                                         MM                                      #
###################################################################################################################
# Output the largest dataset.
# (1) combine wi_unip and wo_unip -> (2) groupby smiles & sequence -> (3) screen SD < 1.0 [METHOD #2]
#%% Taking Average of Km values
data_assgn_seqs_df_1 = copy.deepcopy(data_assgn_seqs_df_0)
data_assgn_seqs_df_1[data_name] = pd.to_numeric(data_assgn_seqs_df_1[data_name], downcast = "float")

# Processing compound, sequence pairs that have duplicates. Taking the average and the SD.
data_avg_value_df = data_assgn_seqs_df_1.groupby(['smiles', 'Sequence']).agg(
                                                     val_avg = pd.NamedAgg(column = data_name, aggfunc = 'mean'), 
                                                     val_SD  = pd.NamedAgg(column = data_name, aggfunc = 'std' ), 
                                                    )

# Need to keep those without a val_SD, which means there is only one record for that sequence and smiles.
# Fillna in order to keep those rows when applying the val_SD < val_SD_threshold filter.
data_avg_value_df["val_SD"] = data_avg_value_df["val_SD"].fillna(value=0)                                                         
data_avg_value_df.reset_index(inplace=True)

# Removing any datapoints where Km SD is more than or equal to 1.
data_avg_value_df_screened = data_avg_value_df[data_avg_value_df['val_SD'] < val_SD_threshold]


data_assgn_seqs_df_screened = pd.merge(data_assgn_seqs_df_1, data_avg_value_df_screened, on=['smiles', 'Sequence'], how = 'left')
data_assgn_seqs_df_screened.drop(columns = [data_name, ], inplace = True)


data_assgn_seqs_df_screened = data_assgn_seqs_df_screened.dropna(subset=['val_avg'])
data_assgn_seqs_df_screened.rename(columns={'val_avg': data_name}, inplace=True)
data_assgn_seqs_df_screened.drop_duplicates(subset=["smiles", "Sequence"], keep='first', inplace=True)

# Change column names for output.
data_assgn_seqs_df_screened = data_assgn_seqs_df_screened[["smiles", "Sequence", data_name]]
data_assgn_seqs_df_screened.rename(columns={'smiles': "CMPD_SMILES", "Sequence": "SEQ"}, inplace=True)

# Add back all wi_unip data.
data_all_seqs_df_screened = pd.concat([data_assgn_seqs_df_screened, data_wi_seqs_avg_df_screened], axis=0)
data_all_seqs_df_screened.drop_duplicates(subset=['CMPD_SMILES', 'SEQ', ], keep='first', inplace=True)

# Get rid of `index` columns.
data_all_seqs_df_screened.reset_index(inplace=True)
data_all_seqs_df_screened = data_all_seqs_df_screened[["CMPD_SMILES", "SEQ", data_name]]


# Output.
#data_all_seqs_df_screened.to_csv(output_folder / ("BRENDA_" + data_name + '_wi_wo_unip_avg_val_screened_0.csv'))
data_all_seqs_df_screened.to_csv(output_folder / ("BRENDA_" + data_name + output_2_tail))



###################################################################################################################
#                                        ,...  ,,                         ,,                                      #
#                                      .d' ""  db                       `7MM                                      #
#                                      dM`                                MM                                      #
#                                     mMMmm  `7MM  `7MMpMMMb.   ,6"Yb.    MM                                      #
#                                      MM      MM    MM    MM  8)   MM    MM                                      #
#                                      MM      MM    MM    MM   ,pm9MM    MM                                      #
#                                      MM      MM    MM    MM  8M   MM    MM                                      #
#                                    .JMML.  .JMML..JMML  JMML.`Moo9^Yo..JMML.                                    #
###################################################################################################################

data_all_seqs_combined_df_X = copy.deepcopy(data_all_seqs_combined_df_0)
data_all_seqs_combined_df_X[data_name] = pd.to_numeric(data_all_seqs_combined_df_X[data_name], downcast = "float")

# Processing organism, EC, cmpd combinations that have duplicates. Taking the average and the SD.
data_avg_value_df = data_all_seqs_combined_df_X.groupby(['organism', 'EC', 'cmpd']).agg(
                        val_avg     = pd.NamedAgg(column = data_name,    aggfunc = 'mean' ), 
                        val_SD      = pd.NamedAgg(column = data_name,    aggfunc = 'std'  ), 
                        uniprot_id  = pd.NamedAgg(column = 'uniprot_id', aggfunc = 'count'), 
                                                                                )
# Need to keep those without a val_SD, which means there is only one record for that sequence and smiles.
# Fillna in order to keep those rows when applying the val_SD < val_SD_threshold filter.
data_avg_value_df["val_SD"] = data_avg_value_df["val_SD"].fillna(value=0)                                                         
data_avg_value_df.reset_index(inplace=True)

# Creating a dictionary df that has Organism, EC number, and substrate triplicate entries that don't have large variation between Km values and different sequences
data_avg_value_df_screened = data_avg_value_df[data_avg_value_df['val_SD'] < val_SD_threshold]

# Merging datasets and filtering the pangenomic dataset for the Km SD condition above
data_all_seqs_combined_avg_df_screened_X = pd.merge(data_all_seqs_combined_df_X, data_avg_value_df_screened, on=['organism', 'EC', 'cmpd'], how = 'left')    
data_all_seqs_combined_avg_df_screened_X.drop(columns = [data_name, ], inplace = True)

data_all_seqs_combined_avg_df_screened_X = data_all_seqs_combined_avg_df_screened_X.dropna(subset=['val_avg'])
data_all_seqs_combined_avg_df_screened_X.rename(columns={'val_avg': data_name}, inplace=True)
data_all_seqs_combined_avg_df_screened_X.drop_duplicates(subset=['organism', 'EC', 'cmpd'], keep='first', inplace=True)

print("\n\n"+"="*90+"\n#5.X.X Removing any datapoints where Km SD is more than or equal to 1, data_all_seqs_combined_avg_df_screened_X: ")
beautiful_print(data_all_seqs_combined_avg_df_screened_X)

# Change column names for output.
data_all_seqs_combined_avg_df_screened_X = data_all_seqs_combined_avg_df_screened_X[["smiles", "Sequence", data_name]]
data_all_seqs_combined_avg_df_screened_X.rename(columns={'smiles': "CMPD_SMILES", "Sequence": "SEQ"}, inplace=True)

# Add back all wi_unip data.
data_all_seqs_combined_avg_df_screened_X = pd.concat([data_all_seqs_combined_avg_df_screened_X, data_wi_seqs_avg_df_screened], axis=0)
data_all_seqs_combined_avg_df_screened_X.drop_duplicates(subset=['CMPD_SMILES', 'SEQ', ], keep='first', inplace=True)

# Get rid of `index` columns.
data_all_seqs_combined_avg_df_screened_X.reset_index(inplace=True)
data_all_seqs_combined_avg_df_screened_X = data_all_seqs_combined_avg_df_screened_X[["CMPD_SMILES", "SEQ", data_name]]

# Output.
#data_all_seqs_combined_avg_df_screened_X.to_csv(output_folder / ("BRENDA_" + data_name + '_wi_wo_unip_avg_val_screened_X.csv'))
data_all_seqs_combined_avg_df_screened_X.to_csv(output_folder / ("BRENDA_" + data_name + output_3_tail))

print("Done!")

###################################################################################################################
###################################################################################################################
# Additional Check or Analysis.
'''
cmpd_list = data_all_seqs_avg_df_screened["cmpd"    ].tolist()
seqs_list = data_all_seqs_avg_df_screened["Sequence"].tolist()
smls_list = data_all_seqs_avg_df_screened["smiles"  ].tolist()
valu_list = data_all_seqs_avg_df_screened["Km"      ].tolist()

seqs_smls_valu_list = []
for i, (seqs, smls, valu, cmpd) in enumerate(zip(seqs_list, smls_list, valu_list, cmpd_list)):
    seqs_smls_valu_list.append((seqs, seqs))

del seqs_list, smls_list, valu_list



data_all_final_df_rana = pd.read_csv(data_folder / "brenda_km_pangenomic_chkpt.csv", on_bad_lines = 'skip', header = 0, sep=',', encoding='cp1252')
cmpd_list = data_all_final_df_rana["Substrate"].tolist()
seqs_list = data_all_final_df_rana["Sequence" ].tolist()
smls_list = data_all_final_df_rana["SMILES"   ].tolist()
valu_list = data_all_final_df_rana["Km [mM]"  ].tolist()

seqs_smls_valu_list_rana = []
for i, (seqs, smls, valu, cmpd) in enumerate(zip(seqs_list, smls_list, valu_list, cmpd_list)):
    seqs_smls_valu_list_rana.append((seqs, seqs))

del seqs_list, smls_list, valu_list

questionable_set = set(seqs_smls_valu_list_rana) - set(seqs_smls_valu_list)

print(len(set(seqs_smls_valu_list_rana)))
print(len(set(seqs_smls_valu_list)))

print(len(questionable_set))
for i in questionable_set:
    print(i)
'''


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

