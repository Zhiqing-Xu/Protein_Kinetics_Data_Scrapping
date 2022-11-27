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
import copy
import pickle
import argparse
import numpy as np
import pandas as pd
#--------------------------------------------------#
from urllib.parse import quote
from urllib.request import urlopen
from unidecode import unidecode
###################################################################################################################
###################################################################################################################
# Basic Functions
# Print the DataFrame obtained.
def beautiful_print(df):
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'      , 20, 
                           'display.min_rows'      , 20, 
                           'display.max_columns'   , 10, 
                           #"display.max_colwidth" , None,
                           "display.width"         , None,
                           "expand_frame_repr"     , True,
                           "max_seq_items"         , None,):  # more options can be specified
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

print(unidecode("Ã¤ndern"))


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
data_folder_name    = ["Ki_BRENDA"        , "KM_3_BRENDA"      , "kcat_BRENDA"        , "kcat_KM_BRENDA"        ][0]
data_file_name      = ["brenda_Ki_raw.csv", "brenda_KM_raw.csv", "brenda_kcat_raw.csv", "brenda_kcat_KM_raw.csv"][0]
output_file_name    = ["Ki_BRENDA.csv"    , "KM_BRENDA.csv"    , "kcat_BRENDA.csv"    , "kcat_KM_BRENDA.csv"    ][0]
data_name           = ["Ki"               , "Km"               , "kcat"               , "kcat_KM"               ][0]
data_name_0         = ["max_Ki"           , "max_Km"           , "max_kcat"           , "kcat_KM_0"             ][0]

data_folder      = Path("X_DataProcessing/X00_enzyme_datasets/" + data_folder_name)
data_file        = data_file_name

output_folder    = Path("X_DataProcessing/X00_enzyme_datasets_processed/")
output_file      = output_file_name


#====================================================================================================#
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


cmpd_list = raw_df_0["cmpd"].tolist()

###################################################################################################################
###################################################################################################################
# Get dict for reading KEGG compound IDs.

KEGG_id_SMILES_by_ChEBI_CAS_df = pd.read_csv("./CMPD_RXN_tools/KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_ChEBI_CAS.csv", header = 0)

KEGGid_list     = KEGG_id_SMILES_by_ChEBI_CAS_df["KEGGid"].tolist()
SMILES_by_ChEBI = KEGG_id_SMILES_by_ChEBI_CAS_df["SMILES_by_ChEBI"].tolist()
SMILES_by_CAS   = KEGG_id_SMILES_by_ChEBI_CAS_df["SMILES_by_CAS"].tolist()

KEGG_id_SMILES_by_ChEBI_CAS_dict = dict([])
for i, id in enumerate(KEGGid_list):
    KEGG_id_SMILES_by_ChEBI_CAS_dict[id] = [SMILES_by_ChEBI[i], SMILES_by_CAS[i]]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Get dict for reading KEGG compound IDs.

KEGG_id_SMILES_by_mol_df = pd.read_csv("./CMPD_RXN_tools/KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_mol_file_no_uniq.csv", header = 0)

KEGGid_list     = KEGG_id_SMILES_by_mol_df["KEGGid"].tolist()
SMILES_by_mol   = KEGG_id_SMILES_by_mol_df["SMILES_by_mol"].tolist()


KEGG_id_SMILES_by_mol_dict = dict([])
for i, id in enumerate(KEGGid_list):
    KEGG_id_SMILES_by_mol_dict[id] = SMILES_by_mol[i]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Get dict for reading ALL_KEGG_DATA_DICT.

folder = Path("./CMPD_RXN_tools/KEGGScrapSavings")
filename = "CPD02_ALL_KEGG_DATA_DICT_bk.p"

pickle_in1 = open( folder / filename, "rb")
ALL_KEGG_DATA_DICT = pickle.load(pickle_in1)
pickle_in1.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Get dict for reading ALL_KEGG_DATA_DICT.

folder = Path("./CMPD_RXN_tools/MNXScrapSavings")
filename = "CPD01_nme_smiles_MNXid_dict.p"

pickle_in1 = open( folder / filename, "rb")
nme_smiles_MNXid_dict = pickle.load(pickle_in1)
pickle_in1.close()


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                  .g8"""bgd       db        .g8"""bgd MMP""MM""YMM `7MMF'   `7MF' .M"""bgd 
#    __,         .dP'     `M      ;MM:     .dP'     `M P'   MM   `7   MM       M  ,MI    "Y 
#   `7MM         dM'       `     ,V^MM.    dM'       `      MM        MM       M  `MMb.     
#     MM         MM             ,M  `MM    MM               MM        MM       M    `YMMNq. 
#     MM         MM.            AbmmmqMA   MM.              MM        MM       M  .     `MM 
#     MM  ,,     `Mb.     ,'   A'     VML  `Mb.     ,'      MM        YM.     ,M  Mb     dM 
#   .JMML.db       `"bmmmd'  .AMA.   .AMMA.  `"bmmmd'     .JMML.       `bmmmmd"'  P"Ybmmd"  
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Use Cactus
def cactus_name_to_smls(cmpd_name):
    try: 
        url = 'http://cactus.nci.nih.gov/chemical/structure/' + quote(cmpd_name) + '/smiles' 
        cmpd_smiles = urlopen(url).read().decode('utf8') 
        return cmpd_smiles 
    except: 
        cmpd_smiles = "None"
        return cmpd_smiles

#====================================================================================================#
# rctt_list & cmpd_list.
rctt_list = copy.deepcopy(cmpd_list)
uniq_rctt_list = list(set(rctt_list))
num_name = len(uniq_rctt_list)
print("Number of unique compounds, len(list(set(rctt_list))): ", num_name)

cmpd_smls_dict = dict([])
with open(data_folder / 'cmpd_smls_cactus.txt', 'w') as f:
    count_x = 0
    f.write(";CMPD;CMPD_SMILES"+"\n")
    for i, one_cmpd_name in enumerate(uniq_rctt_list[0:]):
        if i > -1:
            one_cmpd_smls = cactus_name_to_smls(one_cmpd_name)
            cmpd_smls_dict[one_cmpd_name] = one_cmpd_smls
            if i % 1 == 0:
                print(i, one_cmpd_name , one_cmpd_smls)
            one_cmpd_name = unidecode(one_cmpd_name)
            # write a row to the csv file
            f.write(str(count_x)+";"+str(one_cmpd_name)+";"+str(one_cmpd_smls)+"\n")
            count_x+=1


#====================================================================================================#
# Use cactus again for modified names
cmpd_smls_df = pd.read_csv(data_folder / "cmpd_smls_cactus.txt",
                           encoding  = "ISO-8859-1", 
                           index_col = None        , 
                           header    = 0           , 
                           sep       = ';'         , 
                           engine    = 'python'    , 
                           )

cmpd_list = cmpd_smls_df["CMPD"].tolist()
smls_list = cmpd_smls_df["CMPD_SMILES"].tolist()


with open(data_folder / "cmpd_smls_cactus.csv", 'w') as f:
    f.write(";CMPD;CMPD_SMILES"+"\n")
    for i, (cmpd, smls) in enumerate(zip(cmpd_list, smls_list)):
        if type(cmpd) != str:
            cmpd = str(cmpd)

        if cmpd.find("/in") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("/in", "")
            one_cmpd_smls = cactus_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)

        elif cmpd.find("/out") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("/out", "")
            one_cmpd_smls = cactus_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)

        elif cmpd.find("[side 1]") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("[side 1]", "")
            one_cmpd_smls = cactus_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)

        elif cmpd.find("[side 2]") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("[side 2]", "")
            one_cmpd_smls = cactus_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)

        else:
            f.write(str(i)+";"+cmpd+";"+str(smls)+"\n")
            #print(i, cmpd, smls) 
            


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                  `7MM"""Mq. `7MMF'   `7MF'`7MM"""Yp,   .g8"""bgd `7MMF'  `7MMF'`7MM"""YMM  `7MMM.     ,MMF'
#                    MM   `MM.  MM       M    MM    Yb .dP'     `M   MM      MM    MM    `7    MMMb    dPMM  
#   pd*"*b.          MM   ,M9   MM       M    MM    dP dM'       `   MM      MM    MM   d      M YM   ,M MM  
#  (O)   j8          MMmmdM9    MM       M    MM"""bg. MM            MMmmmmmmMM    MMmmMM      M  Mb  M' MM  
#      ,;j9          MM         MM       M    MM    `Y MM.           MM      MM    MM   Y  ,   M  YM.P'  MM  
#   ,-='    ,,       MM         YM.     ,M    MM    ,9 `Mb.     ,'   MM      MM    MM     ,M   M  `YM'   MM  
#  Ammmmmmm db     .JMML.        `bmmmmd"'  .JMMmmmd9    `"bmmmd'  .JMML.  .JMML..JMMmmmmMMM .JML. `'  .JMML.
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Use Pubchempy
from pubchempy import get_compounds, Compound

def pubchem_name_to_smls(cmpd_name):
    try: 
        cmpd_list = get_compounds(cmpd_name, 'name')
        cmpd_smiles = cmpd_list[0].isomeric_smiles
        return cmpd_smiles 
    except: 
        cmpd_smiles = "None"
        return cmpd_smiles

uniq_rctt_list = []
for one_rctt in rctt_list:
    if one_rctt not in uniq_rctt_list:
        uniq_rctt_list.append(one_rctt)

num_name = len(uniq_rctt_list)
print("len(list(set(rctt_list))): ", num_name)

cmpd_smls_dict = dict([])

with open(data_folder / "cmpd_smls_pubchem.txt", 'w') as f:
    count_x=0
    f.write(";CMPD;CMPD_SMILES"+"\n")
    for i, one_cmpd_name in enumerate(uniq_rctt_list[0:]):
        one_cmpd_smls = pubchem_name_to_smls(one_cmpd_name)
        cmpd_smls_dict[one_cmpd_name] = one_cmpd_smls
        if i % 1 == 0:
            print(i, one_cmpd_name , one_cmpd_smls)
        one_cmpd_name = unidecode(one_cmpd_name)
        # write a row to the csv file
        f.write(str(count_x)+";"+str(one_cmpd_name)+";"+str(one_cmpd_smls)+"\n")
        count_x+=1
        


#====================================================================================================#
# Use pubchem again for modified names

cmpd_smls_df = pd.read_csv(data_folder / "cmpd_smls_pubchem.txt",
                       encoding  = "ISO-8859-1",
                       index_col = None, 
                       header    = 0, 
                       sep       = ';',
                       engine    = 'python')


cmpd_list = cmpd_smls_df["CMPD"].tolist()
smls_list = cmpd_smls_df["CMPD_SMILES"].tolist()


with open(data_folder / "cmpd_smls_pubchem.csv", 'w') as f:
    f.write(";CMPD;CMPD_SMILES"+"\n")
    for i, (cmpd, smls) in enumerate(zip(cmpd_list, smls_list)):
        print(cmpd)
        if type(cmpd) != str:
            cmpd = str(cmpd)

        if cmpd.find("/in") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("/in", "")
            one_cmpd_smls = pubchem_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)

        elif cmpd.find("/out") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("/out", "")
            one_cmpd_smls = pubchem_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)
            
        elif cmpd.find("[side 1]") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("[side 1]", "")
            one_cmpd_smls = pubchem_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)

        elif cmpd.find("[side 2]") != -1 and smls == "None":
            cmpd_modified = cmpd.replace("[side 2]", "")
            one_cmpd_smls = pubchem_name_to_smls(cmpd_modified)
            f.write(str(i)+";"+cmpd+";"+str(one_cmpd_smls)+"\n")
            print(i, cmpd, cmpd_modified, one_cmpd_smls)

        else:
            f.write(str(i)+";"+cmpd+";"+str(smls)+"\n")
            print(i, cmpd, smls)
            

###################################################################################################################
###################################################################################################################

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




