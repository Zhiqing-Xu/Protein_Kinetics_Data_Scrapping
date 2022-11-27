#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
###################################################################################################################
###################################################################################################################
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
import ast
import pickle
import scipy.io
import subprocess

from tqdm import tqdm
from random import shuffle

#--------------------------------------------------#
import numpy as np


#--------------------------------------------------#

from AP_convert import *
from AP_convert import Get_Unique_SMILES

GetUnqSmi = Get_Unique_SMILES(isomericSmiles = True, SMARTS_bool = False)




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# These codes read through the tsv files from MetaNetX database, 
# trying to obtain canonical SMILES string that keeps isomerism,
# and keeps those SMILES strings that is potentially invalid. 


def smiles_hard_coding_process(one_smiles, skip = False):
    if skip == True:
        return one_smiles

    smiles_processed=max(one_smiles.split("."), key = len)
    #smiles_processed=smiles_processed.replace("([*])", "")
    #smiles_processed=smiles_processed.replace("[*]", "")
    smiles_processed=smiles_processed.replace("[O-]", "O")
    smiles_processed=smiles_processed.replace("[N+]", "N")
    smiles_processed=smiles_processed.replace("[NH+]", "N")
    smiles_processed=smiles_processed.replace("[NH2+]", "N")
    smiles_processed=smiles_processed.replace("[NH3+]", "N")
    smiles_processed=smiles_processed.replace("[N-]", "N")
    #smiles_processed=smiles_processed.replace("*-", "N")
    #smiles_processed=smiles_processed.replace("-*", "N")
    #smiles_processed=smiles_processed.replace("(-*)", "N")
    return smiles_processed

def get_longest_smiles_from_aggregates(one_smiles_string):
    return max(one_smiles_string.split("."), key = len)

#====================================================================================================#
def readlines_chem_prop(file, mark="", skip_header=False):
    count_x=0
    smiles_id_nme_dict=dict([])
    # skip lines not containing information
    # There are 388 lines of header (trivial information about the DB)
    if skip_header==True:
        for i in range(388): 
            line = file.readline()

    # Process 10^6 at a time in case there are problematic lines in the data file.
    # It is easy to combine the dictionaries with a few lines of code.
    for i in tqdm(range(100000)): 
        count_x+=1
        line = file.readline()
        if line != "" and line != "\n" :
            cmpd_info_list = line.split("\t")
            #print (cmpd_info_list[0], count_x+1)

            if cmpd_info_list[1].find("compound")==-1:
                one_smiles = cmpd_info_list[6]

                smiles_processed = get_longest_smiles_from_aggregates(one_smiles)
                smiles_processed = GetUnqSmi.UNQSMI(smiles_processed)

                smiles_id_nme_dict[smiles_processed] = [cmpd_info_list[1], cmpd_info_list[0],]
        else:
            return smiles_id_nme_dict
    return smiles_id_nme_dict

#====================================================================================================#
def readlines_chem_prop2(file, mark="", skip_header=False):
    count_x=0
    smiles_id_nme_dict=dict([])
    # skip lines not containing information
    # There are 388 lines of header (trivial information about the DB)
    if skip_header==True:
        for i in range(388): 
            line = file.readline()

    # Process 10^6 at a time in case there are problematic lines in the data file.
    # It is easy to combine the dictionaries with a few lines of code.
    for i in tqdm(range(100000)): 
        count_x+=1
        line = file.readline()
        if line != "" and line != "\n" :
            cmpd_info_list=line.split("\t")
            #print (cmpd_info_list[0], count_x+1)

            if cmpd_info_list[1].find("compound")==-1:
                one_smiles=cmpd_info_list[6]

                smiles_processed = get_longest_smiles_from_aggregates(one_smiles)
                smiles_processed = GetUnqSmi.UNQSMI(smiles_processed)

                smiles_id_nme_dict[cmpd_info_list[1].lower()] = [smiles_processed, cmpd_info_list[0],]
        else:
            return smiles_id_nme_dict
    return smiles_id_nme_dict


###################################################################################################################
###################################################################################################################
def main():
    #====================================================================================================#
    # These codes read through the tsv files from MetaNetX database, 
    # trying to obtain canonical SMILES string that keeps isomerism,
    # and keeps those SMILES strings that is potentially invalid. 
    #====================================================================================================#
    # Process 10^6 at a time in case there are problematic lines in the data file.
    # Write the compound info into 7 dictionaries, each info of with up to 10^6 compounds

    


    '''
    file_MNX_cmpd_address="./MNXScrapSavings/raw_data/chem_prop.tsv"
    file_MNX_cmpd=open(file_MNX_cmpd_address)
    smiles_id_nme_dict1=readlines_chem_prop(file_MNX_cmpd,"",True)
    smiles_id_nme_dict2=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict3=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict4=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict5=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict6=readlines_chem_prop(file_MNX_cmpd,"",False)
    smiles_id_nme_dict7=readlines_chem_prop(file_MNX_cmpd,"",False)

    pickle_out1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict1","wb")
    pickle.dump(smiles_id_nme_dict1, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict2","wb")
    pickle.dump(smiles_id_nme_dict2, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict3","wb")
    pickle.dump(smiles_id_nme_dict3, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict4","wb")
    pickle.dump(smiles_id_nme_dict4, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict5","wb")
    pickle.dump(smiles_id_nme_dict5, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict6","wb")
    pickle.dump(smiles_id_nme_dict6, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict7","wb")
    pickle.dump(smiles_id_nme_dict7, pickle_out1)
    pickle_out1.close()
    
    #====================================================================================================#
    
    file_MNX_cmpd_address="./MNXScrapSavings/raw_data/chem_prop.tsv"
    file_MNX_cmpd=open(file_MNX_cmpd_address)
    nme_smiles_id_dict1=readlines_chem_prop2(file_MNX_cmpd,"",True)
    nme_smiles_id_dict2=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict3=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict4=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict5=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict6=readlines_chem_prop2(file_MNX_cmpd,"",False)
    nme_smiles_id_dict7=readlines_chem_prop2(file_MNX_cmpd,"",False)

    pickle_out1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict1","wb")
    pickle.dump(nme_smiles_id_dict1, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict2","wb")
    pickle.dump(nme_smiles_id_dict2, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict3","wb")
    pickle.dump(nme_smiles_id_dict3, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict4","wb")
    pickle.dump(nme_smiles_id_dict4, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict5","wb")
    pickle.dump(nme_smiles_id_dict5, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict6","wb")
    pickle.dump(nme_smiles_id_dict6, pickle_out1)
    pickle_out1.close()
    pickle_out1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict7","wb")
    pickle.dump(nme_smiles_id_dict7, pickle_out1)
    pickle_out1.close()    
    '''





    #====================================================================================================#
    #====================================================================================================#
    # Open all pickles and retrieve compound/reaction info.
    pickle_in1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict1","rb")
    smiles_id_nme_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict2","rb")
    smiles_id_nme_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict3","rb")
    smiles_id_nme_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict4","rb")
    smiles_id_nme_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict5","rb")
    smiles_id_nme_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict6","rb")
    smiles_id_nme_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_smiles_id_nme_dict7","rb")
    smiles_id_nme_dict7=pickle.load(pickle_in1)
    pickle_in1.close()

    smiles_id_nme_dict=smiles_id_nme_dict1.copy()
    smiles_id_nme_dict.update(smiles_id_nme_dict2)
    smiles_id_nme_dict.update(smiles_id_nme_dict3)
    smiles_id_nme_dict.update(smiles_id_nme_dict4)
    smiles_id_nme_dict.update(smiles_id_nme_dict5)
    smiles_id_nme_dict.update(smiles_id_nme_dict6)
    smiles_id_nme_dict.update(smiles_id_nme_dict7)

    # Add meaningless keys to the dict to avoid errors
    smiles_id_nme_dict["BIOMASS"] = ["",""]
    smiles_id_nme_dict["MNXM0"]   = ["",""]
    #====================================================================================================#
    # Open all pickles and retrieve compound/reaction info.
    pickle_in1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict1","rb")
    nme_smiles_id_dict1=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict2","rb")
    nme_smiles_id_dict2=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict3","rb")
    nme_smiles_id_dict3=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict4","rb")
    nme_smiles_id_dict4=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict5","rb")
    nme_smiles_id_dict5=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict6","rb")
    nme_smiles_id_dict6=pickle.load(pickle_in1)
    pickle_in1.close()
    pickle_in1=open("./MNXScrapSavings/temp/temp_nme_smiles_id_dict7","rb")
    nme_smiles_id_dict7=pickle.load(pickle_in1)
    pickle_in1.close()

    nme_smiles_id_dict=nme_smiles_id_dict1.copy()
    nme_smiles_id_dict.update(nme_smiles_id_dict2)
    nme_smiles_id_dict.update(nme_smiles_id_dict3)
    nme_smiles_id_dict.update(nme_smiles_id_dict4)
    nme_smiles_id_dict.update(nme_smiles_id_dict5)
    nme_smiles_id_dict.update(nme_smiles_id_dict6)
    nme_smiles_id_dict.update(nme_smiles_id_dict7)

    #====================================================================================================#

    pickle_out1=open("./MNXScrapSavings/CPD01_smiles_MNXid_nme_dict.p","wb")
    pickle.dump(smiles_id_nme_dict, pickle_out1)
    pickle_out1.close()

    pickle_out1=open("./MNXScrapSavings/CPD01_nme_smiles_MNXid_dict.p","wb")
    pickle.dump(nme_smiles_id_dict, pickle_out1)
    pickle_out1.close()

    cofactor_mnx_id=["MNXM0","MNXM1","MNXM2","MNXM3","MNXM4","MNXM5","MNXM6","MNXM7","MNXM8",\
                     "MNXM9","MNXM10","MNXM11","MNXM12","MNXM13","MNXM14","MNXM15","MNXM15"]

    #====================================================================================================#
    print ("Done Merge!")
    print (smiles_id_nme_dict[ GetUnqSmi.UNQSMI( "CC=O" ) ])
    print (nme_smiles_id_dict["acetyl-coa"])
    print ()
    print ()
    print (nme_smiles_id_dict["Malonyl-CoA".lower()])
    print (nme_smiles_id_dict["Succinyl-CoA".lower()])
    print (smiles_id_nme_dict[ GetUnqSmi.UNQSMI( "CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1OP(=O)([O-])[O-]" ) ])

###################################################################################################################
###################################################################################################################
if __name__ == '__main__':
    main()



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



