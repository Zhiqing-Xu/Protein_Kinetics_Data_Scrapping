#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Microsoft VS header
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
print("="*80)
#--------------------------------------------------#
if __name__ == "__main__":
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")
#--------------------------------------------------#
    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")
#--------------------------------------------------#
    if "__file__" in globals().keys():
        try:
            os.chdir(os.path.dirname(__file__))
            print('CurrentDir: ', os.getcwd())
        except:
            print("Problems with navigating to the file dir.")
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
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
#--------------------------------------------------#
import pickle
import argparse
import numpy as np
import pandas as pd


#--------------------------------------------------#
import csv
import pickle
import urllib
import requests

#--------------------------------------------------#
# bioservices.KEGG
import cirpy
import pubchempy
from bioservices import * 
from chemspipy import *
cs = ChemSpider("85cfd898-cc63-4347-9dec-0b4964a387c6")
#--------------------------------------------------#
from AP_convert import *

GetUnqSmi = Get_Unique_SMILES(isomericSmiles = True, SMARTS_bool = False)
unis = GetUnqSmi.UNQSMI



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#  `7MM"""Yb. `7MM"""Yp,      .g8"""bgd   .g8""8q. `7MN.   `7MF'`7MMF'   `7MF'`7MM"""YMM  `7MM"""Mq.   .M"""bgd `7MMF' .g8""8q. `7MN.   `7MF'.M"""bgd  #
#    MM    `Yb. MM    Yb    .dP'     `M .dP'    `YM. MMN.    M    `MA     ,V    MM    `7    MM   `MM. ,MI    "Y   MM .dP'    `YM. MMN.    M ,MI    "Y  #
#    MM     `Mb MM    dP    dM'       ` dM'      `MM M YMb   M     VM:   ,V     MM   d      MM   ,M9  `MMb.       MM dM'      `MM M YMb   M `MMb.      #
#    MM      MM MM"""bg.    MM          MM        MM M  `MN. M      MM.  M'     MMmmMM      MMmmdM9     `YMMNq.   MM MM        MM M  `MN. M   `YMMNq.  #
#    MM     ,MP MM    `Y    MM.         MM.      ,MP M   `MM.M      `MM A'      MM   Y  ,   MM  YM.   .     `MM   MM MM.      ,MP M   `MM.M .     `MM  #
#    MM    ,dP' MM    ,9    `Mb.     ,' `Mb.    ,dP' M     YMM       :MM;       MM     ,M   MM   `Mb. Mb     dM   MM `Mb.    ,dP' M     YMM Mb     dM  #
#  .JMMmmmdP' .JMMmmmd9       `"bmmmd'    `"bmmd"' .JML.    YM        VF      .JMMmmmmMMM .JMML. .JMM.P"Ybmmd"  .JMML. `"bmmd"' .JML.    YM P"Ybmmd"   #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#--------------------------------------------------#
# 001 : KEGG -> CAS, KEGG -> ChEBI, KEGG -> PubChem ID
# Use saved KEGG_DBlinks, Link KEGG ID to CAS No., ChEBI ID, and PubChem ID.
# format: dict( "KEGG_id" : [ ("CAS", ), ("ChEBI", ), ("PubChem", ), ] )
pickle_in1 = open("./KEGGScrapSavings/CPD02_KEGG_CAS_ChEBI_PubChem_dict.p","rb")
KEGG_CAS_ChEBI_PubChem_dict = pickle.load(pickle_in1)
pickle_in1.close()
#print (KEGG_CAS_ChEBI_PubChem_dict)

#--------------------------------------------------#
# 002 : ChEBI -> SMILES
# Use bioservices.ChEBI(), to convert ChEBI to SMILES.
def ChEBI_to_SMILES(ChEBI_a):
    ch = ChEBI()
    res = ch.getCompleteEntity("CHEBI:" + ChEBI_a)
    return res.smiles

#--------------------------------------------------#
# 003 : SMILES -> PubChem ID
# Use pubchempy to convert smiles to PubChem_id???
def smiles_to_PubChem_id(SMILES_a):
    # Example SMILES_a='Cn1cnc2n(C)c(=O)n(C)c(=O)c12'
    a = pubchempy.get_compounds(SMILES_a, 'smiles')
    return a[0]

#--------------------------------------------------#
# 004 : PubChem_id -> properties_str_list
# Use pubchempy, convert PubChem_id to a list of molecule info
def PubChem_id_to_info_dict(PubChem_id):
    c = pubchempy.Compound.from_cid(PubChem_id)
    return c.to_dict().keys()
    # return c.to_dict(properties=properties_str.lower().split(', ')) # cannot remember why return this before

properties_str_list=['bond_stereo_count'            , 
                     'defined_atom_stereo_count'    , 
                     'mmff94_partial_charges_3d'    , 
                     'inchi'                        , 
                     'h_bond_donor_count'           , 
                     'feature_selfoverlap_3d'       , 
                     'canonical_smiles'             , 
                     'shape_fingerprint_3d'         , 
                     'isotope_atom_count'           , 
                     'molecular_weight'             , 
                     'coordinate_type'              , 
                     'charge'                       , 
                     'conformer_rmsd_3d'            , 
                     'isomeric_smiles'              , 
                     'exact_mass'                   , 
                     'rotatable_bond_count'         , 
                     'xlogp'                        , 
                     'defined_bond_stereo_count'    , 
                     'iupac_name'                   , 
                     'monoisotopic_mass'            , 
                     'tpsa'                         , 
                     'volume_3d'                    , 
                     'inchikey'                     , 
                     'elements'                     , 
                     'bonds'                        , 
                     'mmff94_energy_3d'             , 
                     'conformer_id_3d'              , 
                     'atoms'                        , 
                     'fingerprint'                  , 
                     'covalent_unit_count'          , 
                     'shape_selfoverlap_3d'         , 
                     'undefined_atom_stereo_count'  , 
                     'cid'                          , 
                     'cactvs_fingerprint'           , 
                     'pharmacophore_features_3d'    , 
                     'effective_rotor_count_3d'     , 
                     'record'                       , 
                     'complexity'                   , 
                     'heavy_atom_count'             , 
                     'undefined_bond_stereo_count'  , 
                     'h_bond_acceptor_count'        , 
                     'molecular_formula'            ,  
                     'atom_stereo_count'            , 
                     'multipoles_3d'
                     ]



#--------------------------------------------------#
# 005 - pubchempy, converting PubChem_id to SMILES
def PubChem_id_to_SMILES(PubChem_id):
    c = pubchempy.Compound.from_cid(PubChem_id)
    info_dict = c.to_dict()
    return info_dict['canonical_smiles']
    # return c.to_dict(properties=properties_str.lower().split(', ')) # cannot remember why return this before

#print PubChem_id_to_SMILES('1234')

#------------------------------
# 005-1 - pubchempy, converting SMILES to canonical, non-big-pi SMILES
def SMILES_to_PubChem_canonical_SMILES_(SMILES_a):
    temp_a = pubchempy.get_properties('CanonicalSMILES', SMILES_a , 'smiles')
    a= temp_a[0]['CanonicalSMILES']
    return a.encode('ascii','ignore')


#--------------------------------------------------#
# 006 - cirpy, converting CAS to SMILES
def CAS_to_SMILES(CAS_a):
    return cirpy.resolve(CAS_a, 'smiles')


#--------------------------------------------------#
# 007, Use 001 and 002 to convert KEGG_id to SMILES
def KEGG_id_to_SMILES_by_ChEBI(KEGG_id_a):

    ChEBI_a = KEGG_CAS_ChEBI_PubChem_dict[KEGG_id_a][1][0]
    return ChEBI_to_SMILES(ChEBI_a)

#--------------------------------------------------#
# 008, Use 001 and 006 to convert KEGG_id to SMILES
def KEGG_id_to_SMILES_by_CAS(KEGG_id_a):
    
    CAS_a = KEGG_CAS_ChEBI_PubChem_dict[KEGG_id_a][0][0]
    return CAS_to_SMILES(CAS_a)

#--------------------------------------------------#
# 009 KEGG id to compound name. (Getting the first name in the KEGG DB)
def KEGG_id_to_cmpd_name(one_KEGG_id):
    KEGG_cmpd_info_bioservices=str(KEGG().get(one_KEGG_id))
    kcib_list = KEGG_cmpd_info_bioservices.split('\n')

    if KEGG_cmpd_info_bioservices.find("NAME")!=-1:
        for one_row in kcib_list:
            if one_row.find("NAME")!=-1:
                str_info=one_row.replace("NAME",'')
                str_info=str_info.replace(" ",'')
                nme = str_info.replace(";",'')
                break
    else:
        nme = "None"

    return nme




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#   `7MMF'  `7MMF'     db     `7MM"""Mq. `7MM"""Yb.     .g8"""bgd   .g8""8q. `7MM"""Yb. `7MMF'`7MN.   `7MF' .g8"""bgd                                  #
#     MM      MM      ;MM:      MM   `MM.  MM    `Yb. .dP'     `M .dP'    `YM. MM    `Yb. MM    MMN.    M .dP'     `M                                  #
#     MM      MM     ,V^MM.     MM   ,M9   MM     `Mb dM'       ` dM'      `MM MM     `Mb MM    M YMb   M dM'       `                                  #
#     MMmmmmmmMM    ,M  `MM     MMmmdM9    MM      MM MM          MM        MM MM      MM MM    M  `MN. M MM                                           #
#     MM      MM    AbmmmqMA    MM  YM.    MM     ,MP MM.         MM.      ,MP MM     ,MP MM    M   `MM.M MM.    `7MMF'                                #
#     MM      MM   A'     VML   MM   `Mb.  MM    ,dP' `Mb.     ,' `Mb.    ,dP' MM    ,dP' MM    M     YMM `Mb.     MM                                  #
#   .JMML.  .JMML.AMA.   .AMMA.JMML. .JMM.JMMmmmdP'     `"bmmmd'    `"bmmd"' .JMMmmmdP' .JMML..JML.    YM   `"bmmmdPY                                  #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

def smiles_hard_coding_process(one_smiles, skip = False):
    if skip == True:
        return one_smiles

    smiles_processed=max(one_smiles.split("."), key = len)
    smiles_processed=smiles_processed.replace("([*])", "")
    smiles_processed=smiles_processed.replace("[*]", "")
    smiles_processed=smiles_processed.replace("[O-]", "O")
    smiles_processed=smiles_processed.replace("[N+]", "N")
    smiles_processed=smiles_processed.replace("[NH+]", "N")
    smiles_processed=smiles_processed.replace("[NH2+]", "N")
    smiles_processed=smiles_processed.replace("[NH3+]", "N")
    smiles_processed=smiles_processed.replace("[N-]", "N")
    #smiles_processed=smiles_processed.replace("*-", "N")
    #smiles_processed=smiles_processed.replace("-*", "N")
    #smiles_processed=smiles_processed.replace("(-*)", "N")
    smiles_processed=smiles_processed.replace("*", "")
    return smiles_processed

def get_longest_smiles_from_aggregates(one_smiles_string):
    return max(one_smiles_string.split("."), key = len)


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    .M"""bgd   .g8"""bgd `7MM"""Mq.        db      `7MM"""Mq.         `7MMF' `YMM' `7MM"""YMM    .g8"""bgd    .g8"""bgd                               #
#   ,MI    "Y .dP'     `M   MM   `MM.      ;MM:       MM   `MM.          MM   .M'     MM    `7  .dP'     `M  .dP'     `M                               #
#   `MMb.     dM'       `   MM   ,M9      ,V^MM.      MM   ,M9           MM .d"       MM   d    dM'       `  dM'       `                               #
#     `YMMNq. MM            MMmmdM9      ,M  `MM      MMmmdM9            MMMMM.       MMmmMM    MM           MM                                        #
#   .     `MM MM.           MM  YM.      AbmmmqMA     MM                 MM  VMA      MM   Y  , MM.    `7MMF'MM.    `7MMF'                             #
#   Mb     dM `Mb.     ,'   MM   `Mb.   A'     VML    MM                 MM   `MM.    MM     ,M `Mb.     MM  `Mb.     MM                               #
#   P"Ybmmd"    `"bmmmd'  .JMML. .JMM..AMA.   .AMMA..JMML.             .JMML.   MMb..JMMmmmmMMM   `"bmmmdPY    `"bmmmdPY                               #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

def initialize_KEGG_id():
    # All_KEGG_IDs_list: ["C00001", ...]
    return ["C" + str(cpd_idx+100000).replace("1", "", 1)  for cpd_idx in range(23000)[1:]]







#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#   `7MM"""Yb.  `7MM"""Yp,      `7MM"""YMM  `7MN.   `7MF'MMP""MM""YMM `7MM"""Mq.  `7MMF'`7MM"""YMM   .M"""bgd                                          #
#     MM    `Yb.  MM    Yb        MM    `7    MMN.    M  P'   MM   `7   MM   `MM.   MM    MM    `7  ,MI    "Y                                          #
#     MM     `Mb  MM    dP        MM   d      M YMb   M       MM        MM   ,M9    MM    MM   d    `MMb.                                              #
#     MM      MM  MM"""bg.        MMmmMM      M  `MN. M       MM        MMmmdM9     MM    MMmmMM      `YMMNq.                                          #
#     MM     ,MP  MM    `Y        MM   Y  ,   M   `MM.M       MM        MM  YM.     MM    MM   Y  , .     `MM                                          #
#     MM    ,dP'  MM    ,9        MM     ,M   M     YMM       MM        MM   `Mb.   MM    MM     ,M Mb     dM                                          #
#   .JMMmmmdP'  .JMMmmmd9       .JMMmmmmMMM .JML.    YM     .JMML.    .JMML. .JMM..JMML..JMMmmmmMMM P"Ybmmd"                                           #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

# Script for scrapping the CPD02_KEGG_CAS_ChEBI_PubChem_dict. 
def parse_KEGG_text_API(KEGG_id, attribute): 
# to (mainly) get the CAS and ChEBI of the compound with known KEGG_id
# attribute can only be one string type variable in the list below,
# [FORMULA,EXACT_MASS,MOL_WEIGHT,CAS,PubChem,ChEBI,ChEMBL,KNApSAcK,PDB-CCD,3DMET,NIKKAJI,ATOM]
# NAME,REACTION,PATHWAY,ENZYME will not be parsed for now

    KEGG_API_url = "http://rest.kegg.jp/get/"+KEGG_id
    html = urllib.request.urlopen(KEGG_API_url)
    KEGG_cmpd_info=html.read()
    if KEGG_cmpd_info=='':
        print ('Cannot find this KEGG_id!')
        return 'Not Exist'
    KEGG_cmpd_info_list=KEGG_cmpd_info.decode('utf-8').split('\n')

    Assert(attribute in ['FORMULA','EXACT_MASS','MOL_WEIGHT','CAS','PubChem','ChEBI','ChEMBL','KNApSAcK','PDB-CCD','3DMET','NIKKAJI','ATOM'], f"attribute Error!")

    if attribute == 'CAS':
        attribute = 'DBLINKS     CAS:'
    if attribute in ['PubChem','ChEBI','ChEMBL','KNApSAcK','PDB-CCD','3DMET','NIKKAJI']:
        attribute=attribute+':'


    for str_info in KEGG_cmpd_info_list:
        if str_info.find(attribute)!=-1:
            str_info = str_info.replace(attribute,'')
            str_info = str_info.replace('DBLINKS','')
            str_info_set = set(str_info.split(' '))
            str_info_set.discard('')
            str_info_tuple = tuple(str_info_set)
                

            return str_info_tuple
    #print 'Cannot find this attribute!'
    return 'None'

#====================================================================================================#
def generate_CPD02_KEGG_CAS_ChEBI_PubChem_dict(folder = "./", filename = "KEGG_CAS_ChEBI_PubChem_dict.p"):
    def parse_KEGG_DBlinks(KEGG_id):
        # Use parse_KEGG_text_API to get CAS, ChEBI and PubChem ID.
        DBlinks=[]
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'CAS'))
        if DBlinks[0]=='Not Exist':
            return DBlinks
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'ChEBI'))
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'PubChem'))
        return DBlinks

    #--------------------------------------------------#
    # Use the following code to generate the dictionary if the .p file doesnt exist.
    '''
    All_KEGG_IDs_list=initialize_KEGG_id()
    CPD02_KEGG_CAS_ChEBI_PubChem_dict=dict([])
    for KEGG_id in All_KEGG_IDs_list:
        DBlinks=parse_KEGG_DBlinks(KEGG_id)
        try:
            DBlinks=parse_KEGG_DBlinks(KEGG_id)
        except Exception as error_message:
            print('Error raised when parsing', KEGG_id, ': ', error_message, '!!!')
            continue
        if DBlinks==['Not Exist']:
            continue
        CPD02_KEGG_CAS_ChEBI_PubChem_dict[KEGG_id]=DBlinks
        print(KEGG_id, DBlinks)
    print(CPD02_KEGG_CAS_ChEBI_PubChem_dict)
    
    pickle_out1 = open(Path(folder) / filename,"wb")
    pickle.dump(CPD02_KEGG_CAS_ChEBI_PubChem_dict, pickle_out1)
    pickle_out1.close()
    
    '''
    #--------------------------------------------------#
    pickle_in1 = open(Path(folder) / filename,"rb")
    CPD02_KEGG_CAS_ChEBI_PubChem_dict = pickle.load(pickle_in1)
    pickle_in1.close()


    return CPD02_KEGG_CAS_ChEBI_PubChem_dict




#====================================================================================================#
def update_CPD02_KEGG_CAS_ChEBI_PubChem_dict(folder = "./", filename = "KEGG_CAS_ChEBI_PubChem_dict.p"):
    def parse_KEGG_DBlinks(KEGG_id):
        # Use parse_KEGG_text_API to get CAS, ChEBI and PubChem ID.
        DBlinks=[]
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'CAS'))
        if DBlinks[0]=='Not Exist':
            return DBlinks
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'ChEBI'))
        DBlinks.append(parse_KEGG_text_API(KEGG_id,'PubChem'))
        return DBlinks

    #--------------------------------------------------#
    pickle_in1 = open(Path(folder) / filename,"rb")
    CPD02_KEGG_CAS_ChEBI_PubChem_dict = pickle.load(pickle_in1)
    pickle_in1.close()

    #--------------------------------------------------#
    # Use the following code to generate the dictionary if the .p file doesnt exist.
    
    Update_KEGG_IDs_list = ["C07138" ,"C21719" ,"C21954" ,"C21889" ,"C21889" ,"C21737" ,"C21896" ,"C21895" ,"C22029" ,"C22127" ,\
                            "C22128" ,"C21649" ,"C21649" ,"C21593" ,"C22031" ,"C21785" ,"C21669" ,"C22034" ,"C22034" ,"C22036" ,\
                            "C22036" ,"C22181" ,"C22180" ,"C22182" ,"C21775" ,"C21858" ,"C21858" ,"C22033" ,"C22084" ,"C22224" ]

    for KEGG_id in Update_KEGG_IDs_list:
        DBlinks=parse_KEGG_DBlinks(KEGG_id)
        try:
            DBlinks=parse_KEGG_DBlinks(KEGG_id)
        except Exception as error_message:
            print('Error raised when parsing', KEGG_id, ': ', error_message, '!!!')
            continue
        if DBlinks==['Not Exist']:
            continue
        CPD02_KEGG_CAS_ChEBI_PubChem_dict[KEGG_id]=DBlinks
        print(KEGG_id, DBlinks)
    print(CPD02_KEGG_CAS_ChEBI_PubChem_dict)
    
    pickle_out1 = open(Path(folder) / filename,"wb")
    pickle.dump(CPD02_KEGG_CAS_ChEBI_PubChem_dict, pickle_out1)
    pickle_out1.close()
    
    #--------------------------------------------------#
    pickle_in1 = open(Path(folder) / filename,"rb")
    CPD02_KEGG_CAS_ChEBI_PubChem_dict = pickle.load(pickle_in1)
    pickle_in1.close()


    return CPD02_KEGG_CAS_ChEBI_PubChem_dict





#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    `7MMF' `YMM'`7MM"""YMM    .g8"""bgd    .g8"""bgd                         .M"""bgd `7MMM.     ,MMF'`7MMF'`7MMF'      `7MM"""YMM   .M"""bgd         #
#      MM   .M'    MM    `7  .dP'     `M  .dP'     `M            `MM.        ,MI    "Y   MMMb    dPMM    MM    MM          MM    `7  ,MI    "Y         #
#      MM .d"      MM   d    dM'       `  dM'       `              `Mb.      `MMb.       M YM   ,M MM    MM    MM          MM   d    `MMb.             #
#      MMMMM.      MMmmMM    MM           MM                MMMMMMMMMMMM;      `YMMNq.   M  Mb  M' MM    MM    MM          MMmmMM      `YMMNq.         #
#      MM  VMA     MM   Y  , MM.    `7MMF'MM.    `7MMF'             ,M'      .     `MM   M  YM.P'  MM    MM    MM      ,   MM   Y  , .     `MM         #
#      MM   `MM.   MM     ,M `Mb.     MM  `Mb.     MM             .M'        Mb     dM   M  `YM'   MM    MM    MM     ,M   MM     ,M Mb     dM         #
#    .JMML.   MMb.JMMmmmmMMM   `"bmmmdPY    `"bmmmdPY                        P"Ybmmd"  .JML. `'  .JMML..JMML..JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"          #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

def generate_KEGG_id_SMILES_by_ChEBI_CAS_csv(folder = "./", filename = "KEGG_id_SMILES_by_ChEBI_CAS.csv"): # Step 01

    # Write a csv
    text_file = open(Path(folder) / filename, "w")

    for one_KEGGid in list(KEGG_CAS_ChEBI_PubChem_dict.keys()):
        #--------------------------------------------------#
        try:
            smiles_by_ChEBI = KEGG_id_to_SMILES_by_ChEBI(one_KEGGid)
        except Exception:
            smiles_by_ChEBI = "None" 
        #--------------------------------------------------#
        try:
            smiles_by_CAS = KEGG_id_to_SMILES_by_CAS(one_KEGGid)
        except Exception:
            smiles_by_CAS = "None" 
        #--------------------------------------------------#
        print (str(smiles_by_ChEBI) + "\n" + str(smiles_by_CAS) + ", " + one_KEGGid)

        text_file.write(one_KEGGid)
        text_file.write(",")
        text_file.write(str(smiles_by_ChEBI))
        text_file.write(",")
        text_file.write(str(smiles_by_CAS))
        text_file.write("\n")

    text_file.close()
    
    return




#====================================================================================================#
def update_KEGG_id_SMILES_by_ChEBI_CAS_csv(folder = "./", filename = "KEGG_id_SMILES_by_ChEBI_CAS.csv"): # Step 01

    # Write a csv
    text_file = open(Path(folder) / filename, "a")

    Update_KEGG_IDs_list = ["C07138" ,"C21719" ,"C21954" ,"C21889" ,"C21889" ,"C21737" ,"C21896" ,"C21895" ,"C22029" ,"C22127" ,\
                            "C22128" ,"C21649" ,"C21649" ,"C21593" ,"C22031" ,"C21785" ,"C21669" ,"C22034" ,"C22034" ,"C22036" ,\
                            "C22036" ,"C22181" ,"C22180" ,"C22182" ,"C21775" ,"C21858" ,"C21858" ,"C22033" ,"C22084" ,"C22224" ]

    for one_KEGGid in Update_KEGG_IDs_list:
        #--------------------------------------------------#
        try:
            smiles_by_ChEBI = KEGG_id_to_SMILES_by_ChEBI(one_KEGGid)
        except Exception:
            smiles_by_ChEBI = "None" 
        #--------------------------------------------------#
        try:
            smiles_by_CAS = KEGG_id_to_SMILES_by_CAS(one_KEGGid)
            smiles_by_CAS = smiles_by_CAS if smiles_by_CAS != "N" else "None"
        except Exception:
            smiles_by_CAS = "None" 
        #--------------------------------------------------#
        print (str(smiles_by_ChEBI) + "\n" + str(smiles_by_CAS) + ", " + one_KEGGid)

        text_file.write(one_KEGGid)
        text_file.write(",")
        text_file.write(str(smiles_by_ChEBI))
        text_file.write(",")
        text_file.write(str(smiles_by_CAS))
        text_file.write("\n")

    text_file.close()
    
    return





#====================================================================================================#
def generate_KEGG_id_SMILES_by_mol_file(folder = "./", filename = "KEGG_id_SMILES_by_mol_file.csv"):

    All_KEGG_IDs_list = initialize_KEGG_id()

    # Write a csv
    text_file = open(Path(folder) / filename, "w")
    text_file.write("KEGGid,SMILES_by_mol\n")

    for one_KEGGid in All_KEGG_IDs_list:
        if os.path.exists("./mol-files/" + one_KEGGid + ".mol"):
            mol_x = Chem.MolFromMolFile("./mol-files/" + one_KEGGid + ".mol")
            smiles_x  = MolToSmiles_ZX(mol_x, bad_ss_dict = {}, 
                                        isomericSmiles     = True  ,
                                        kekuleSmiles       = False ,  
                                        canonical          = False  ,
                                        SMARTS_bool        = False , )
            print(one_KEGGid, smiles_x, )
        
            text_file.write(one_KEGGid)
            text_file.write(",")
            text_file.write(str(smiles_x))
            text_file.write("\n")

    text_file.close()
    return

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    `7MMF' `YMM'`7MM"""YMM    .g8"""bgd    .g8"""bgd                                    db     `7MMF'      `7MMF'                                     #
#      MM   .M'    MM    `7  .dP'     `M  .dP'     `M              `MM.                 ;MM:      MM          MM                                       #
#      MM .d"      MM   d    dM'       `  dM'       `                `Mb.              ,V^MM.     MM          MM                                       #
#      MMMMM.      MMmmMM    MM           MM                  MMMMMMMMMMMMD           ,M  `MM     MM          MM                                       #
#      MM  VMA     MM   Y  , MM.    `7MMF'MM.    `7MMF'               ,M'             AbmmmqMA    MM      ,   MM      ,                                #
#      MM   `MM.   MM     ,M `Mb.     MM  `Mb.     MM               .M'              A'     VML   MM     ,M   MM     ,M                                #
#    .JMML.   MMb.JMMmmmmMMM   `"bmmmdPY    `"bmmmdPY                              .AMA.   .AMMA.JMMmmmmMMM .JMMmmmmMMM                                #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


def KEGGid_to_ALL(folder = "./KEGGScrapSavings", filename = "CPD02_ALL_KEGG_DATA_DICT.p"): # Step 03


    All_KEGG_IDs_list = initialize_KEGG_id()
    ALL_KEGG_DATA_DICT = dict([]) # {"KEGGid" : dict([]), "KEGGid" : dict([]), ... }
    bad_KEGG_id_list = []
    for one_id in All_KEGG_IDs_list:
        print(one_id)
        ALL_KEGG_DATA_DICT[one_id] = dict([])
        try: 
            #--------------------------------------------------#
            # 
            KEGG_cmpd_info_bioservices = str(KEGG().get(one_id))
            KEGG_info_lines = KEGG_cmpd_info_bioservices.split('\n')
            text_file = open(Path(folder) / "KEGG_DB" / (one_id + "_all_data.txt"), "w")
            text_file.write(KEGG_cmpd_info_bioservices)
            text_file.close()

            #--------------------------------------------------#
            one_key = ""
            for one_line in KEGG_info_lines:
                if one_line != "":
                    if one_line[0] != " ":
                        one_key = one_line.split(" ")[0]
                        values_part = one_line.replace(one_key, "", 1)
                        values_part_list = values_part.split(" ")
                        values_part_list = [x for x in values_part_list if x != ""]
                        if values_part_list != []:
                            first_char_idx = values_part.find(values_part_list[0][0])
                            if first_char_idx != -1:
                                values_part = values_part.replace(" " * first_char_idx, "", 1)
                            else:
                                values_part = ""
                        else:
                            values_part = ""
                        ALL_KEGG_DATA_DICT[one_id][one_key] = []
                        ALL_KEGG_DATA_DICT[one_id][one_key].append(values_part)
                    else:
                        one_line_list = one_line.split(" ")
                        one_line_list = [x for x in one_line_list if x != ""]
                        first_char_idx = one_line.find(one_line_list[0][0])

                        if first_char_idx != -1:
                            values_part = one_line.replace(" " * first_char_idx, "", 1)
                        else:
                            values_part = ""
                        ALL_KEGG_DATA_DICT[one_id][one_key].append(values_part)
            print("Successfully parse ", one_id)
        except Exception:
            bad_KEGG_id_list.append(one_id)
            print( "Error returned when parsing ", one_id, str(KEGG().get(one_id)) )

    print(ALL_KEGG_DATA_DICT)

    pickle_out1 = open(Path(folder) / filename,"wb")
    pickle.dump(ALL_KEGG_DATA_DICT, pickle_out1)
    pickle_out1.close()

    return 

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    `7MMF' `YMM'`7MM"""YMM    .g8"""bgd    .g8"""bgd                             `7MN.   `7MF'      db      `7MMM.     ,MMF'`7MM"""YMM                #
#      MM   .M'    MM    `7  .dP'     `M  .dP'     `M              `MM.             MMN.    M       ;MM:       MMMb    dPMM    MM    `7                #
#      MM .d"      MM   d    dM'       `  dM'       `                `Mb.           M YMb   M      ,V^MM.      M YM   ,M MM    MM   d                  #
#      MMMMM.      MMmmMM    MM           MM                  MMMMMMMMMMMMD         M  `MN. M     ,M  `MM      M  Mb  M' MM    MMmmMM                  #
#      MM  VMA     MM   Y  , MM.    `7MMF'MM.    `7MMF'               ,M'           M   `MM.M     AbmmmqMA     M  YM.P'  MM    MM   Y  ,               #
#      MM   `MM.   MM     ,M `Mb.     MM  `Mb.     MM               .M'             M     YMM    A'     VML    M  `YM'   MM    MM     ,M               #
#    .JMML.   MMb.JMMmmmmMMM   `"bmmmdPY    `"bmmmdPY                             .JML.    YM  .AMA.   .AMMA..JML. `'  .JMML..JMMmmmmMMM               #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


































#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    `7MM"""Mq. `7MM"""YMM        db     `7MM"""Yb.       .M"""bgd `7MMM.     ,MMF'`7MMF'`7MMF'      `7MM"""YMM   .M"""bgd                             #
#      MM   `MM.  MM    `7       ;MM:      MM    `Yb.    ,MI    "Y   MMMb    dPMM    MM    MM          MM    `7  ,MI    "Y                             #
#      MM   ,M9   MM   d        ,V^MM.     MM     `Mb    `MMb.       M YM   ,M MM    MM    MM          MM   d    `MMb.                                 #
#      MMmmdM9    MMmmMM       ,M  `MM     MM      MM      `YMMNq.   M  Mb  M' MM    MM    MM          MMmmMM      `YMMNq.                             #
#      MM  YM.    MM   Y  ,    AbmmmqMA    MM     ,MP    .     `MM   M  YM.P'  MM    MM    MM      ,   MM   Y  , .     `MM                             #
#      MM   `Mb.  MM     ,M   A'     VML   MM    ,dP'    Mb     dM   M  `YM'   MM    MM    MM     ,M   MM     ,M Mb     dM                             #
#    .JMML. .JMM.JMMmmmmMMM .AMA.   .AMMA.JMMmmmdP'      P"Ybmmd"  .JML. `'  .JMML..JMML..JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"                              #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Have not test this yet.

def Test_():

    def read_smiles(dictionary, one_smiles):
    # Test if the dictionary contains certain SMILES strings.
        try:    
            print ([unis(one_smiles),] + dictionary[ unis( one_smiles ) ]," ,")
        except Exception:
            print ("!!!!! Cannot Find : ", unis(one_smiles), " !!!!!" )
        return 

    pickle_in1=open("./KEGGScrapSavings/KEGG_nme_canonical_SMILES_dict.p","rb")
    KEGG_nme_canonical_SMILES_dict=pickle.load(pickle_in1)
    pickle_in1.close()
    d1=KEGG_nme_canonical_SMILES_dict
    #print unis("CC(=O)I")
    




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#    `7MM"""YMM  `YMM'   `MP'       db      `7MMM.     ,MMF'`7MM"""Mq. `7MMF'      `7MM"""YMM   .M"""bgd                                               #
#      MM    `7    VMb.  ,P        ;MM:       MMMb    dPMM    MM   `MM.  MM          MM    `7  ,MI    "Y                                               #
#      MM   d       `MM.M'        ,V^MM.      M YM   ,M MM    MM   ,M9   MM          MM   d    `MMb.                                                   #
#      MMmmMM         MMb        ,M  `MM      M  Mb  M' MM    MMmmdM9    MM          MMmmMM      `YMMNq.                                               #
#      MM   Y  ,    ,M'`Mb.      AbmmmqMA     M  YM.P'  MM    MM         MM      ,   MM   Y  , .     `MM                                               #
#      MM     ,M   ,P   `MM.    A'     VML    M  `YM'   MM    MM         MM     ,M   MM     ,M Mb     dM                                               #
#    .JMMmmmmMMM .MM:.  .:MMa..AMA.   .AMMA..JML. `'  .JMML..JMML.     .JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"                                                #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

if __name__ == '__main__':

    mol_x = Chem.MolFromMolFile("./mol-files/" + "C02202" + ".mol")
    smiles_x = MolToSmiles_ZX(mol_x)
    print(smiles_x)


    #====================================================================================================#
    # Test 001 : KEGG -> CAS, KEGG -> ChEBI, KEGG -> PubChem ID

    print("\nTest KEGG_CAS_ChEBI_PubChem_dict.p: " )

    KEGG_CAS_ChEBI_PubChem_dict_new = \
        generate_CPD02_KEGG_CAS_ChEBI_PubChem_dict(folder   = "./KEGGScrapSavings", 
                                                   filename = "CPD02_KEGG_CAS_ChEBI_PubChem_dict.p")

    print("len(KEGG_CAS_ChEBI_PubChem_dict.keys()): ", len(KEGG_CAS_ChEBI_PubChem_dict.keys()))

    print("KEGG_CAS_ChEBI_PubChem_dict[\"C00001\"]: ", KEGG_CAS_ChEBI_PubChem_dict["C00001"])
    print("KEGG_CAS_ChEBI_PubChem_dict[\"C00002\"]: ", KEGG_CAS_ChEBI_PubChem_dict["C00002"])

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Test 001 : KEGG -> CAS, KEGG -> ChEBI, KEGG -> PubChem ID

    # print("\nTest KEGG_CAS_ChEBI_PubChem_dict.p: " )
    # KEGG_CAS_ChEBI_PubChem_dict_new = \
    #     update_CPD02_KEGG_CAS_ChEBI_PubChem_dict(folder   = "./KEGGScrapSavings", 
    #                                                filename = "CPD02_KEGG_CAS_ChEBI_PubChem_dict.p")


    #====================================================================================================#
    # Test 007, Use 001 and 002 to convert KEGG_id to SMILES
    print("\nTest KEGG_id_to_SMILES_by_ChEBI (007, Use 001 and 002 to convert KEGG_id to SMILES): ")

    print("KEGG_id_to_SMILES_by_ChEBI(\"C00005\"): ", unis(KEGG_id_to_SMILES_by_ChEBI("C00005")))
    print("KEGG_id_to_SMILES_by_ChEBI(\"C00006\"): ", unis(KEGG_id_to_SMILES_by_ChEBI("C00006")))
    print("KEGG_id_to_SMILES_by_ChEBI(\"C00007\"): ", unis(KEGG_id_to_SMILES_by_ChEBI("C00007")))


    #====================================================================================================#
    # Test 008, Use 001 and 006 to convert KEGG_id to SMILES
    print("\nTest KEGG_id_to_SMILES_by_CAS (007, Use 001 and 002 to convert KEGG_id to SMILES): ")

    print("KEGG_id_to_SMILES_by_CAS(\"C00005\"): ", unis(KEGG_id_to_SMILES_by_CAS("C00005")))
    print("KEGG_id_to_SMILES_by_CAS(\"C00006\"): ", unis(KEGG_id_to_SMILES_by_CAS("C00006")))
    print("KEGG_id_to_SMILES_by_CAS(\"C00007\"): ", unis(KEGG_id_to_SMILES_by_CAS("C00007")))

    #====================================================================================================#
    # Test 009, Convert KEGG id to compound name.

    # print (str(KEGG().get("C30105")))
    # print (str(KEGG().get("C00106")))
    # print (str(KEGG().get("C00107")))


    print("\nTest 009, Convert KEGG id to compound name: ")
    print("KEGG_id_to_cmpd_name(\"C00105\"): ", KEGG_id_to_cmpd_name("C00105"))
    print("KEGG_id_to_cmpd_name(\"C00106\"): ", KEGG_id_to_cmpd_name("C00106"))
    print("KEGG_id_to_cmpd_name(\"C00107\"): ", KEGG_id_to_cmpd_name("C00107"))




    #====================================================================================================#
    # Test generate_KEGG_id_SMILES_by_ChEBI_CAS_csv (Writing the csv file takes hours)
    print("\nTest generate_KEGG_id_SMILES_by_ChEBI_CAS_csv() : ")
    # Won't run, cuz it takes sooooo long.
    '''
    generate_KEGG_id_SMILES_by_ChEBI_CAS_csv(folder = "./KEGGScrapSavings", filename = "CPD02_KEGG_id_SMILES_by_ChEBI_CAS_new.csv")
    '''
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    print("\nTest update_KEGG_id_SMILES_by_ChEBI_CAS_csv() : ")
    # Won't run, cuz it will be adding duplicated lines to the file.
    '''
    update_KEGG_id_SMILES_by_ChEBI_CAS_csv(folder = "./KEGGScrapSavings", filename = "CPD02_KEGG_id_SMILES_by_ChEBI_CAS.csv")
    '''

    #====================================================================================================#
    # Test CPD02_KEGG_id_SMILES_by_ChEBI_CAS.csv (Writing the csv file takes hours)

    KEGG_id_SMILES_by_ChEBI_CAS_df = pd.read_csv("./KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_ChEBI_CAS.csv", header = 0)

    KEGGid_list     = KEGG_id_SMILES_by_ChEBI_CAS_df["KEGGid"].tolist()
    SMILES_by_ChEBI = KEGG_id_SMILES_by_ChEBI_CAS_df["SMILES_by_ChEBI"].tolist()
    SMILES_by_CAS   = KEGG_id_SMILES_by_ChEBI_CAS_df["SMILES_by_CAS"].tolist()

    KEGG_id_SMILES_by_ChEBI_CAS_dict = dict([])
    for i, id in enumerate(KEGGid_list):
        KEGG_id_SMILES_by_ChEBI_CAS_dict[id] = [SMILES_by_ChEBI[i], SMILES_by_CAS[i]]



    print("KEGG_id_SMILES_by_ChEBI_CAS_dict[\"C00105\"]: ", KEGG_id_SMILES_by_ChEBI_CAS_dict["C00105"])
    print("KEGG_id_SMILES_by_ChEBI_CAS_dict[\"C00106\"]: ", KEGG_id_SMILES_by_ChEBI_CAS_dict["C00106"])
    print("KEGG_id_SMILES_by_ChEBI_CAS_dict[\"C00107\"]: ", KEGG_id_SMILES_by_ChEBI_CAS_dict["C00107"])


    #====================================================================================================#
    # Test KEGGid_to_ALL
    # KEGGid_to_ALL(folder = "./KEGGScrapSavings", filename = "CPD02_ALL_KEGG_DATA_DICT.p")
    folder = Path("./KEGGScrapSavings")
    filename = "CPD02_ALL_KEGG_DATA_DICT_bk.p"

    pickle_in1 = open( folder / filename, "rb")
    ALL_KEGG_DATA_DICT = pickle.load(pickle_in1)
    pickle_in1.close()

    print(ALL_KEGG_DATA_DICT["C00001"]["NAME"])


    #====================================================================================================#
    # Test generate_KEGG_id_SMILES_by_mol_file
    # KEGGid_to_ALL(folder = "./KEGGScrapSavings", filename = "CPD02_ALL_KEGG_DATA_DICT.p")

    generate_KEGG_id_SMILES_by_mol_file(folder = "./KEGGScrapSavings", filename = "CPD02_KEGG_id_SMILES_by_mol_file_no_uniq.csv")






















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







#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#   `7MM"""Mq. `7MM"""YMM    .g8"""bgd `YMM'   `MM' .g8"""bgd `7MMF'      `7MMF'`7MN.   `7MF' .g8"""bgd                                                #
#     MM   `MM.  MM    `7  .dP'     `M   VMA   ,V .dP'     `M   MM          MM    MMN.    M .dP'     `M                                                #
#     MM   ,M9   MM   d    dM'       `    VMA ,V  dM'       `   MM          MM    M YMb   M dM'       `                                                #
#     MMmmdM9    MMmmMM    MM              VMMP   MM            MM          MM    M  `MN. M MM                                                         #
#     MM  YM.    MM   Y  , MM.              MM    MM.           MM      ,   MM    M   `MM.M MM.    `7MMF'                                              #
#     MM   `Mb.  MM     ,M `Mb.     ,'      MM    `Mb.     ,'   MM     ,M   MM    M     YMM `Mb.     MM                                                #
#   .JMML. .JMM.JMMmmmmMMM   `"bmmmd'     .JMML.    `"bmmmd'  .JMMmmmmMMM .JMML..JML.    YM   `"bmmmdPY                                                #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# 004 - ChemSpider, converting CS_id to mol files (2D & 3D)
'''
cs = ChemSpider("PW2qJbvGaAK1gAthkR6n4jDxoUJ3xw2D")
c1 = cs.get_compound(2424)  # Specify compound by ChemSpider ID
c2 = cs.search('benzene')  # S

info = cs.get_extended_compound_info(2424)
print info
mol=cs.get_record_mol(2424, calc3d=True)
print mol
'''
#--------------------------------------------------#
'''
c3=cs.search('CCCCC')
print 'caonima'
print c3
a=[]
b=[]
for result in c3:
    a.append(result.mol2D)
    b.append(result.mol3D)
print a
print b
'''
#--------------------------------------------------#
# 004 - 
'''
from bioservices import *
s = ChemSpider("PW2qJbvGaAK1gAthkR6n4jDxoUJ3xw2D")
s.find("Pyridine")
results = s.GetExtendedCompoundInfo(1020)
print results['averagemass']
'''
#--------------------------------------------------#
# 004 - 
'''
from chemspipy import *
cs = ChemSpider("PW2qJbvGaAK1gAthkR6n4jDxoUJ3xw2D")
c1 = cs.get_compound(2424)  # Specify compound by ChemSpider ID
c2 = cs.search('benzene')  # S
info = cs.get_extended_compound_info(2424)
print info
mol=cs.get_record_mol(2424, calc3d=True)
print mol
c3=cs.search('Cn1cnc2n(C)c(=O)n(C)c(=O)c12')
print 'caonima'
for result in c3:
    print(result.csid)
    '''




#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#   `7MM"""Yb. `7MM"""YMM  `7MM"""Mq.`7MM"""Mq. `7MM"""YMM    .g8"""bgd     db  MMP""MM""YMM `7MM"""YMM  `7MM"""Yb.                                    #
#     MM    `Yb. MM    `7    MM   `MM. MM   `MM.  MM    `7  .dP'     `M    ;MM: P'   MM   `7   MM    `7    MM    `Yb.                                  #
#     MM     `Mb MM   d      MM   ,M9  MM   ,M9   MM   d    dM'       `   ,V^MM.     MM        MM   d      MM     `Mb                                  #
#     MM      MM MMmmMM      MMmmdM9   MMmmdM9    MMmmMM    MM           ,M  `MM     MM        MMmmMM      MM      MM                                  #
#     MM     ,MP MM   Y  ,   MM        MM  YM.    MM   Y  , MM.          AbmmmqMA    MM        MM   Y  ,   MM     ,MP                                  #
#     MM    ,dP' MM     ,M   MM        MM   `Mb.  MM     ,M `Mb.     ,' A'     VML   MM        MM     ,M   MM    ,dP'                                  #
#   .JMMmmmdP' .JMMmmmmMMM .JMML.    .JMML. .JMM.JMMmmmmMMM   `"bmmmd'.AMA.   .AMMA.JMML.    .JMMmmmmMMM .JMMmmmdP'                                    #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# Reformatting old data files so that dont need hours to scrap data again. 
'''
smiles_df_1 = pd.read_csv("./KEGGScrapSavings/KEGG_id_SMILES_1_bk.csv", header = None)
smiles_df_1.columns = ["KEGGid", "SMILES"]

KEGGid_list_1 = smiles_df_1["KEGGid"].tolist()
SMILES_list_1 = smiles_df_1["SMILES"].tolist()

dict_1 = dict([])
for one_KEGG in KEGGid_list_1:
    dict_1[one_KEGG.replace(" ", "")] = SMILES_list_1[KEGGid_list_1.index(one_KEGG)].replace(" ", "")



smiles_df_2 = pd.read_csv("./KEGGScrapSavings/KEGG_id_SMILES_2_bk.csv", header = None)
smiles_df_2.columns = ["KEGGid", "SMILES"]

KEGGid_list_2 = smiles_df_2["KEGGid"].tolist()
SMILES_list_2 = smiles_df_2["SMILES"].tolist()

dict_2 = dict([])
for one_KEGG in KEGGid_list_2:
    dict_2[one_KEGG.replace(" ", "")] = SMILES_list_2[KEGGid_list_2.index(one_KEGG)].replace(" ", "")



KEGGid_all_list = sorted(list(set(KEGGid_list_1).union(set(KEGGid_list_2))))
print(KEGGid_all_list[:50])



text_file = open('./KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_ChEBI_CAS.csv', "w") 
text_file.write("KEGGid,SMILES_by_ChEBI,SMILES_by_CAS\n")
for one_KEGG_id in KEGGid_all_list:
    text_file.write(one_KEGG_id.replace(" ", ""))
    text_file.write(",")
    text_file.write("None" if dict_1[one_KEGG_id.replace(" ", "")] == "NAN" else dict_1[one_KEGG_id.replace(" ", "")]  )
    text_file.write(",")
    text_file.write("None" if dict_2[one_KEGG_id.replace(" ", "")] == "N" else dict_2[one_KEGG_id.replace(" ", "")]  )
    text_file.write("\n")




text_file = open('./KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_ChEBI.csv', "w") 
text_file.write("KEGGid,SMILES_by_ChEBI\n")
for one_KEGG_id in KEGGid_all_list:
    text_file.write(one_KEGG_id.replace(" ", ""))
    text_file.write(",")
    text_file.write("None" if dict_1[one_KEGG_id.replace(" ", "")] == "NAN" else dict_1[one_KEGG_id.replace(" ", "")]  )
    text_file.write("\n")



text_file = open('./KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_CAS.csv', "w") 
text_file.write("KEGGid,SMILES_by_CAS\n")
for one_KEGG_id in KEGGid_all_list:
    text_file.write(one_KEGG_id.replace(" ", ""))
    text_file.write(",")
    text_file.write("None" if dict_2[one_KEGG_id.replace(" ", "")] == "N" else dict_2[one_KEGG_id.replace(" ", "")]  )
    text_file.write("\n")
    '''