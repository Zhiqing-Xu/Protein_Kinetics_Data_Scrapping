# Protein kinetics data scrapping from BRENDA database.

This repository contains raw data of protein kinetics (Ki, kcat, KM and kcat/KM) downloaded from BRENDA as well as the entire process of cleaning the data for compound protein interaction learning.

## Data Processing:
### W00_Data_BRENDA_cmpd_smiles_dict.py 
Search all compound names found in the data file on pubchem and cactus and obtain two dictionaries of compound -> smiles.

### W01_Data_BRENDA_Prep.py
Process the data and output three different files for each dataset.
All the outputs are in [./X_DataProcessing/X00_enzyme_datasets_processed/](https://github.com/Zhiqing-Xu/protein_kinetics_data_scrapping/tree/main/X_DataProcessing/X00_enzyme_datasets_processed).

## Outputs: 
### Output #1: _wi_unip_avg_val_screened.csv
> smallest dataset, most safe and confident\
(1) Get data with UniProt IDs (wi_unip);\
(2) Groupby smiles & sequence;\
(3) Screen SD < 1.0.

### Output #2: _wi_wo_unip_avg_val_screened_0.csv
> pan-genomic approach, largest dataset, least safe or unambiguous\
(1) For data without UniProt IDs (wo_unip), Groupby smiles & sequence;\
(2) Screen SD < 1.0;\
(3) Combine the processed wo_unip and wi_unip (in Output #1).

### Output #3: _wi_wo_unip_avg_val_screened_0.csv
> include some data via pan-genomic approach, strictly screen the dataset to ensure reliability\
(1) combine wi_unip (known sequences) and wo_unip (assigned sequences);\
(2) groupby org, EC & cmpd;\
(3) Screen SD < 1.0;\
(4) Add back all wi_unip data (in Output #1).
