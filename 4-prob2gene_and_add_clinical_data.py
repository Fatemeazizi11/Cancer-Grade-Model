import re
import os
import sys
import time
import pickle 
import requests
import numpy as np
import pandas as pd
from ipywidgets import widgets
import GEOparse

MAIN_DIR = '.'
DATASETS = [
            'GSE11121',
            'GSE18864',
            'GSE20711',
            'GSE23593',
            'GSE27120',
            'GSE32646',
            'GSE36771', 
            'GSE42568', 
            'GSE50948', 
            'GSE5460', 
            'GSE11001',
            'GSE87007',
            'GSE88770',
            'GSE7390', 
            'GSE78958', 
            'GSE45255',
            'GSE61304',
            'GSE63471',
            'GSE21653',
            'GSE26639', 
            'GSE17907',
            'GSE10810',
            'GSE25066', 
            'GSE47109',
            'GSE95700',
            'GSE5327',
            'GSE48390',
            'GSE58984',
            'GSE103091',
            'GSE45827',
            'GSE65194', 
            'GSE1456',
            'GSE102484']

MAIN_DIR_DATASET = os.path.join(MAIN_DIR, 'dataset')

annot = pd.read_csv('prob2gene.csv')


probe_to_symbol = pd.Series(annot['V2'].values, index=annot['V1']).to_dict()df = pd.read_csv(f'{MAIN_DIR}/merged_COMBAT_rma.matrix.csv', index_col=0).transpose()    
df.shape



result = []
for id, row in df.iterrows():
    genes = {}
    
    # Iterate over columns
    for probe, value in row.items():
        if probe_to_symbol.get(probe):
            if probe_to_symbol[probe] in genes:
                genes[probe_to_symbol[probe]].append(float(value))  # if we alread have the gene in our dictionary
            else:
                genes[probe_to_symbol[probe]] = [float(value)] # create a new list for the values of this gene
                
        # Just for the case of probes with 'X' at the beginning        
        if probe_to_symbol.get(probe[1:]):
            if probe_to_symbol[probe[1:]] in genes:
                genes[probe_to_symbol[probe[1:]]].append(float(value))
            else:
                genes[probe_to_symbol[probe[1:]]] = [float(value)]
                
    res_genes = {'sample_id': id.split('.')[0]}
    for k, v in genes.items():
        res_genes[k] = sum(v) / float(len(v))

    result.append(res_genes)
    
df_genes = pd.DataFrame(result)
# df_genes.head(5)




# deleting AFF 
unwanted =df_genes.filter(regex='^AFFX')
df_genes.drop(unwanted, axis=1, inplace=True)

# deleting -at
unwanted =df_genes.filter(regex='_at')
df_genes.drop(unwanted, axis=1, inplace=True)
df_genes.head(5) 

unwanted = df_genes.filter(regex='---')
df_genes.drop(unwanted, axis=1, inplace=True)
df_genes.head(5) 



# selecting 'sample_id' as index
df_genes =df_genes.set_index('sample_id')
df_genes.head(5)



os.mkdir(os.path.join(MAIN_DIR, 'dataset'))
for dataset in DATASETS:
    os.mkdir(os.path.join(MAIN_DIR, 'dataset', dataset))

pd.options.display.max_columns = 10
# pd.options.display.max_rows = None

MAIN_DIR_DATASET = os.path.join(MAIN_DIR, 'dataset')

# Specify the GSE ID for the dataset you're interested in
gse_id = "GSE11121"  # Replace with your GSE ID

# Load the GSE dataset
gse = GEOparse.get_GEO(geo=gse_id)

# Extract clinical data (phenotype data) from the GSE object
clinical_data = gse.phenotype_data

clinical_data.head(500)

columns_of_interest = {
    "sample_id": "geo_accession",
    "tissue": "source_name_ch1",
    "grade": "characteristics_ch1.2.grade",
    "stage": "characteristics_ch1.4.Stage",
    "subtype": "characteristics_ch1.6.SUBTYPE",
    "tumor-size(mm)": "characteristics_ch1.3.tumor size",
    "age-at-diagnosis(yrs)": "characteristics_ch1.2.age at diagnosis",
    "death-event-time(yrs)": "characteristics_ch1.5.SURV_DEATH",
    "death-event": "characteristics_ch1.4.DEATH_BC",
    "er": "characteristics_ch1.5.er status",
    "pr": "characteristics_ch1.6.pgr status",
    "her2": "characteristics_ch1.3.her2",
    "t.rfs": "characteristics_ch1.2.SURV_RELAPSE",
    "e.rfs": "characteristics_ch1.1.RELAPSE",
    "t.tdm": "characteristics_ch1.3.mfs (days)",
    "e.tdm": "characteristics_ch1.7.event_metastasis",
}

# Create a new DataFrame to store selected columns
df_selected = pd.DataFrame()

# Loop through the desired columns
for desired_column, actual_column in columns_of_interest.items():
    if actual_column in clinical_data.columns:
        # If the actual column exists, add it to the new DataFrame
        df_selected[desired_column] = clinical_data[actual_column]
    else:
        # If the actual column does not exist, create the column and fill it with NaN
        df_selected[desired_column] = np.nan






output_csv_path = os.path.join(MAIN_DIR_DATASET, gse_id, f'{gse_id}_annotation.csv')
df_selected.to_csv(output_csv_path, index=False)


df_total = None
df_clinical_data = pd.DataFrame()

for dataset in DATASETS:
    print(f'Loading {dataset} cancer data...')
    
    # loading clinical data
    df_cancer_ann = pd.read_csv(f'{MAIN_DIR_DATASET}/{dataset}/{dataset}_annotation.csv') 
    #df_cancer_ann['Dataset'] = cancer
    #df_cancer_ann.rename(columns={'sample':'sample_id'}, inplace=True) 
    
    # select important annotations
    df_important_ann = df_cancer_ann[['sample_id', 'tissue', 'grade', 'stage', 'subtype', 'tumor-size(mm)',
                                      'age-at-diagnosis(yrs)', 'death-event-time(yrs)', 'death-event', 'er',
                                      'pr', 'her2', 't.rfs', 'e.rfs',  't.tdm', 'e.tdm']]
    
    # all clinical data
    # df_clinical_data = df_clinical_data.append(df_important_ann, ignore_index = True)
    df_clinical_data = pd.concat([df_clinical_data, df_important_ann], ignore_index=True)


df_clinical_data = df_clinical_data.replace({'subtype': {' Basal': 'Basal', 'lumB': 'LumB', 'lumA': 'LumA', ' LumA': 'LumA', ' LumB':'LumB', ' HER2':'HER2', 'Normal': -1, '-1':-1}})
df_clinical_data = df_clinical_data.replace({'subtype': {'Her2':'HER2', ' Her2': 'HER2'}})
df_clinical_data = df_clinical_data.replace({'subtype': {'Basal':3, 'HER2':2, 'LumB':1, 'LumA':0 }})



df_genes_copy = df_genes.copy()
def convert_sample_id(sample):
    return sample.split('_')[0]
    
df_genes_copy.reset_index(inplace=True)
df_genes_copy.columns
df_genes_copy['sample_id'] = df_genes_copy['sample_id'].apply(convert_sample_id)



cancer_GE_Clinical_data = pd.DataFrame()
cancer_GE_Clinical_data = df_genes_copy.merge(df_clinical_data, on=['sample_id'])
cancer_GE_Clinical_data.shape


# save data as dict 
with open(f'{MAIN_DIR}/breast_clinical_data.pkl', 'wb') as fp:
    pickle.dump(df_clinical_data, fp)


# merge clinical data and gene experession
# cancer_GE_Clinical_data = pd.DataFrame()
# cancer_GE_Clinical_data = df_genes.merge(df_clinical_data, on=['sample_id'])

# remove dublicated rowes
cancer_GE_Clinical_data = cancer_GE_Clinical_data.drop_duplicates(['sample_id'], keep='last')
cancer_GE_Clinical_data.shape


# save data as dict 
with open(f'{MAIN_DIR}/breast_merged_clinical_gene_data.pkl', 'wb') as fp:
    pickle.dump(cancer_GE_Clinical_data, fp)


df_all_data = pd.read_pickle(f'{MAIN_DIR}/breast_merged_clinical_gene_data.pkl') 
df_all_data.shape   


df_all_cancer = df_all_data[df_all_data['tissue'] == 1]
df_all_cancer.to_pickle('df_tumour_expression.pkl')
df_all_cancer.to_csv('df_tumour_expression.csv')


df_all_normal = df_all_data[df_all_data['tissue'] == 0]
df_all_normal.to_pickle('df_normal_sample_expression.pkl')
df_all_normal.to_csv('df_normal_sample_expression.csv')


df_total = pd.read_pickle(f'{MAIN_DIR}/breast_merged_clinical_gene_data.pkl')
df_total.drop([
    'tissue', 'stage', 'subtype', 'tumor-size(mm)',
                                      'age-at-diagnosis(yrs)', 'death-event-time(yrs)', 'death-event', 'er',
                                      'pr', 'her2', 't.rfs', 'e.rfs',  't.tdm', 'e.tdm'
], axis=1, inplace=True)
df_total.head(6000)


df_total = df_total.replace({"grade": {'1': 1}})
df_total = df_total.replace({"grade": {'2': 2}})
df_total = df_total.replace({"grade": {'3': 3}})
df_total = df_total.replace({"grade": {'--': np.nan}})


df_total = df_total.replace({"grade": {3.0: 3}})
df_total = df_total.replace({"grade": {4.0: 3}})
df_total = df_total.replace({"grade": {1.0: 1}})

df_total['grade'].value_counts()








