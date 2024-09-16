library(affy)
library(dplyr)
library(data.table)
library(R.utils)



# CEL files  
DATASETS = list('GSE11121',
                'GSE18864',
                'GSE20711',
                'GSE23593',
                'GSE27120',
                'GSE32646',
                'GSE36771',
                'GSE42568',
                'GSE50948',
                'GSE5460_1',
                'GSE5460_2',
                'GSE5460_3',
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
                'GSE102484_1',
                'GSE102484_2',
                'GSE102484_3'
                )

MAIN_DIR = ''

for (dataset in DATASETS){ 
    gc()
    rm(list = apropos("cancer_affy"))

    # set working directory
    celpath = sprintf('%s/%s', MAIN_DIR, dataset)
    setwd(celpath)
    
    
    print(sprintf('untaring %s ...', dataset))
    untar(sprintf('%s_RAW.tar', dataset), exdir = "data")
    print(sprintf('%s untared', dataset))
     
    cels = list.files("data/", pattern = "CEL|cel")
     
    print(sprintf('unzipping %s ...', dataset))
    sapply(paste("data", cels, sep = "/"), gunzip)
    print(sprintf('%s unzipped', dataset))
    
    cels = list.files("data", pattern = "CEL|cel")
    
    # Changing the cells data to Affy to use in rma (ReadAffy)
    cancer_affy <- affy::ReadAffy(celfile.path = sprintf('%s/%s', celpath, 'data'))
    
    # Save the file
    save(cancer_affy, file=sprintf('%s/%s/%s_Affy.RData', MAIN_DIR, dataset, dataset))
    print(sprintf('%s affy saved', dataset))
    
    }
