library(affy)      # For handling Affymetrix microarray data
library(sva)       # For the COMBAT method
library(limma)     # For various microarray data processing functions
MAIN_DIR = ''

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

# Changeing the cels data to Affy to use in rma (ReadAffy)
esets = list()
for(dataset in DATASETS) {
    env <- new.env()
    filename = sprintf('%s/%s/rma/%s.rma.RData', MAIN_DIR, dataset, dataset)
    print(sprintf('Loading file %s', filename))
    nm <- load(filename, env)[1]
    cancer_data <- env[[nm]]
    pData(cancer_data)$dataset <- dataset
    esets[[dataset]] <- cancer_data
}
# Get the common genes across all datasets
common_genes <- Reduce(intersect, lapply(esets, function(x) rownames(exprs(x))))

# Subset each dataset to only the common genes
esets_common <- lapply(esets, function(x) x[common_genes, ])

# Combine the expression matrices from each dataset
combined_expr <- do.call(cbind, lapply(esets_common, exprs))

# Combine the phenotype data (pData) from each dataset
combined_pheno <- do.call(rbind, lapply(esets_common, pData))

pheno_sample_names <- rownames(combined_pheno)
# pheno_sample_names <- sub("^.*\\.(GSM\\d+(_.*)?\\.CEL)$", "\\1", pheno_sample_names)
pheno_sample_names <- sub("^.*\\.(GSM\\d+(_.*)?\\.[cC][eE][lL]{1,2})$", "\\1", pheno_sample_names)

rownames(combined_pheno) <- pheno_sample_names


# Create a batch variable based on the dataset origin
batch <- factor(rep(names(esets_common), sapply(esets_common, ncol)))

# Create the combined ExpressionSet
combined_eset <- ExpressionSet(assayData=combined_expr, phenoData=AnnotatedDataFrame(combined_pheno))
save(combined_eset, file=sprintf('%s/analysed_datasets/merged/eset_merged_NON_COMBAT_rma.RData', MAIN_DIR))
save(batch, file=sprintf('%s/analysed_datasets/merged/eset_merged_batch.RData', MAIN_DIR))


# Apply COMBAT to remove batch effects
load(sprintf('%s/analysed_datasets/merged/eset_merged_NON_COMBAT_rma.RData', MAIN_DIR))
load(sprintf('%s/analysed_datasets/merged/eset_merged_batch.RData', MAIN_DIR))

combat_corrected <- ComBat(dat=exprs(combined_eset), batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
save(combat_corrected, file=sprintf('%s/analysed_datasets/merged/eset_merged_COMBAT_CORRECTED_rma.RData', MAIN_DIR))
load(sprintf('%s/analysed_datasets/merged/eset_merged_COMBAT_CORRECTED_rma.RData', MAIN_DIR))
# Update the ExpressionSet with the COMBAT corrected data
exprs(combined_eset) <- combat_corrected

# Save the final ExpressionSet
save(combined_eset, file=sprintf('%s/analysed_datasets/merged/eset_merged_COMBAT_rma.RData', MAIN_DIR))

# Multidimensional Scaling (MDS) plot
plotMDS(combined_eset, col=batch)

cancer_data.matrix <- affy::exprs(combined_eset)
# Assuming cancer_data.matrix is the expression matrix
sample_names <- colnames(cancer_data.matrix)

# Detect duplicate samples based on simplified naming (e.g., "GSM1234")
simplified_names <- sapply(strsplit(sample_names, "\\."), function(x) x[1])

# Identify duplicates
duplicates <- duplicated(simplified_names)
duplicated_samples <- sample_names[duplicates]

print("List of duplicated samples:")
print(duplicated_samples)

# Remove duplicates
cancer_data_uniq_matrix <- cancer_data.matrix[, !duplicates]

# Assign simplified unique names to the columns
colnames(cancer_data_uniq_matrix) <- simplified_names[!duplicates]

# Save the matrix to a CSV file
write.csv(cancer_data_uniq_matrix, file = sprintf('%s/analysed_datasets/merged/merged_COMBAT_rma.matrix.csv', MAIN_DIR), row.names = TRUE)



#######################################################################3



eset_COMBAT = inSilicoMerging::merge(esets, method="COMBAT")

# saveing affy file 
print(sprintf('Saving affy eset merged file'))
dir.create(sprintf('%s/analysed_datasets/merged', MAIN_DIR))
save(eset_COMBAT, file=sprintf('%s/analysed_datasets/merged/eset_merged_COMBAT_rma.RData', MAIN_DIR))

cancer_data.matrix <- affy::exprs(eset_COMBAT)
head(cancer_data.matrix)

# list of dublicated samples
print(sprintf('list of duplicated'))
len = ncol(cancer_data.matrix)
duplicated_names = list()
sample_names = colnames(cancer_data.matrix)
j = 0
for(i in 1:len){
    item = unlist(strsplit(unlist(strsplit(colnames(cancer_data.matrix)[i], split='.', fixed=TRUE))[1], split='_', fixed=TRUE))[1]
    
    if(!item %in% sample_names){
      sample_names[i] = item
    }
    else {
        duplicated_names[j] = sample_names[i]
        j = j + 1
   }
    
}

duplicated_names

# remove dublicate samples
print(sprintf('removing duplicated'))

len = ncol(cancer_data.matrix)
duplicated_names = list()
sample_names = colnames(cancer_data.matrix)
j = 0

for(i in 1:len){
    item = unlist(strsplit(unlist(strsplit(colnames(cancer_data.matrix)[i], split='.', fixed=TRUE))[1], split='_', fixed=TRUE))[1]
    sample_names[i] = item  
}

colnames(cancer_data.matrix)= sample_names

# save exptession as csv (.matrix.csv)
write.csv(cancer_data_uniq.matrix, file= sprintf('%s/analysed_datasets/merged/merged_COMBAT_rma.matrix.csv',MAIN_DIR), row.names=TRUE)

plotMDS(eset_COMBAT)
