# As there are many different tools and methods for getting alignments of reads, we will look at starting the process with two common input types. We will be using a
# count table, like that we would have if we were loading from a text file and we will use an ExpressionSet(eset) object, which is an object type common in
# Bioconductor. The prepared dataset will be the 'modencodefly' data from NHGRI encyclopedia of DNA elements project for the model organism, 'Drosophila melanogaster'.
# The dataset contains 147 different samples for 'D. melanogaster', a fruit fly with approximately 110 Mbp genome, annotated with about 15,000 gene features.

# edgeR ==> edgeR is a widely used and powerful package that implements negative binomial models suitable for sparse count data such as RNA-seq data in a general linear
# model framework, which are poweful for describing and comprehending count relationships and exact tests for multi-group experiments. It uses a weighted style
# normalization called "TMM", which is the weighted mean of log ratio between sample and control, after removal of genes with high counts and outlying log ratios.
# The TMM value should be close to one, but can be used as a corretion factor to be applied to overall library sizes. 
  
## Estimating differential expression with edgeR. For estimating differential expressions with edgeR from a count table (i.e. in a text file), we use the following
# steps:

## Dependencies.
library(forcats)
library(dplyr)
library(ggplot2)
library(magrittr)
library(edgeR)

## Count table
## 1. Load cont data.
count_dataframe <- readr::read_tsv(file.path(getwd(), "datasets", "ch1", "modencodefly_count_table.txt"))
genes <- count_dataframe[['gene']]
counts_dataframe[['gene']] <- NULL
count_matrix <- as.matrix(count_dataframe)
rownames(count_matrix) <- genes
pheno_data <- readr::read_table2(file.path(getwd(), "datasets", "ch1", "modencodefly_phenodata.txt"))

# read_table2--> read whitespace-seperated columns into a tibble; this function is deprecated because we renamed it to 'read_table()' and removed the old 'read_table' function, which was too strict for most cases and was analogous to just using 'read_fwf()'. 
  
## 2. Specify experiments of interest.
experiments_of_interest <- c("L1Larvae", "L2Larvae")
columns_of_interest <- which(pheno_data[['stage']]%in%experiments_of_interest)

## 3. For the grouping factor. 
grouping <- pheno_data[['stage']][columns_of_interest]%>%forcats::as_factor() # %>% is the 'infix operator' defined by magrittr and heavily used by dplyr. The function
# passes left-hand side of the operator to the first argument of the right-hand side of the operator. 

## 4. Form the subset of count data.
counts_of_interest <- count_matrix[,columns_of_interest]

## 5. Create the DGE object. 
count_dge <- edgeR::DGEList(counts = counts_of_interest, group = grouping)
# **DGEList() object holds the dataset to be analyzed by edgeR and the subsequent calculations performed on the dataset. 

## 6. Perform differential analysis.
design <- model.matrix(~grouping) # model.matrix() creates a design (or model), e.g., by expanding factors to a set of dummy variables (depending on the contrasts)
# and expanding interactions similarily.
eset_dge <- edgeR::estimateDisp(eset_dge, design) # estimateDisp() maximizes the negative binomial likelihood to give the estimate of the common, trended, and tagwise
# dispersions across all tags.
fit <- edgeR::glmQLFit(eset_dge, design)
result<- edgeR::glmQLFTet(fit, coef=2)
topTags(result)

