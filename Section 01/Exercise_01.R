# As there are many different tools and methods for getting alignments of reads, we will look at starting the process with two common input types. We will be using a 
# count table, like that we would have if we were loading from a text file and we will use an ExpressionSet(eset) object, which is an object type common in Bioconductor.
# The prepared dataset will be the 'modencodefly' data from NHGRI encyclopedia of DNA elements project for the model organism, _Drosophila melanogaster'. The dataset contains
# 147 different samples for 'D. melanogaster', a fruit fly with approximately 110 Mbp genome, annotated with about 15,000 gene features.

# edgeR ==> edgeR is a widely used and powerful package that implements negative binomial models suitable for sparse count data such as RNA-seq data in a general linear model framework, 
which are poweful for describing and comprehending count relationships and exact tests for multi-group experiments. It uses a weighted style normalization called "TMM", 
which is the weighted mean of log ratio between sample and control, after removal of genes with high counts and outlying log ratios. The TMM value should be close to 
one, but can be used as a corretion factor to be applied to overall library sizes. 
  
## Estimating differential expression with edgeR

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
columns_of_interest <- which( pheno_data[['stage']] %in% experiments_of_interest ) #gives the numeric indices of the above experiments in the matrix

### Form grouping factor
library(magrittr)
grouping <- pheno_data[['stage']][columns_of_interest] %>% 
  forcats::as_factor()

### Form subset count data
counts_of_interest <-  count_matrix[,columns_of_interest]

### Form DGE
library(edgeR)
count_dge <- edgeR::DGEList(counts = counts_of_interest, group = grouping)


### Perform differential expression 
design <- model.matrix(~ grouping)
eset_dge <- edgeR::estimateDisp(count_dge, design)
fit <- edgeR::glmQLFit(eset_dge, design)
result <- edgeR::glmQLFTest(fit, coef=2)
topTags(result)


##eset
load(file.path(getwd(), "datasets/ch1/modencodefly_eset.RData"))

experiments_of_interest <- c("L1Larvae", "L2Larvae")
columns_of_interest <- which( phenoData(modencodefly.eset)[['stage']] %in% experiments_of_interest )

grouping <- droplevels(phenoData(modencodefly.eset)[['stage']][columns_of_interest] )

counts_of_interest <- exprs(modencodefly.eset)[, columns_of_interest]

eset_dge <- edgeR::DGEList(
  counts = counts_of_interest,
  group = grouping 
  )

design <- model.matrix(~ grouping)
eset_dge <- edgeR::estimateDisp(eset_dge, design)
fit <- edgeR::glmQLFit(eset_dge, design)
result <- edgeR::glmQLFTest(fit, coef=2)
topTags(result)
