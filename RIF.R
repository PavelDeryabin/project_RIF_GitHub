
########################################################################################################

#                                    RIF endometrium mRNA analysis

########################################################################################################

set.seed(42)

library(GEOquery)
library(psych)
library(sva)
library(preprocessCore)
library(endest)
library(limma)
library(splines)
library(org.Hs.eg.db)
library(tximeta)
library(DESeq2)
library(biomaRt)
library(gprofiler2)
library(tximport)
library(Seurat)
library(SeuratWrappers)
library(DoubletFinder)
library(harmony)
library(UCell)
library(clusterProfiler)
library(msigdbr)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggExtra)
library(patchwork)
library(RColorBrewer)
library(ggVennDiagram)
library(lemon)
library(pals) 

setwd('/home/pd/projects/article_2025_RIF/project_RIF_GitHub/')

writeLines(capture.output(sessionInfo()), 'sessionInfo.txt')




### Prepare metadata and count matrices of GEO datasets

## GSE103465

# Get the data
dataset.GSE103465 <- getGEO('GSE103465')
dataset.GSE103465
dataset.GSE103465 <- dataset.GSE103465[[1]]

# Extract and set mdata
mdata.GSE103465 <- pData(dataset.GSE103465)
View(mdata.GSE103465) # Look for possible batch variables
table(mdata.GSE103465$contact_country)
t.test(as.numeric(`age:ch1`) ~ `subject status:ch1`, mdata.GSE103465)
table(mdata.GSE103465$`sampling time:ch1`)

mdata.GSE103465$names <- rownames(mdata.GSE103465)
mdata.GSE103465$GEO.accession <- 'GSE103465'
mdata.GSE103465$GSM.accession <- mdata.GSE103465$geo_accession
mdata.GSE103465$Sample.name <- mdata.GSE103465$title
mdata.GSE103465$Phenotype <- mdata.GSE103465$`subject status:ch1`
mdata.GSE103465$Menstrual.cycle.time <- mdata.GSE103465$`sampling time:ch1`
mdata.GSE103465$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE103465$Age.years <- mdata.GSE103465$`age:ch1`
mdata.GSE103465 <- mdata.GSE103465[, c(40:47)]

table(mdata.GSE103465$Phenotype)
mdata.GSE103465 <- mdata.GSE103465 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('health control') ~ 'Fertile',
    Phenotype %in% c('recurrent implantation failure (RIF) patient') ~ 'RIF'))

# Extract exp mtx
exp.GSE103465 <- exprs(dataset.GSE103465)
class(exp.GSE103465)
dim(exp.GSE103465)
exp.GSE103465[1:5, 1:5]
boxplot(exp.GSE103465) # as indicated at GEO, the dataset was background corrected and log-transformed using RMA
exp.GSE103465 <- normalize.quantiles(exp.GSE103465, keep.names = TRUE) # quantile-normalize

# Extract features data
fdata.GSE103465 <- fData(dataset.GSE103465)
head(fdata.GSE103465)

# Annotate exp mtx by genes symbols, select probes with max mean value across the samples and set rownames to symbols
exp.GSE103465 <- as.data.frame(exp.GSE103465)
exp.GSE103465[1:5, ]
exp.GSE103465$GENE_SYMBOL <- as.character(fdata.GSE103465$`Gene Symbol`)
exp.GSE103465$rowMean <- rowMeans(exp.GSE103465[, 1:6])
head(exp.GSE103465)

exp.GSE103465.max <- exp.GSE103465 %>% arrange(GENE_SYMBOL, -rowMean) %>% group_by(GENE_SYMBOL) %>% slice_head(n = 1) %>% as.data.frame()
head(exp.GSE103465.max)
exp.GSE103465[exp.GSE103465$GENE_SYMBOL %in% 'HIST1H3G', ]
exp.GSE103465.max[exp.GSE103465.max$GENE_SYMBOL %in% 'HIST1H3G', ]

sum(duplicated(fdata.GSE103465$`Gene Symbol`))
sum(duplicated(exp.GSE103465$GENE_SYMBOL))
sum(duplicated(exp.GSE103465.max$GENE_SYMBOL))
rownames(exp.GSE103465.max) <- exp.GSE103465.max$GENE_SYMBOL
rownames(exp.GSE103465.max)[1:50]
exp.GSE103465.max <- exp.GSE103465.max[-c(1:28), c(1:6)] # remove uninformative features and technical columns
head(exp.GSE103465.max) # clean exp mtx
boxplot(exp.GSE103465.max) # normalized and log transformed

# EDA
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE103465.max))), 500))
pc <- prcomp(t(exp.GSE103465.max[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE103465$PC1 <- pca$PC1
mdata.GSE103465$PC2 <- pca$PC2
ggplot(data = mdata.GSE103465, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE103465[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE103465 <- estimate_cycle_time(exprs = as.matrix(exp.GSE103465.max),
                                                entrez_ids = mapIds(org.Hs.eg.db,
                                                                    keys = rownames(exp.GSE103465.max),
                                                                    column = 'ENTREZID',
                                                                    keytype = 'SYMBOL',
                                                                    multiVals = first))
mdata.GSE103465$EndEst <- endest.results.GSE103465$estimated_time

ggplot(data = mdata.GSE103465, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE103465[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype') + scale_color_viridis_c()

# Save results
write.csv(exp.GSE103465.max, 'processed.data/exp.GSE103465.max.csv')
write.csv(mdata.GSE103465, file = 'processed.data/mdata.GSE103465.csv')



## GSE111974

# Get the data
dataset.GSE111974 <- getGEO('GSE111974')
dataset.GSE111974
dataset.GSE111974 <- dataset.GSE111974[[1]]

# Extract and set mdata
mdata.GSE111974 <- pData(dataset.GSE111974)
View(mdata.GSE111974) # Look for possible batch variables
mdata.GSE111974$names <- rownames(mdata.GSE111974)
mdata.GSE111974$GEO.accession <- 'GSE111974'
mdata.GSE111974$GSM.accession <- mdata.GSE111974$geo_accession
mdata.GSE111974$Sample.name <- mdata.GSE111974$title
mdata.GSE111974$Phenotype <- sub('\\_.*', '', 
                                 sub('Endometrial_Tissue*.', '', 
                                     sub('Endometrial_Tissue_Fertile*.', '', mdata.GSE111974$title)))
mdata.GSE111974$Menstrual.cycle.time <- NA
mdata.GSE111974$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE111974$Age.years <- NA
mdata.GSE111974 <- mdata.GSE111974[, c(32:39)]

table(mdata.GSE111974$Phenotype)
mdata.GSE111974 <- mdata.GSE111974 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('Control') ~ 'Fertile',
    Phenotype %in% c('RIF') ~ 'RIF'))

# Extract exp mtx
exp.GSE111974 <- exprs(dataset.GSE111974)
class(exp.GSE111974)
dim(exp.GSE111974)
exp.GSE111974[1:5, 1:5]
boxplot(exp.GSE111974) # as indicated in GEO, exp values are QN normalized, and RMA corrected and log transformed

# Extract features data
fdata.GSE111974 <- fData(dataset.GSE111974)
head(fdata.GSE111974)

# Annotate exp mtx by genes symbols, select probes with max mean value across the samples and set rownames to symbols
# This approach was applied in the present study to process any microarray data
exp.GSE111974 <- as.data.frame(exp.GSE111974)
exp.GSE111974$GENE_SYMBOL <- as.character(fdata.GSE111974$GENE_SYMBOL)
exp.GSE111974$rowMean <- rowMeans(exp.GSE111974[, 1:48])
head(exp.GSE111974)
exp.GSE111974.max <- exp.GSE111974 %>% arrange(GENE_SYMBOL, -rowMean) %>% group_by(GENE_SYMBOL) %>% slice_head(n = 1) %>% as.data.frame()
head(exp.GSE111974.max)
rownames(exp.GSE111974.max) <- exp.GSE111974.max$GENE_SYMBOL
rownames(exp.GSE111974.max)[1:50]
exp.GSE111974.max <- exp.GSE111974.max[-c(1:2), c(1:48)] # remove uninformative features and technical columns
head(exp.GSE111974.max) # clean exp mtx
boxplot(exp.GSE111974.max) 

# Annotate exp mtx by genes symbols, average probes signals for a gene () and set rownames to symbols 
# This approach was implemented only ones to replicate analysis of the original study
exp.GSE111974.avrg <- exp.GSE111974[, -50] %>% group_by(GENE_SYMBOL) %>% 
  summarise(across(starts_with('GSM'), mean, na.rm = TRUE)) %>% as.data.frame()
rownames(exp.GSE111974.avrg) <- exp.GSE111974.avrg$GENE_SYMBOL
rownames(exp.GSE111974.avrg)[1:50]
exp.GSE111974.avrg <- exp.GSE111974.avrg[-c(1:2), c(2:49)] # remove uninformative features and technical columns
head(exp.GSE111974.avrg) # clean exp mtx
boxplot(exp.GSE111974.avrg)

# EDA
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE111974.max))), 500))
pc <- prcomp(t(exp.GSE111974.max[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE111974$PC1 <- pca$PC1
mdata.GSE111974$PC2 <- pca$PC2

ggplot(data = mdata.GSE111974, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE111974[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE111974 <- estimate_cycle_time(exprs = as.matrix(exp.GSE111974.max),
                                                entrez_ids = mapIds(org.Hs.eg.db,
                                                                    keys = rownames(exp.GSE111974.max),
                                                                    column = 'ENTREZID',
                                                                    keytype = 'SYMBOL',
                                                                    multiVals = first))
mdata.GSE111974$EndEst <- endest.results.GSE111974$estimated_time

ggplot(data = mdata.GSE111974, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE111974[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype') + scale_color_viridis_c()

# Save results
write.csv(exp.GSE111974.max, 'processed.data/exp.GSE111974.max.csv')
write.csv(exp.GSE111974.avrg, 'processed.data/exp.GSE111974.avrg.csv')
write.csv(mdata.GSE111974, file = 'processed.data/mdata.GSE111974.csv')



## GSE188409

# Get the data
dataset.GSE188409 <- getGEO('GSE188409')
dataset.GSE188409
dataset.GSE188409 <- dataset.GSE188409[[1]]

# Extract and set mdata
mdata.GSE188409 <- pData(dataset.GSE188409)
View(mdata.GSE188409) # Look for possible batch variables
mdata.GSE188409$names <- rownames(mdata.GSE188409)
mdata.GSE188409$GEO.accession <- 'GSE188409'
mdata.GSE188409$GSM.accession <- mdata.GSE188409$geo_accession
mdata.GSE188409$Sample.name <- mdata.GSE188409$title
mdata.GSE188409$Phenotype <- mdata.GSE188409$`disease state:ch1`
mdata.GSE188409$Menstrual.cycle.time <- NA
mdata.GSE188409$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE188409$Age.years <- NA
mdata.GSE188409 <- mdata.GSE188409[, c(35:42)]

table(mdata.GSE188409$Phenotype)
mdata.GSE188409 <- mdata.GSE188409 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('normal') ~ 'Fertile',
    Phenotype %in% c('RIF') ~ 'RIF'))

# Extract exp mtx
exp.GSE188409 <- exprs(dataset.GSE188409)
class(exp.GSE188409)
dim(exp.GSE188409)
exp.GSE188409[1:5, 1:5]
boxplot(exp.GSE188409) # as indicated in GEO, Aligent global normalization was applied - using this values directly

# Extract features data
fdata.GSE188409 <- fData(dataset.GSE188409)
head(fdata.GSE188409)

# Annotate exp mtx by genes symbols, select probes with max mean value across the samples and set rownames to symbols
exp.GSE188409 <- as.data.frame(exp.GSE188409)
exp.GSE188409[1:5, ]
exp.GSE188409$GENE_SYMBOL <- as.character(fdata.GSE188409$ORF)
exp.GSE188409$rowMean <- rowMeans(exp.GSE188409[, 1:10])
head(exp.GSE188409)

exp.GSE188409.max <- exp.GSE188409 %>% arrange(GENE_SYMBOL, -rowMean) %>% group_by(GENE_SYMBOL) %>% slice_head(n = 1) %>% as.data.frame()
head(exp.GSE188409.max)

rownames(exp.GSE188409.max) <- exp.GSE188409.max$GENE_SYMBOL
rownames(exp.GSE188409.max)[1:50]
exp.GSE188409.max <- exp.GSE188409.max[-c(1), c(1:10)] # remove uninformative features and technical columns
head(exp.GSE188409.max) # clean exp mtx
boxplot(exp.GSE188409.max) # normalized and log transformed

# EDA
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE188409.max))), 500))
pc <- prcomp(t(exp.GSE188409.max[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE188409$PC1 <- pca$PC1
mdata.GSE188409$PC2 <- pca$PC2

ggplot(data = mdata.GSE188409, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE188409[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE188409 <- estimate_cycle_time(exprs = as.matrix(exp.GSE188409.max),
                                                entrez_ids = mapIds(org.Hs.eg.db,
                                                                    keys = rownames(exp.GSE188409.max),
                                                                    column = 'ENTREZID',
                                                                    keytype = 'SYMBOL',
                                                                    multiVals = first))
mdata.GSE188409$EndEst <- endest.results.GSE188409$estimated_time

ggplot(data = mdata.GSE188409, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE188409[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype') + scale_color_viridis_c()

# Save results
write.csv(exp.GSE188409.max, 'processed.data/exp.GSE188409.max.csv')
write.csv(mdata.GSE188409, file = 'processed.data/mdata.GSE188409.csv')



## GSE26787

# Get the data
dataset.GSE26787 <- getGEO('GSE26787')
dataset.GSE26787
dataset.GSE26787 <- dataset.GSE26787[[1]]

# Extract and set mdata
mdata.GSE26787 <- pData(dataset.GSE26787)
View(mdata.GSE26787) # Look for possible batch variables
mdata.GSE26787$names <- rownames(mdata.GSE26787)
mdata.GSE26787$GEO.accession <- 'GSE26787'
mdata.GSE26787$GSM.accession <- mdata.GSE26787$geo_accession
mdata.GSE26787$Sample.name <- mdata.GSE26787$title
mdata.GSE26787$Phenotype <- mdata.GSE26787$source_name_ch1
mdata.GSE26787$Menstrual.cycle.time <- NA
mdata.GSE26787$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE26787$Age.years <- NA
mdata.GSE26787 <- mdata.GSE26787[, c(37:44)]

table(mdata.GSE26787$Phenotype)
mdata.GSE26787 <- mdata.GSE26787 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('fertile') ~ 'Fertile',
    Phenotype %in% c('implantation failure') ~ 'RIF',
    Phenotype %in% c('recurrent spontaneous abortion') ~ 'RPL'))

# Extract exp mtx
exp.GSE26787 <- exprs(dataset.GSE26787)
class(exp.GSE26787)
dim(exp.GSE26787)
exp.GSE26787[1:5, 1:5]
boxplot(exp.GSE26787) # as indicated in GEO, GC-RMA normalization was applied - using this values directly

# Extract features data
fdata.GSE26787 <- fData(dataset.GSE26787)
head(fdata.GSE26787)

# Annotate exp mtx by genes symbols, select probes with max mean value across the samples and set rownames to symbols
exp.GSE26787 <- as.data.frame(exp.GSE26787)
exp.GSE26787[1:5, ]
exp.GSE26787$GENE_SYMBOL <- gsub(' ///.*', '', fdata.GSE26787$`Gene Symbol`)
exp.GSE26787$rowMean <- rowMeans(exp.GSE26787[, 1:15])
head(exp.GSE26787)

exp.GSE26787.max <- exp.GSE26787 %>% arrange(GENE_SYMBOL, -rowMean) %>% group_by(GENE_SYMBOL) %>% slice_head(n = 1) %>% as.data.frame()
head(exp.GSE26787.max)

sum(duplicated(fdata.GSE26787$`Gene Symbol`))
sum(duplicated(exp.GSE26787$GENE_SYMBOL))
sum(duplicated(exp.GSE26787.max$GENE_SYMBOL))
rownames(exp.GSE26787.max) <- exp.GSE26787.max$GENE_SYMBOL
rownames(exp.GSE26787.max)[1:50]
exp.GSE26787.max <- exp.GSE26787.max[-c(1), c(1:15)] # remove uninformative features and technical columns
head(exp.GSE26787.max) # clean exp mtx
boxplot(exp.GSE26787.max) # normalized and log transformed

# EDA
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE26787.max))), 500))
pc <- prcomp(t(exp.GSE26787.max[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE26787$PC1 <- pca$PC1
mdata.GSE26787$PC2 <- pca$PC2
ggplot(data = mdata.GSE26787, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE26787[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE26787 <- estimate_cycle_time(exprs = as.matrix(exp.GSE26787.max),
                                               entrez_ids = mapIds(org.Hs.eg.db,
                                                                   keys = rownames(exp.GSE26787.max),
                                                                   column = 'ENTREZID',
                                                                   keytype = 'SYMBOL',
                                                                   multiVals = first))
mdata.GSE26787$EndEst <- endest.results.GSE26787$estimated_time

ggplot(data = mdata.GSE26787, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE26787[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype') + scale_color_viridis_c()

# Save results
write.csv(subset(mdata.GSE26787, Phenotype != 'RPL'), file = 'processed.data/mdata.GSE26787.csv')
write.csv(exp.GSE26787.max[, rownames(subset(mdata.GSE26787, Phenotype != 'RPL'))], 'processed.data/exp.GSE26787.max.csv')



## GSE58144

# Get the data
dataset.GSE58144 <- getGEO('GSE58144')
dataset.GSE58144
dataset.GSE58144 <- dataset.GSE58144[[1]]

# Extract and set mdata
mdata.GSE58144 <- pData(dataset.GSE58144)
View(mdata.GSE58144) # Look for possible batch variables
mdata.GSE58144$names <- rownames(mdata.GSE58144)
mdata.GSE58144$GEO.accession <- 'GSE58144'
mdata.GSE58144$GSM.accession <- mdata.GSE58144$geo_accession
mdata.GSE58144$Sample.name <- mdata.GSE58144$title
mdata.GSE58144$Phenotype <- mdata.GSE58144$characteristics_ch1.1
mdata.GSE58144$Menstrual.cycle.time <- NA
mdata.GSE58144$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE58144$Age.years <- mdata.GSE58144$characteristics_ch1.6
mdata.GSE58144$Batch <- mdata.GSE58144$`batch:ch1`
mdata.GSE58144 <- mdata.GSE58144[, c(63:71)]

table(mdata.GSE58144$Phenotype)
mdata.GSE58144 <- mdata.GSE58144 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('condition: control') ~ 'Fertile',
    Phenotype %in% c('condition: recurrent implantation failure after IVF') ~ 'RIF'))

mdata.GSE58144$Age.years <- sub('age: *', '', mdata.GSE58144$Age.years)

# Extract exp mtx
exp.GSE58144 <- exprs(dataset.GSE58144)
class(exp.GSE58144)
dim(exp.GSE58144)
exp.GSE58144[1:5, 1:5]
boxplot(exp.GSE58144[, 1:10]) # as.indicated in GEO, these are log2 ratio (sample/reference) values, download soft.gz file

if (! file.exists('processed.data/GSE58144_family.soft.gz')) {
  archive_url <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58144/soft/GSE58144_family.soft.gz'
  download.file(archive_url, destfile = 'processed.data/GSE58144_family.soft.gz', method = 'curl')
}
dataset.GSE58144.soft_file <- getGEO('GSE58144', filename = 'processed.data/GSE58144_family.soft.gz') 
length(dataset.GSE58144.soft_file@gsms[['GSM1402321']]@dataTable@table[['ID_REF']])
sum(rownames(exp.GSE58144) == dataset.GSE58144.soft_file@gsms[['GSM1402321']]@dataTable@table[['ID_REF']])
sum(rownames(exp.GSE58144) %in% dataset.GSE58144.soft_file@gsms[['GSM1402321']]@dataTable@table[['ID_REF']])
exp.GSE58144 <- exp.GSE58144[as.character(dataset.GSE58144.soft_file@gsms[['GSM1402321']]@dataTable@table[['ID_REF']]), ]
exp.GSE58144[1:10, 'GSM1402321']
dataset.GSE58144.soft_file@gsms[['GSM1402321']]@dataTable@table[['VALUE']][1:10] # same log2 ratio (sample/reference) values
for (sample in colnames(exp.GSE58144)) {
  exp.GSE58144[, sample] <- dataset.GSE58144.soft_file@gsms[[sample]]@dataTable@table[['A_VALUE']]
}
boxplot(exp.GSE58144[, 1:10]) # Quantile normalized and corrected 1/2log2(sample*reference) - using these values

# Extract features data
fdata.GSE58144 <- fData(dataset.GSE58144)
head(fdata.GSE58144)
fdata.GSE58144 <- fdata.GSE58144[rownames(exp.GSE58144), ]

# Annotate exp mtx by genes symbols, select probes with max mean value across the samples and set rownames to symbols
exp.GSE58144 <- as.data.frame(exp.GSE58144)
exp.GSE58144[1:5, ]
exp.GSE58144$GENE_SYMBOL <- as.character(fdata.GSE58144$GENE_SYMBOL)
exp.GSE58144$rowMean <- rowMeans(exp.GSE58144[, 1:115])
head(exp.GSE58144)

exp.GSE58144.max <- exp.GSE58144 %>% arrange(GENE_SYMBOL, -rowMean) %>% group_by(GENE_SYMBOL) %>% slice_head(n = 1) %>% as.data.frame()
head(exp.GSE58144.max)

rownames(exp.GSE58144.max) <- exp.GSE58144.max$GENE_SYMBOL
rownames(exp.GSE58144.max)[1:50]
exp.GSE58144.max <- exp.GSE58144.max[-c(1), c(1:115)] # remove uninformative features and technical columns
head(exp.GSE58144.max) # clean exp mtx
boxplot(exp.GSE58144.max[, 1:10])

# EDA
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE58144.max))), 500))
pc <- prcomp(t(exp.GSE58144.max[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE58144$PC1 <- pca$PC1
mdata.GSE58144$PC2 <- pca$PC2
ggplot(data = mdata.GSE58144, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE58144[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

ggplot(data = mdata.GSE58144, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE58144[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Batch') # significant batch effect - considering it in DEA model further

# Remove batch 2 cohort 2 samples and check whether the Phenotype variable is balanced
table(mdata.GSE58144$Phenotype, mdata.GSE58144$Batch)
mdata.GSE58144.subset <- mdata.GSE58144 %>% filter(Batch != 'cohort 2, batch 2')
exp.GSE58144.max.subset <- exp.GSE58144.max[, rownames(mdata.GSE58144.subset)]

# Estimate relative menstrual cycle progression
endest.results.GSE58144.subset <- estimate_cycle_time(exprs = as.matrix(exp.GSE58144.max.subset),
                                                      entrez_ids = mapIds(org.Hs.eg.db,
                                                                    keys = rownames(exp.GSE58144.max.subset),
                                                                    column = 'ENTREZID',
                                                                    keytype = 'SYMBOL',
                                                                    multiVals = first))
mdata.GSE58144.subset$EndEst <- endest.results.GSE58144.subset$estimated_time

hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE58144.max.subset))), 500))
pc <- prcomp(t(exp.GSE58144.max.subset[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE58144.subset$PC1 <- pca$PC1
mdata.GSE58144.subset$PC2 <- pca$PC2

ggplot(data = mdata.GSE58144.subset, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE58144.subset[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Batch') # significant batch effect - considering it in DEA model further

ggplot(data = mdata.GSE58144.subset, aes(x = PC1, y = PC2, color = EndEst)) + # now main variance is here
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE58144[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
table(mdata.GSE58144.subset$Batch)
mdata.GSE58144.subset <- mdata.GSE58144.subset %>%
  mutate(Batch = case_when(
    Batch == 'cohort 1' ~ 'Cohort 1',
    Batch == 'cohort 2, batch 1' ~ 'Cohort 2, batch 1'))

write.csv(mdata.GSE58144.subset, file = 'processed.data/mdata.GSE58144.subset.csv')
write.csv(exp.GSE58144.max.subset, 'processed.data/exp.GSE58144.max.subset.csv')



## GSE71331

# Get the data
dataset.GSE71331 <- getGEO('GSE71331')
dataset.GSE71331
dataset.GSE71331 <- dataset.GSE71331[[1]]

# Extract and set mdata
mdata.GSE71331 <- pData(dataset.GSE71331)
View(mdata.GSE71331) # Look for possible batch variables
mdata.GSE71331$names <- rownames(mdata.GSE71331)
mdata.GSE71331$GEO.accession <- 'GSE71331'
mdata.GSE71331$GSM.accession <- mdata.GSE71331$geo_accession
mdata.GSE71331$Sample.name <- mdata.GSE71331$title
mdata.GSE71331$Phenotype <- mdata.GSE71331$characteristics_ch1
mdata.GSE71331$Menstrual.cycle.time <- NA
mdata.GSE71331$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE71331$Age.years <- NA
mdata.GSE71331 <- mdata.GSE71331[, c(39:46)]

table(mdata.GSE71331$Phenotype)
mdata.GSE71331 <- mdata.GSE71331 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('patients: Control') ~ 'Fertile',
    Phenotype %in% c('patients: RIF (recurrent implantation failure)') ~ 'RIF'))

# Extract exp mtx
exp.GSE71331 <- exprs(dataset.GSE71331)
class(exp.GSE71331)
dim(exp.GSE71331)
exp.GSE71331[1:5, 1:5]
boxplot(exp.GSE71331[, 1:10]) # as indicated at GEO, dataset was corrected, transformed and normalized 
# by percentile normalization using Agilent Technologies software
exp.GSE71331 <- exp.GSE71331 - min(exp.GSE71331) # making the mtx positive and using it for analysis

# Extract features data
fdata.GSE71331 <- fData(dataset.GSE71331)
# fdata df does not contain probes annotation, so retriving it manually from biomaRt as decribed here
# https://support.bioconductor.org/p/9135846/#9135885
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
tables <- listAttributes(ensembl)
tables[grep('agilent', tables[, 1]), ] # looking for GPL19072	Agilent-052909 CBC_lncRNAmRNA_V3
annot <- getBM(attributes = c('agilent_gpl19072', 
                 'wikigene_description',
                 'ensembl_gene_id',
                 'entrezgene_id',
                 'gene_biotype',
                 'external_gene_name'),
  filters = 'agilent_gpl19072',
  values = fdata.GSE71331$ID,
  mart = ensembl)
annot.to.fdata <- merge(x = fdata.GSE71331, y = annot, by.x = 'ID', by.y = 'agilent_gpl19072', all.x = T)
annot.to.fdata <- annot.to.fdata[!(duplicated(annot.to.fdata$ID)), ]
rownames(annot.to.fdata) <- annot.to.fdata$ID
annot.to.fdata <- annot.to.fdata[rownames(fdata.GSE71331), ]
head(annot.to.fdata, 50)

# Annotate exp mtx by genes symbols, select probes with max mean value across the samples and set rownames to symbols
exp.GSE71331 <- as.data.frame(exp.GSE71331)
exp.GSE71331[1:5, ]
exp.GSE71331$GENE_SYMBOL <- as.character(annot.to.fdata$external_gene_name)
exp.GSE71331$rowMean <- rowMeans(exp.GSE71331[, 1:12])
head(exp.GSE71331)

exp.GSE71331.max <- exp.GSE71331 %>% arrange(GENE_SYMBOL, -rowMean) %>% group_by(GENE_SYMBOL) %>% slice_head(n = 1) %>% as.data.frame()
head(exp.GSE71331.max)

exp.GSE71331.max <- exp.GSE71331.max[!(is.na(exp.GSE71331.max$GENE_SYMBOL)), ]
rownames(exp.GSE71331.max) <- exp.GSE71331.max$GENE_SYMBOL
rownames(exp.GSE71331.max)[1:50]
exp.GSE71331.max <- exp.GSE71331.max[-c(1), c(1:12)] # remove uninformative features and technical columns
head(exp.GSE71331.max) # clean exp mtx
boxplot(exp.GSE71331.max) # normalized and log transformed

# EDA
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE71331.max))), 500))
pc <- prcomp(t(exp.GSE71331.max[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE71331$PC1 <- pca$PC1
mdata.GSE71331$PC2 <- pca$PC2
ggplot(data = mdata.GSE71331, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE71331[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE71331 <- estimate_cycle_time(exprs = as.matrix(exp.GSE71331.max),
                                               entrez_ids = mapIds(org.Hs.eg.db,
                                                                   keys = rownames(exp.GSE71331.max),
                                                                   column = 'ENTREZID',
                                                                   keytype = 'SYMBOL',
                                                                   multiVals = first))
mdata.GSE71331$EndEst <- endest.results.GSE71331$estimated_time

ggplot(data = mdata.GSE71331, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE71331[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
write.csv(exp.GSE71331.max, 'processed.data/exp.GSE71331.max.csv')
write.csv(mdata.GSE71331, file = 'processed.data/mdata.GSE71331.csv')



## GSE92324

# Get the data
dataset.GSE92324 <- getGEO('GSE92324')
dataset.GSE92324
dataset.GSE92324 <- dataset.GSE92324[[1]]

# Extract and set mdata
mdata.GSE92324 <- pData(dataset.GSE92324)
View(mdata.GSE92324) # Look for possible batch variables
mdata.GSE92324$names <- rownames(mdata.GSE92324)
mdata.GSE92324$GEO.accession <- 'GSE92324'
mdata.GSE92324$GSM.accession <- mdata.GSE92324$geo_accession
mdata.GSE92324$Sample.name <- mdata.GSE92324$title
mdata.GSE92324$Phenotype <- sub('\\_.*', '', mdata.GSE92324$`sample group and identifier:ch1`)
mdata.GSE92324$Menstrual.cycle.time <- mdata.GSE92324$`sample collection day:ch1`
mdata.GSE92324$Menstrual.cycle.phase <- 'Stimulated WOI'
mdata.GSE92324$Age.years <- mdata.GSE92324$`age:ch1`
mdata.GSE92324 <- mdata.GSE92324[, c(39:46)]

table(mdata.GSE92324$Phenotype)
mdata.GSE92324 <- mdata.GSE92324 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('control') ~ 'Fertile',
    Phenotype %in% c('case') ~ 'RIF'))

# Extract exp mtx
exp.GSE92324 <- exprs(dataset.GSE92324)
class(exp.GSE92324)
dim(exp.GSE92324)
exp.GSE92324[1:5, 1:5]
boxplot(exp.GSE92324) # GEO: The data was analysed using Genome studio software 2011. QN was used to normalize the data
exp.GSE92324 <- log2(exp.GSE92324 + abs(min(exp.GSE92324)) + 1) # adjusting values to become positive and log2 transforming
boxplot(exp.GSE92324) 

# Extract features data
fdata.GSE92324 <- fData(dataset.GSE92324)
head(fdata.GSE92324)

# Annotate exp mtx by genes symbols, select probes with max mean value across the samples and set rownames to symbols
exp.GSE92324 <- as.data.frame(exp.GSE92324)
exp.GSE92324[1:5, ]
exp.GSE92324$GENE_SYMBOL <- as.character(fdata.GSE92324$ILMN_Gene)
exp.GSE92324$rowMean <- rowMeans(exp.GSE92324[, 1:18])
head(exp.GSE92324)

exp.GSE92324.max <- exp.GSE92324 %>% arrange(GENE_SYMBOL, -rowMean) %>% group_by(GENE_SYMBOL) %>% slice_head(n = 1) %>% as.data.frame()
head(exp.GSE92324.max)

sum(duplicated(fdata.GSE92324$ILMN_Gene))
sum(duplicated(exp.GSE92324$GENE_SYMBOL))
sum(duplicated(exp.GSE92324.max$GENE_SYMBOL))
rownames(exp.GSE92324.max) <- exp.GSE92324.max$GENE_SYMBOL
rownames(exp.GSE92324.max)[1:50]
exp.GSE92324.max <- exp.GSE92324.max[-c(1:13), c(1:18)] # remove uninformative features and technical columns
head(exp.GSE92324.max) # clean exp mtx
boxplot(exp.GSE92324.max)

# EDA
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE92324.max))), 500))
pc <- prcomp(t(exp.GSE92324.max[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE92324$PC1 <- pca$PC1
mdata.GSE92324$PC2 <- pca$PC2
ggplot(data = mdata.GSE92324, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE92324[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE92324 <- estimate_cycle_time(exprs = as.matrix(exp.GSE92324.max),
                                               entrez_ids = mapIds(org.Hs.eg.db,
                                                                   keys = rownames(exp.GSE92324.max),
                                                                   column = 'ENTREZID',
                                                                   keytype = 'SYMBOL',
                                                                   multiVals = first))
mdata.GSE92324$EndEst <- endest.results.GSE92324$estimated_time

ggplot(data = mdata.GSE92324, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE92324[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
write.csv(mdata.GSE92324, file = 'processed.data/mdata.GSE92324.csv')
write.csv(exp.GSE92324.max, 'processed.data/exp.GSE92324.max.csv')



### PRJNA314429 

# Prepare metadata
mdata.PRJNA314429 <- read.csv('/home/pd/datasets/hs_endometrium/bulk_rna_seq/PRJNA314429/SraRunTable.txt', header = T, sep = ',', row.names = 1)
head(mdata.PRJNA314429)

mdata.PRJNA314429$names <- rownames(mdata.PRJNA314429)
mdata.PRJNA314429$GEO.accession <- 'PRJNA314429'
mdata.PRJNA314429$GSM.accession <- mdata.PRJNA314429$BioSample
mdata.PRJNA314429$Sample.name <- mdata.PRJNA314429$Sample.Name
mdata.PRJNA314429$Phenotype <- gsub('[0-9]+', '', mdata.PRJNA314429$Sample.Name)
mdata.PRJNA314429$Menstrual.cycle.time <- NA
mdata.PRJNA314429$Menstrual.cycle.phase <- 'Secretory mid'
mdata.PRJNA314429$Age.years <- mdata.PRJNA314429$Age
mdata.PRJNA314429$files <- file.path('/home/pd/datasets/hs_endometrium/bulk_rna_seq/PRJNA314429/quants/', mdata.PRJNA314429$names, 'quant.sf')
file.exists(mdata.PRJNA314429$files)
mdata.PRJNA314429 <- mdata.PRJNA314429[, c(31:39)]

table(mdata.PRJNA314429$Phenotype)
mdata.PRJNA314429 <- subset(mdata.PRJNA314429, Phenotype != 'RM') # keep only Fertile and RIF samples
mdata.PRJNA314429 <- mdata.PRJNA314429 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('C') ~ 'Fertile',
    Phenotype %in% c('RIF') ~ 'RIF'))

# Read counts and create DEseq2 object
se.PRJNA314429 <- tximeta(mdata.PRJNA314429)
gse.PRJNA314429 <- summarizeToGene(se.PRJNA314429)
gse.PRJNA314429 <- addIds(gse.PRJNA314429, 'SYMBOL', gene = TRUE)
ddsTC.PRJNA314429 <- DESeqDataSet(gse.PRJNA314429, design = ~ 1) # no adjustment of variation

# Keep features having maximum expression per known gene symbol identifier and at least 1 count in all samples
ddsTC.PRJNA314429 <- ddsTC.PRJNA314429[order(rowMeans(assay(ddsTC.PRJNA314429)), decreasing = TRUE), ]
ddsTC.PRJNA314429 <- ddsTC.PRJNA314429[
  !duplicated(rowData(ddsTC.PRJNA314429)$SYMBOL) & !(is.na(rowData(ddsTC.PRJNA314429)$SYMBOL)), ]
rownames(ddsTC.PRJNA314429) <- rowData(ddsTC.PRJNA314429)$SYMBOL
ddsTC.PRJNA314429 <- ddsTC.PRJNA314429[rowSums(counts(ddsTC.PRJNA314429) >= 1) >= ncol(counts(ddsTC.PRJNA314429)), ]

# EDA
vst.PRJNA314429.mtx <- assay(DESeq2::vst(ddsTC.PRJNA314429, blind = TRUE))
boxplot(vst.PRJNA314429.mtx)

hvg <- names(tail(sort(rowVars(vst.PRJNA314429.mtx)), 500))
pc <- prcomp(t(vst.PRJNA314429.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.PRJNA314429$PC1 <- pca$PC1
mdata.PRJNA314429$PC2 <- pca$PC2

ggplot(data = mdata.PRJNA314429, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.PRJNA314429[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.PRJNA314429 <- estimate_cycle_time(exprs = vst.PRJNA314429.mtx,
                                                  entrez_ids = mapIds(org.Hs.eg.db,
                                                                   keys = rownames(vst.PRJNA314429.mtx),
                                                                   column = 'ENTREZID',
                                                                   keytype = 'SYMBOL',
                                                                   multiVals = first))
mdata.PRJNA314429$EndEst <- endest.results.PRJNA314429$estimated_time

ggplot(data = mdata.PRJNA314429, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.PRJNA314429[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
write.csv(mdata.PRJNA314429, file = 'processed.data/mdata.PRJNA314429.csv')
write.csv(vst.PRJNA314429.mtx, 'processed.data/vst.PRJNA314429.mtx.csv')



### GSE106602 

# Prepare metadata
mdata.GSE106602 <- read.csv('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE106602/SraRunTable.txt', header = T, sep = ',', row.names = 1)
head(mdata.GSE106602)
table(mdata.GSE106602$tissue_state, mdata.GSE106602$receptivity_time.point)
table(mdata.GSE106602$country_of_origin, mdata.GSE106602$receptivity_time.point)

mdata.GSE106602$names <- rownames(mdata.GSE106602)
mdata.GSE106602$GEO.accession <- 'GSE106602'
mdata.GSE106602$GSM.accession <- mdata.GSE106602$GEO_Accession..exp.
mdata.GSE106602$Sample.name <- mdata.GSE106602$donor
mdata.GSE106602$Phenotype <- mdata.GSE106602$tissue_state
mdata.GSE106602$Menstrual.cycle.time <- mdata.GSE106602$receptivity_time.point
mdata.GSE106602$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE106602$Age.years <- NA
mdata.GSE106602$files <- file.path('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE106602/quants/', mdata.GSE106602$names, 'quant.sf')
file.exists(mdata.GSE106602$files)
mdata.GSE106602$Batch <- mdata.GSE106602$Sequencing_batch
table(mdata.GSE106602$Phenotype, mdata.GSE106602$Menstrual.cycle.time)
table(subset(mdata.GSE106602, Menstrual.cycle.time == 'LH+7')$Batch)
mdata.GSE106602 <- subset(mdata.GSE106602, Menstrual.cycle.time == 'LH+7') # interesting samples are within one batch
mdata.GSE106602 <- mdata.GSE106602[ , c(9:17)]

table(mdata.GSE106602$Phenotype)
mdata.GSE106602 <- mdata.GSE106602 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('Healthy') ~ 'Fertile',
    Phenotype %in% c('IVF') ~ 'RIF'))

# Read counts and create DEseq2 object
se.GSE106602 <- tximeta(mdata.GSE106602)
gse.GSE106602 <- summarizeToGene(se.GSE106602)
gse.GSE106602 <- addIds(gse.GSE106602, 'SYMBOL', gene = TRUE)
ddsTC.GSE106602 <- DESeqDataSet(gse.GSE106602, design = ~ 1) # no adjustment of variation

# Keep features having maximum expression per known gene symbol identifier and at least 1 count in all samples
ddsTC.GSE106602 <- ddsTC.GSE106602[order(rowMeans(assay(ddsTC.GSE106602)), decreasing = TRUE), ]
ddsTC.GSE106602 <- ddsTC.GSE106602[
  !duplicated(rowData(ddsTC.GSE106602)$SYMBOL) & !(is.na(rowData(ddsTC.GSE106602)$SYMBOL)), ]
rownames(ddsTC.GSE106602) <- rowData(ddsTC.GSE106602)$SYMBOL
ddsTC.GSE106602 <- ddsTC.GSE106602[rowSums(counts(ddsTC.GSE106602) >= 1) >= ncol(counts(ddsTC.GSE106602)), ]

# EDA
vst.GSE106602.mtx <- assay(DESeq2::vst(ddsTC.GSE106602, blind = TRUE))
boxplot(vst.GSE106602.mtx)

hvg <- names(tail(sort(rowVars(vst.GSE106602.mtx)), 500))
pc <- prcomp(t(vst.GSE106602.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE106602$PC1 <- pca$PC1
mdata.GSE106602$PC2 <- pca$PC2

ggplot(data = mdata.GSE106602, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE106602[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE106602 <- estimate_cycle_time(exprs = vst.GSE106602.mtx,
                                                  entrez_ids = mapIds(org.Hs.eg.db,
                                                                      keys = rownames(vst.GSE106602.mtx),
                                                                      column = 'ENTREZID',
                                                                      keytype = 'SYMBOL',
                                                                      multiVals = first))
mdata.GSE106602$EndEst <- endest.results.GSE106602$estimated_time

ggplot(data = mdata.GSE106602, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE106602[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
write.csv(mdata.GSE106602, file = 'processed.data/mdata.GSE106602.csv')
write.csv(vst.GSE106602.mtx, 'processed.data/vst.GSE106602.mtx.csv')



### GSE205398 

# Prepare metadata
mdata.GSE205398 <- read.csv('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE205398/SraRunTable.txt', header = T, sep = ',', row.names = 1)
head(mdata.GSE205398)

mdata.GSE205398$names <- rownames(mdata.GSE205398)
mdata.GSE205398$GEO.accession <- 'GSE205398'
mdata.GSE205398$GSM.accession <- mdata.GSE205398$Sample.Name
mdata.GSE205398$Sample.name <- mdata.GSE205398$Library.Name
mdata.GSE205398$Phenotype <- mdata.GSE205398$diagnosis
mdata.GSE205398$Menstrual.cycle.time <- mdata.GSE205398$source_name
mdata.GSE205398$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE205398$Age.years <- NA
mdata.GSE205398$files <- file.path('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE205398/quants/', mdata.GSE205398$names, 'quant.sf')
file.exists(mdata.GSE205398$files)
mdata.GSE205398 <- mdata.GSE205398[ , c(29:37)]

table(mdata.GSE205398$Phenotype)
mdata.GSE205398 <- mdata.GSE205398 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('control') ~ 'Fertile',
    Phenotype %in% c('RIF') ~ 'RIF'))

# Read counts and create DEseq2 object
se.GSE205398 <- tximeta(mdata.GSE205398)
gse.GSE205398 <- summarizeToGene(se.GSE205398)
gse.GSE205398 <- addIds(gse.GSE205398, 'SYMBOL', gene = TRUE)
ddsTC.GSE205398 <- DESeqDataSet(gse.GSE205398, design = ~ 1) # no adjustment of variation

# Keep features having maximum expression per known gene symbol identifier and at least 1 count in all samples
ddsTC.GSE205398 <- ddsTC.GSE205398[order(rowMeans(assay(ddsTC.GSE205398)), decreasing = TRUE), ]
ddsTC.GSE205398 <- ddsTC.GSE205398[
  !duplicated(rowData(ddsTC.GSE205398)$SYMBOL) & !(is.na(rowData(ddsTC.GSE205398)$SYMBOL)), ]
rownames(ddsTC.GSE205398) <- rowData(ddsTC.GSE205398)$SYMBOL
ddsTC.GSE205398 <- ddsTC.GSE205398[rowSums(counts(ddsTC.GSE205398) >= 1) >= ncol(counts(ddsTC.GSE205398)), ]

# EDA
vst.GSE205398.mtx <- assay(DESeq2::vst(ddsTC.GSE205398, blind = TRUE))
boxplot(vst.GSE205398.mtx)

hvg <- names(tail(sort(rowVars(vst.GSE205398.mtx)), 500))
pc <- prcomp(t(vst.GSE205398.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE205398$PC1 <- pca$PC1
mdata.GSE205398$PC2 <- pca$PC2

ggplot(data = mdata.GSE205398, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE205398[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE205398 <- estimate_cycle_time(exprs = vst.GSE205398.mtx,
                                                entrez_ids = mapIds(org.Hs.eg.db,
                                                                    keys = rownames(vst.GSE205398.mtx),
                                                                    column = 'ENTREZID',
                                                                    keytype = 'SYMBOL',
                                                                    multiVals = first))
mdata.GSE205398$EndEst <- endest.results.GSE205398$estimated_time

ggplot(data = mdata.GSE205398, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE205398[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
write.csv(mdata.GSE205398, file = 'processed.data/mdata.GSE205398.csv')
write.csv(vst.GSE205398.mtx, 'processed.data/vst.GSE205398.mtx.csv')



### GSE207362 

# Prepare metadata
mdata.GSE207362 <- read.csv('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE207362/SraRunTable.txt', header = T, sep = ',', row.names = 1)
head(mdata.GSE207362)

mdata.GSE207362.1 <- read.csv('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE207362/SraRunTable.csv', header = T, sep = ',', row.names = 1)
head(mdata.GSE207362.1)

mdata.GSE207362$names <- rownames(mdata.GSE207362)
mdata.GSE207362$GEO.accession <- 'GSE207362'
mdata.GSE207362$GSM.accession <- mdata.GSE207362.1$Library.Name
mdata.GSE207362$Sample.name <- mdata.GSE207362.1$Sample_name
mdata.GSE207362$Phenotype <- mdata.GSE207362.1$Group
mdata.GSE207362$Menstrual.cycle.time <- mdata.GSE207362$genotype.treatment_condition
mdata.GSE207362$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE207362$Age.years <- NA
mdata.GSE207362$files <- file.path('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE207362/quants/', mdata.GSE207362$names, 'quant.sf')
file.exists(mdata.GSE207362$files)
mdata.GSE207362 <- mdata.GSE207362[ , c(33:41)]

table(mdata.GSE207362$Phenotype)
mdata.GSE207362 <- mdata.GSE207362 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('Control') ~ 'Fertile',
    Phenotype %in% c('RIF') ~ 'RIF'))

# Data availability statement https://www.nature.com/articles/s41419-023-05832-x#Sec10
# Among the uploaded human RNA-seq data, we only used the one from patients aged under 40 years without uterine diseases: 
# #2, 3, 5–7, 9–13, 16–18, 20, 21, 23, 24, 26, 31–33, 35, 36, 
# and 38–40 for fertile control; #1, 4–7, 9, 14, 17, 19, 21, 23, and 25 for RIF.

mdata.GSE207362 <- subset(mdata.GSE207362, Sample.name %in% c(
  sort(mdata.GSE207362$Sample.name)[c(39, 50, 63:65, 67, 29:32, 35:37, 40, 41, 43, 44, 46, 52:54, 56, 57, 59, 60, 62,
                                      1, 22:25, 27, 6, 9, 11, 14, 16, 18)]))
sort(mdata.GSE207362$Sample.name)
table(mdata.GSE207362$Phenotype)

# Read counts and create DEseq2 object
se.GSE207362 <- tximeta(mdata.GSE207362)
gse.GSE207362 <- summarizeToGene(se.GSE207362)
gse.GSE207362 <- addIds(gse.GSE207362, 'SYMBOL', gene = TRUE)
ddsTC.GSE207362 <- DESeqDataSet(gse.GSE207362, design = ~ 1) # no adjustment of variation

# Keep features having maximum expression per known gene symbol identifier and at least 1 count in all samples
ddsTC.GSE207362 <- ddsTC.GSE207362[order(rowMeans(assay(ddsTC.GSE207362)), decreasing = TRUE), ]
ddsTC.GSE207362 <- ddsTC.GSE207362[
  !duplicated(rowData(ddsTC.GSE207362)$SYMBOL) & !(is.na(rowData(ddsTC.GSE207362)$SYMBOL)), ]
rownames(ddsTC.GSE207362) <- rowData(ddsTC.GSE207362)$SYMBOL
ddsTC.GSE207362 <- ddsTC.GSE207362[rowSums(counts(ddsTC.GSE207362) >= 1) >= ncol(counts(ddsTC.GSE207362)), ]

# EDA
vst.GSE207362.mtx <- assay(DESeq2::vst(ddsTC.GSE207362, blind = TRUE))
boxplot(vst.GSE207362.mtx)

hvg <- names(tail(sort(rowVars(vst.GSE207362.mtx)), 500))
pc <- prcomp(t(vst.GSE207362.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE207362$PC1 <- pca$PC1
mdata.GSE207362$PC2 <- pca$PC2

ggplot(data = mdata.GSE207362, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE207362[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE207362 <- estimate_cycle_time(exprs = vst.GSE207362.mtx,
                                                entrez_ids = mapIds(org.Hs.eg.db,
                                                                    keys = rownames(vst.GSE207362.mtx),
                                                                    column = 'ENTREZID',
                                                                    keytype = 'SYMBOL',
                                                                    multiVals = first))
mdata.GSE207362$EndEst <- endest.results.GSE207362$estimated_time

ggplot(data = mdata.GSE207362, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE207362[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
write.csv(mdata.GSE207362, file = 'processed.data/mdata.GSE207362.csv')
write.csv(vst.GSE207362.mtx, 'processed.data/vst.GSE207362.mtx.csv')



## GSE243550

# Prepare metadata
mdata.GSE243550 <- read.csv('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE243550/SraRunTable.csv', header = T, sep = ',', row.names = 1)
head(mdata.GSE243550)

mdata.GSE243550$names <- rownames(mdata.GSE243550)
mdata.GSE243550$GEO.accession <- 'GSE243550'
mdata.GSE243550$GSM.accession <- mdata.GSE243550$Sample.Name
mdata.GSE243550$Sample.name <- mdata.GSE243550$Sample.Name
mdata.GSE243550$Phenotype <- mdata.GSE243550$condition
mdata.GSE243550$Menstrual.cycle.time <- mdata.GSE243550$timing_of_endometrial_biopsy
mdata.GSE243550$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE243550$Age.years <- NA
mdata.GSE243550$files <- file.path('/home/pd/datasets/hs_endometrium/bulk_rna_seq/GSE243550/quants/', mdata.GSE243550$names, 'quant.sf')
file.exists(mdata.GSE243550$files)
mdata.GSE243550 <- mdata.GSE243550[ , c(33:41)]

table(mdata.GSE243550$Phenotype)
mdata.GSE243550 <- mdata.GSE243550 %>%
  mutate(Phenotype = case_when(
    Phenotype %in% c('control') ~ 'Fertile',
    Phenotype %in% c('RIF') ~ 'RIF'))

# Read counts and create DEseq2 object
se.GSE243550 <- tximeta(mdata.GSE243550)
gse.GSE243550 <- summarizeToGene(se.GSE243550)
gse.GSE243550 <- addIds(gse.GSE243550, 'SYMBOL', gene = TRUE)
ddsTC.GSE243550 <- DESeqDataSet(gse.GSE243550, design = ~ 1) # no adjustment of variation

# Keep features having maximum expression per known gene symbol identifier and at least 1 count in all samples
ddsTC.GSE243550 <- ddsTC.GSE243550[order(rowMeans(assay(ddsTC.GSE243550)), decreasing = TRUE), ]
ddsTC.GSE243550 <- ddsTC.GSE243550[
  !duplicated(rowData(ddsTC.GSE243550)$SYMBOL) & !(is.na(rowData(ddsTC.GSE243550)$SYMBOL)), ]
rownames(ddsTC.GSE243550) <- rowData(ddsTC.GSE243550)$SYMBOL
ddsTC.GSE243550 <- ddsTC.GSE243550[rowSums(counts(ddsTC.GSE243550) >= 1) >= ncol(counts(ddsTC.GSE243550)), ]

# EDA
vst.GSE243550.mtx <- assay(DESeq2::vst(ddsTC.GSE243550, blind = TRUE))
boxplot(vst.GSE243550.mtx)

hvg <- names(tail(sort(rowVars(vst.GSE243550.mtx)), 500))
pc <- prcomp(t(vst.GSE243550.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pca <- as.data.frame(pc[5]$x)
mdata.GSE243550$PC1 <- pca$PC1
mdata.GSE243550$PC2 <- pca$PC2
# there is a possible outlier
ggplot(data = mdata.GSE243550, aes(x = PC1, y = PC2, color = Phenotype)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE243550[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'Phenotype')

# Estimate relative menstrual cycle progression
endest.results.GSE243550 <- estimate_cycle_time(exprs = vst.GSE243550.mtx,
                                                entrez_ids = mapIds(org.Hs.eg.db,
                                                                    keys = rownames(vst.GSE243550.mtx),
                                                                    column = 'ENTREZID',
                                                                    keytype = 'SYMBOL',
                                                                    multiVals = first))
mdata.GSE243550$EndEst <- endest.results.GSE243550$estimated_time

ggplot(data = mdata.GSE243550, aes(x = PC1, y = PC2, color = EndEst)) +
  geom_point(size=3) +
  labs(title = paste0('PCA ', mdata.GSE243550[1, 2]), 
       x = paste0('PC1: ', round((summary(pc))$importance[2,1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc))$importance[2,2], 4) * 100, '% variance'),
       color = 'EndEst') + scale_color_viridis_c()

# Save results
write.csv(mdata.GSE243550, file = 'processed.data/mdata.GSE243550.csv')
write.csv(vst.GSE243550.mtx, 'processed.data/vst.GSE243550.mtx.csv')







### Figure 1 - Explore GSE111974 dataset as an example

mdata.GSE111974 <- read.csv('processed.data/mdata.GSE111974.csv', sep = ',', header = T, row.names = 1)
exp.GSE111974 <- read.csv('processed.data/exp.GSE111974.max.csv', sep = ',', header = T, row.names = 1)

hvg.RIF.GSE111974 <- names(tail(sort(rowVars(as.matrix(exp.GSE111974))), 500))
pc.RIF.GSE111974 <- prcomp(t(exp.GSE111974[hvg.RIF.GSE111974, ]), center = TRUE, scale. = TRUE)
mdata.GSE111974$PC1 <- as.data.frame(pc.RIF.GSE111974[5]$x)$PC1
mdata.GSE111974$PC2 <- as.data.frame(pc.RIF.GSE111974[5]$x)$PC2

## Figure 1 A - PCA by EndEst and Phenotype
p1a <- ggplot(mdata.GSE111974, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE111974))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE111974))$importance[2, 2], 4) * 100, '% variance'),
       fill = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24))

## Figure 1 B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE111974)
p1b <- ggplot(mdata.GSE111974, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.659', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

## Figure 1 H - Dotplot by EndEst and Phenotype
p1h <- ggplot(mdata.GSE111974, aes(Phenotype, EndEst)) + 
  geom_point(position = position_dodge2(width = 0.3), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.3)) +
  geom_hline(yintercept = 31, color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = 85, color = 'black', linetype = 'dashed') +
  labs(title = '', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13)) +
  coord_flip(ylim = c(10, 88))

# Prepare subset of data to make groups comparable by EndEst
mdata.GSE111974.subset <- subset(mdata.GSE111974, EndEst > 31 & EndEst < 85)
exp.GSE111974.subset <- exp.GSE111974[, rownames(mdata.GSE111974.subset)]

# DEA using simple model using limma (9.2 Two Groups) on full data
Group <- factor(mdata.GSE111974$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(exp.GSE111974, design)
fit <- eBayes(fit)
head(design)
dea.GSE111974.simple.full <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# DEA using simple model using limma (9.2 Two Groups) on the subset of data
Group <- factor(mdata.GSE111974.subset$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(exp.GSE111974.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE111974.simple.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# DEA using mixed model using limma (9.6.2 Many time points) on the subset of data
X <- ns(mdata.GSE111974.subset$EndEst, df = 3)
Group <- factor(mdata.GSE111974.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(exp.GSE111974.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE111974.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# Identify representative genes
dea.GSE111974.simple.subset.up <- dea.GSE111974.simple.subset[
  dea.GSE111974.simple.subset$logFC > 0.667 & dea.GSE111974.simple.subset$P.Value < 0.05, ]
dea.GSE111974.simple.subset.down <- dea.GSE111974.simple.subset[
  dea.GSE111974.simple.subset$logFC < -0.667 & dea.GSE111974.simple.subset$P.Value < 0.05, ]
dea.GSE111974.complex.subset.up <- dea.GSE111974.complex.subset[
  dea.GSE111974.complex.subset$logFC > 0.667 & dea.GSE111974.complex.subset$P.Value < 0.05, ]
dea.GSE111974.complex.subset.down <- dea.GSE111974.complex.subset[
  dea.GSE111974.complex.subset$logFC < -0.667 & dea.GSE111974.complex.subset$P.Value < 0.05, ]

sum(rownames(dea.GSE111974.simple.subset.up) %in% rownames(dea.GSE111974.complex.subset.up)) # become correct with subsetting 
sum(rownames(dea.GSE111974.simple.subset.up) %in% rownames(dea.GSE111974.complex.subset.down))
sum(rownames(dea.GSE111974.simple.subset.down) %in% rownames(dea.GSE111974.complex.subset.down))
sum(rownames(dea.GSE111974.simple.subset.down) %in% rownames(dea.GSE111974.complex.subset.up))

dea.GSE111974.simple.subset.up[!(rownames(dea.GSE111974.simple.subset.up) %in% rownames(dea.GSE111974.complex.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in simple mode, CAPN6 is the top
dea.GSE111974.complex.subset.up[!(rownames(dea.GSE111974.complex.subset.up) %in% rownames(dea.GSE111974.simple.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in complex mode
dea.GSE111974.simple.subset.up[rownames(dea.GSE111974.simple.subset.up) %in% rownames(dea.GSE111974.complex.subset.up), ] %>% 
  arrange(-logFC) %>% head() # intersection, PHF8 is the top

dea.GSE111974.simple.subset.down[!(rownames(dea.GSE111974.simple.subset.down) %in% rownames(dea.GSE111974.complex.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in simple mode, PDZK1 is the top
dea.GSE111974.complex.subset.down[!(rownames(dea.GSE111974.complex.subset.down) %in% rownames(dea.GSE111974.simple.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in complex mode
dea.GSE111974.simple.subset.down[rownames(dea.GSE111974.simple.subset.down) %in% rownames(dea.GSE111974.complex.subset.down), ] %>% 
  arrange(logFC) %>% head() # intersection, SNORD12C is the top

mdata.GSE111974$CAPN6 <- as.numeric(exp.GSE111974['CAPN6', ])
mdata.GSE111974$PDZK1 <- as.numeric(exp.GSE111974['PDZK1', ])

mdata.GSE111974.subset$PHF8 <- as.numeric(exp.GSE111974.subset['PHF8', ])
mdata.GSE111974.subset$SNORD12C <- as.numeric(exp.GSE111974.subset['SNORD12C', ])

## Figure 1 C - Volcano for simple GEA on full dataset
dea.GSE111974.simple.full$symbol <- rownames(dea.GSE111974.simple.full)
dea.GSE111974.simple.full$delabel <- NA
dea.GSE111974.simple.full$delabel[dea.GSE111974.simple.full$symbol %in% c('CAPN6', 'PDZK1')] <- 
  dea.GSE111974.simple.full$symbol[dea.GSE111974.simple.full$symbol %in% c('CAPN6', 'PDZK1')]
p1c <- ggplot(data = dea.GSE111974.simple.full, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dea.GSE111974.simple.full[is.na(dea.GSE111974.simple.full$delabel), ], 
             color = 'black', fill = '#f0f0f0', size = 3, shape = 21, alpha = 0.8) +
  ggrepel::geom_label_repel(aes(label = delabel), fill = 'white', nudge_y = 1, 
                            segment.color = 'transparent', fontface = 'italic', size = 4.2) +
  geom_point(data = dea.GSE111974.simple.full[!(is.na(dea.GSE111974.simple.full$delabel)), ],
             color = 'black', fill = 'red', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = c(-0.667, 0.667), color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), color = 'black', linetype = 'dashed') +
  labs(title = 'DEA model:\nGene expression ~ 1 + Phenotype', x = 'log2 fold change', y = '-log10 p') + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))
  
## Figure 1 D - Down DEG for simple DEA on full dataset, boxplot
dea.GSE111974.simple.full['PDZK1',]
p1d <- ggplot(mdata.GSE111974, aes(Phenotype, PDZK1)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'PDZK1\np = 1.61e-3', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13))

## Figure 1 E - Down DEG for simple DEA on full dataset, EndEst splines
p1e <- ggplot(mdata.GSE111974, aes(EndEst, PDZK1, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'PDZK1', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

## Figure 1 F - Up DEG for simple DEA on full dataset, boxplot
dea.GSE111974.simple.full['CAPN6',]
p1f <- ggplot(mdata.GSE111974, aes(Phenotype, CAPN6)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'CAPN6\np = 1.03e-3', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13))

## Figure 1 G - Up DEG for simple DEA on full dataset, EndEst splines
p1g <- ggplot(mdata.GSE111974, aes(EndEst, CAPN6, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'CAPN6', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2') 

## Figure 1 I - Volcano for simple GEA on the subset
dea.GSE111974.simple.subset$symbol <- rownames(dea.GSE111974.simple.subset)
dea.GSE111974.simple.subset$delabel <- NA
dea.GSE111974.simple.subset$delabel[dea.GSE111974.simple.subset$symbol %in% c('PHF8', 'SNORD12C')] <- 
  dea.GSE111974.simple.subset$symbol[dea.GSE111974.simple.subset$symbol %in% c('PHF8', 'SNORD12C')]
p1i <- ggplot(data = dea.GSE111974.simple.subset, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dea.GSE111974.simple.subset[is.na(dea.GSE111974.simple.subset$delabel), ], 
             color = 'black', fill = '#f0f0f0', size = 3, shape = 21, alpha = 0.8) +
  ggrepel::geom_label_repel(aes(label = delabel), fill = 'white', nudge_y = 1, 
                            segment.color = 'transparent', fontface = 'italic', size = 4.2) +
  geom_point(data = dea.GSE111974.simple.subset[!(is.na(dea.GSE111974.simple.subset$delabel)), ],
             color = 'black', fill = 'red', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = c(-0.667, 0.667), color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), color = 'black', linetype = 'dashed') +
  labs(title = 'DEA model:\nGene expression ~ 1 + Phenotype', x = 'log2 fold change', y = '-log10 p') + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))

## Figure 1 J - Volcano for complex GEA on the subset
dea.GSE111974.complex.subset$symbol <- rownames(dea.GSE111974.complex.subset)
dea.GSE111974.complex.subset$delabel <- NA
dea.GSE111974.complex.subset$delabel[dea.GSE111974.complex.subset$symbol %in% c('PHF8', 'SNORD12C')] <- 
  dea.GSE111974.complex.subset$symbol[dea.GSE111974.complex.subset$symbol %in% c('PHF8', 'SNORD12C')]
p1j <- ggplot(data = dea.GSE111974.complex.subset, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(data = dea.GSE111974.complex.subset[is.na(dea.GSE111974.complex.subset$delabel), ], 
             color = 'black', fill = '#f0f0f0', size = 3, shape = 21, alpha = 0.8) +
  ggrepel::geom_label_repel(aes(label = delabel), fill = 'white', nudge_y = 1, 
                            segment.color = 'transparent', fontface = 'italic', size = 4.2) +
  geom_point(data = dea.GSE111974.complex.subset[!(is.na(dea.GSE111974.complex.subset$delabel)), ],
             color = 'black', fill = 'red', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = c(-0.667, 0.667), color = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), color = 'black', linetype = 'dashed') +
  labs(title = 'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst', x = 'log2 fold change', y = '-log10 p') + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13))

## Figure 1 K - Venn diagram
degs.GSE111974.simple.subset.up <- rownames(dea.GSE111974.simple.subset.up[
  !(rownames(dea.GSE111974.simple.subset.up) %in% rownames(dea.GSE111974.complex.subset.up)), ])
degs.GSE111974.complex.subset.up <- rownames(dea.GSE111974.complex.subset.up[
  !(rownames(dea.GSE111974.complex.subset.up) %in% rownames(dea.GSE111974.simple.subset.up)), ])
degs.GSE111974.intersection.subset.up <- rownames(dea.GSE111974.simple.subset.up[
  rownames(dea.GSE111974.simple.subset.up) %in% rownames(dea.GSE111974.complex.subset.up), ])
degs.GSE111974.simple.subset.down <- rownames(dea.GSE111974.simple.subset.down[
  !(rownames(dea.GSE111974.simple.subset.down) %in% rownames(dea.GSE111974.complex.subset.down)), ])
degs.GSE111974.complex.subset.down <- rownames(dea.GSE111974.complex.subset.down[
  !(rownames(dea.GSE111974.complex.subset.down) %in% rownames(dea.GSE111974.simple.subset.down)), ])
degs.GSE111974.intersection.subset.down <- rownames(dea.GSE111974.simple.subset.down[
  rownames(dea.GSE111974.simple.subset.down) %in% rownames(dea.GSE111974.complex.subset.down), ])

genesets.up <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                      c(degs.GSE111974.simple.subset.up, degs.GSE111974.intersection.subset.up),
                    'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst' = 
                      c(degs.GSE111974.complex.subset.up, degs.GSE111974.intersection.subset.up))
genesets.down <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                        c(degs.GSE111974.simple.subset.down, degs.GSE111974.intersection.subset.down),
                      'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst' = 
                        c(degs.GSE111974.complex.subset.down, degs.GSE111974.intersection.subset.down))

p1k1 <- ggVennDiagram(genesets.down,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.6, 0.01))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787') +
  scale_y_reverse() +
  annotate('text', label = c('DEA model:\nGene expression ~ 1 + Phenotype', 
                             'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst'), 
           x = c(-20, -20), y = c(-1.4, 5.4), size = 4.7)

p1k2 <- ggVennDiagram(genesets.up,
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.6, 0.01))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787') +
  scale_y_reverse() +
  annotate('text', label = c('DEA model:\nGene expression ~ 1 + Phenotype', 
                             'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst'), 
           x = c(-20, -20), y = c(-1.4, 5.4), size = 4.7)

## Figure 1 L - Down and Up DEGs for simple DEA on the subset, boxplots
dea.GSE111974.simple.subset['SNORD12C',]
p1l1 <- ggplot(mdata.GSE111974.subset, aes(Phenotype, SNORD12C)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'SNORD12C\np = 1.78e-13', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 13))

dea.GSE111974.simple.subset['PHF8',]
p1l2 <- ggplot(mdata.GSE111974.subset, aes(Phenotype, PHF8)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'PHF8\np = 9.11e-14', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13))

## Figure 1 K - Down and Up DEG for complex DEA on the subset, EndEst splines
dea.GSE111974.complex.subset['SNORD12C',]
p1m1 <- ggplot(mdata.GSE111974.subset, aes(EndEst, SNORD12C, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'SNORD12C\np = 8.55e-4', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'none') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2') 

dea.GSE111974.complex.subset['PHF8',]
p1m2 <- ggplot(mdata.GSE111974.subset, aes(EndEst, PHF8, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'PHF8\np = 1.03e-4', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2') 

## Figure 1 panel
p1.line1 <- wrap_plots(p1a, p1b, p1c, ncol = 3, widths = c(4, 1, 4))
p1.line2 <- wrap_plots(p1d, p1e, p1f, p1g, ncol = 4, widths = c(1, 3.5, 1, 3.5))
p1.line3 <- wrap_plots(p1h, p1i, p1j, ncol = 3, widths = c(2, 3.5, 3.5))
p1.line4 <- wrap_plots(p1k1, p1l1, p1m1, ncol = 3, widths = c(4.3, 1, 3.7))
p1.line5 <- wrap_plots(p1k2, p1l2, p1m2, ncol = 3, widths = c(4.3, 1, 3.7))
p1 <- wrap_plots(p1.line1, p1.line2, p1.line3, p1.line4, p1.line5, nrow = 5)
ggsave(plot = p1, filename = 'visualization/Figure 1.png', width = 13, height = 18, dpi = 300)





### Supplementary figure 1 - Volcano with top DEGs reported by the authors

# Replication of the original study analysis - Averaging probes signals for a gene and running DEA using simple model
exp.GSE111974.avrg <- read.csv('processed.data/exp.GSE111974.avrg.csv', sep = ',', header = T, row.names = 1)
Group <- factor(mdata.GSE111974$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(as.matrix(exp.GSE111974.avrg), design)
fit <- eBayes(fit)
head(design)
dea.GSE111974.simple.full.avrg <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.GSE111974.simple.full.avrg$symbol <- rownames(dea.GSE111974.simple.full.avrg)
dea.GSE111974.simple.full.avrg$delabel <- 'Not reported'
dea.GSE111974.simple.full.avrg$delabel[dea.GSE111974.simple.full.avrg$symbol %in% 
                                    c('PHF8', 'MARK2', 'SEPT9', 'RBM47', 'FAM21C', 'SMG5', 'LATS1', 'FAM107B', 
                                      'PAPOLA', 'ROCK2', 'HOOK3', 'MLL5', 'NCOA4', 'MAP4', 'CIZ1',
                                      'AGAP1', 'ACTB', 'MOB1A', 'USP33', 'TOP1', 'KIAA1429')] <- 'Up-regulated'
dea.GSE111974.simple.full.avrg$delabel[dea.GSE111974.simple.full.avrg$symbol %in% 
                                    c('XLOC_l2_015397', 'LOC100130557', 'EEF1A1', 'C1orf229')] <- 'Down-regulated'
dea.GSE111974.simple.full.avrg$delabelname <- NA
dea.GSE111974.simple.full.avrg$delabelname[dea.GSE111974.simple.full.avrg$symbol %in% 
                                        c('PHF8', 'MARK2', 'SEPT9', 'RBM47', 'FAM21C', 'SMG5', 'LATS1', 'FAM107B', 
                                          'PAPOLA', 'ROCK2', 'HOOK3', 'MLL5', 'NCOA4', 'MAP4', 'CIZ1',
                                          'AGAP1', 'ACTB', 'MOB1A', 'USP33', 'TOP1', 'KIAA1429',
                                          'XLOC_l2_015397', 'LOC100130557', 'EEF1A1', 'C1orf229')] <- 
  dea.GSE111974.simple.full.avrg$symbol[dea.GSE111974.simple.full.avrg$symbol %in% 
                                     c('PHF8', 'MARK2', 'SEPT9', 'RBM47', 'FAM21C', 'SMG5', 'LATS1', 'FAM107B', 
                                       'PAPOLA', 'ROCK2', 'HOOK3', 'MLL5', 'NCOA4', 'MAP4', 'CIZ1',
                                       'AGAP1', 'ACTB', 'MOB1A', 'USP33', 'TOP1', 'KIAA1429',
                                       'XLOC_l2_015397', 'LOC100130557', 'EEF1A1', 'C1orf229')]

ps1a <- ggplot(data = dea.GSE111974.simple.full.avrg, aes(x = logFC, y = -log10(P.Value), fill = delabel)) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  ggrepel::geom_label_repel(aes(label = delabelname, color = delabel), fill = 'white', nudge_y = 1,
                            segment.color = 'transparent', fontface = 'italic', size = 4.2, max.overlaps = 100) +
  # geom_vline(xintercept = c(-0.667, 0.667), color = 'black', linetype = 'dashed') +
  # geom_hline(yintercept = -log10(0.1), color = 'black', linetype = 'dashed') +
  labs(title = 'Replication of the original study analysis []:\n\nAveraging probes signals for a gene;\nDEA model: Gene expression ~ 1 + Phenotype', 
       x = 'log2 fold change', 
       y = '-log10 p',
       fill = 'DEGs reported in\nthe original study [ ]',
       color = 'DEGs reported in\nthe original study [ ]') + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  scale_fill_manual(
    values = c('#d73027', '#f0f0f0', '#4575b4'),
    breaks = c('Up-regulated', 'Not reported', 'Down-regulated')) +
  scale_color_manual(
    values = c('#d73027', '#f0f0f0', '#4575b4'),
    breaks = c('Up-regulated', 'Not reported', 'Down-regulated'))

# Analysis scheme applied in the present study - Selecting probes with the maximum signal for a gene and running DEA 
dea.GSE111974.simple.full.max <- dea.GSE111974.simple.full
dea.GSE111974.simple.full.max$delabel <- 'Not reported'
dea.GSE111974.simple.full.max$delabel[dea.GSE111974.simple.full.max$symbol %in% 
                                        c('PHF8', 'MARK2', 'SEPT9', 'RBM47', 'FAM21C', 'SMG5', 'LATS1', 'FAM107B', 
                                          'PAPOLA', 'ROCK2', 'HOOK3', 'MLL5', 'NCOA4', 'MAP4', 'CIZ1',
                                          'AGAP1', 'ACTB', 'MOB1A', 'USP33', 'TOP1', 'KIAA1429')] <- 'Up-regulated'
dea.GSE111974.simple.full.max$delabel[dea.GSE111974.simple.full.max$symbol %in% 
                                        c('XLOC_l2_015397', 'LOC100130557', 'EEF1A1', 'C1orf229')] <- 'Down-regulated'
dea.GSE111974.simple.full.max$delabelname <- NA
dea.GSE111974.simple.full.max$delabelname[dea.GSE111974.simple.full.max$symbol %in% 
                                            c('PHF8', 'MARK2', 'SEPT9', 'RBM47', 'FAM21C', 'SMG5', 'LATS1', 'FAM107B', 
                                              'PAPOLA', 'ROCK2', 'HOOK3', 'MLL5', 'NCOA4', 'MAP4', 'CIZ1',
                                              'AGAP1', 'ACTB', 'MOB1A', 'USP33', 'TOP1', 'KIAA1429',
                                              'XLOC_l2_015397', 'LOC100130557', 'EEF1A1', 'C1orf229')] <- 
  dea.GSE111974.simple.full.max$symbol[dea.GSE111974.simple.full.max$symbol %in% 
                                         c('PHF8', 'MARK2', 'SEPT9', 'RBM47', 'FAM21C', 'SMG5', 'LATS1', 'FAM107B', 
                                           'PAPOLA', 'ROCK2', 'HOOK3', 'MLL5', 'NCOA4', 'MAP4', 'CIZ1',
                                           'AGAP1', 'ACTB', 'MOB1A', 'USP33', 'TOP1', 'KIAA1429',
                                           'XLOC_l2_015397', 'LOC100130557', 'EEF1A1', 'C1orf229')]

ps1b <- ggplot(data = dea.GSE111974.simple.full.max, aes(x = logFC, y = -log10(P.Value), fill = delabel)) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  ggrepel::geom_label_repel(aes(label = delabelname, color = delabel), fill = 'white', nudge_y = 1,
                            segment.color = 'transparent', fontface = 'italic', size = 4.2, max.overlaps = 100) +
  # geom_vline(xintercept = c(-0.667, 0.667), color = 'black', linetype = 'dashed') +
  # geom_hline(yintercept = -log10(0.1), color = 'black', linetype = 'dashed') +
  labs(title = 'Analysis scheme applied in the present study:\n\nSelecting probes with the maximum signal for a gene;\nDEA model: Gene expression ~ 1 + Phenotype', 
       x = 'log2 fold change', 
       y = '-log10 p',
       fill = 'DEGs reported in\nthe original study [ ]',
       color = 'DEGs reported in\nthe original study [ ]') + 
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  scale_fill_manual(
    values = c('#d73027', '#f0f0f0', '#4575b4'),
    breaks = c('Up-regulated', 'Not reported', 'Down-regulated')) +
  scale_color_manual(
    values = c('#d73027', '#f0f0f0', '#4575b4'),
    breaks = c('Up-regulated', 'Not reported', 'Down-regulated'))

ps1 <- wrap_plots(ps1a, ps1b, nrow = 2)
ggsave(plot = ps1, filename = 'visualization/Figure S1.png', width = 10, height = 11, dpi = 300)






### Supplementary figure 2 - PCA and boxplots by EndEst and Phenotype for the remaining eleven GEO datasets

mdata.GSE103465 <- read.csv('processed.data/mdata.GSE103465.csv', sep = ',', header = T, row.names = 1)
mdata.GSE106602 <- read.csv('processed.data/mdata.GSE106602.csv', sep = ',', header = T, row.names = 1)
mdata.GSE188409 <- read.csv('processed.data/mdata.GSE188409.csv', sep = ',', header = T, row.names = 1)
mdata.GSE205398 <- read.csv('processed.data/mdata.GSE205398.csv', sep = ',', header = T, row.names = 1)
mdata.GSE207362 <- read.csv('processed.data/mdata.GSE207362.csv', sep = ',', header = T, row.names = 1)
mdata.GSE243550 <- read.csv('processed.data/mdata.GSE243550.csv', sep = ',', header = T, row.names = 1)
mdata.GSE26787 <- read.csv('processed.data/mdata.GSE26787.csv', sep = ',', header = T, row.names = 1)
mdata.GSE58144.subset <- read.csv('processed.data/mdata.GSE58144.subset.csv', sep = ',', header = T, row.names = 1)
mdata.GSE71331 <- read.csv('processed.data/mdata.GSE71331.csv', sep = ',', header = T, row.names = 1)
mdata.GSE92324 <- read.csv('processed.data/mdata.GSE92324.csv', sep = ',', header = T, row.names = 1)
mdata.PRJNA314429 <- read.csv('processed.data/mdata.PRJNA314429.csv', sep = ',', header = T, row.names = 1)

exp.GSE103465 <- read.csv('processed.data/exp.GSE103465.max.csv', sep = ',', header = T, row.names = 1)
exp.GSE106602 <- read.csv('processed.data/vst.GSE106602.mtx.csv', sep = ',', header = T, row.names = 1)
exp.GSE188409 <- read.csv('processed.data/exp.GSE188409.max.csv', sep = ',', header = T, row.names = 1)
exp.GSE205398 <- read.csv('processed.data/vst.GSE205398.mtx.csv', sep = ',', header = T, row.names = 1)
exp.GSE207362 <- read.csv('processed.data/vst.GSE207362.mtx.csv', sep = ',', header = T, row.names = 1)
exp.GSE243550 <- read.csv('processed.data/vst.GSE243550.mtx.csv', sep = ',', header = T, row.names = 1)
exp.GSE26787 <- read.csv('processed.data/exp.GSE26787.max.csv', sep = ',', header = T, row.names = 1)
exp.GSE58144.max.subset <- read.csv('processed.data/exp.GSE58144.max.subset.csv', sep = ',', header = T, row.names = 1)
exp.GSE71331 <- read.csv('processed.data/exp.GSE71331.max.csv', sep = ',', header = T, row.names = 1)
exp.GSE92324 <- read.csv('processed.data/exp.GSE92324.max.csv', sep = ',', header = T, row.names = 1)
exp.PRJNA314429 <- read.csv('processed.data/vst.PRJNA314429.mtx.csv', sep = ',', header = T, row.names = 1)


## GSE103465
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE103465))), 500))
pc.RIF.GSE103465 <- prcomp(t(exp.GSE103465[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE103465$PC1 <- as.data.frame(pc.RIF.GSE103465[5]$x)$PC1
mdata.GSE103465$PC2 <- as.data.frame(pc.RIF.GSE103465[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.1a <- ggplot(mdata.GSE103465, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE103465))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE103465))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE103465)
ps2.1b <- ggplot(mdata.GSE103465, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.793', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.1 <- wrap_plots(ps2.1a, ps2.1b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE103465', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE106602
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE106602))), 500))
pc.RIF.GSE106602 <- prcomp(t(exp.GSE106602[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE106602$PC1 <- as.data.frame(pc.RIF.GSE106602[5]$x)$PC1
mdata.GSE106602$PC2 <- as.data.frame(pc.RIF.GSE106602[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.2a <- ggplot(mdata.GSE106602, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE106602))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE106602))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 2), shape = guide_legend(order = 1))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE106602)
ps2.2b <- ggplot(mdata.GSE106602, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.026', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.2 <- wrap_plots(ps2.2a, ps2.2b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE106602', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE188409
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE188409))), 500))
pc.RIF.GSE188409 <- prcomp(t(exp.GSE188409[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE188409$PC1 <- as.data.frame(pc.RIF.GSE188409[5]$x)$PC1
mdata.GSE188409$PC2 <- as.data.frame(pc.RIF.GSE188409[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.3a <- ggplot(mdata.GSE188409, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE188409))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE188409))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE188409)
ps2.3b <- ggplot(mdata.GSE188409, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.922', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.3 <- wrap_plots(ps2.3a, ps2.3b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE188409', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE205398
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE205398))), 500))
pc.RIF.GSE205398 <- prcomp(t(exp.GSE205398[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE205398$PC1 <- as.data.frame(pc.RIF.GSE205398[5]$x)$PC1
mdata.GSE205398$PC2 <- as.data.frame(pc.RIF.GSE205398[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.4a <- ggplot(mdata.GSE205398, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE205398))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE205398))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE205398)
ps2.4b <- ggplot(mdata.GSE205398, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.456', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.4 <- wrap_plots(ps2.4a, ps2.4b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE205398', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE207362
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE207362))), 500))
pc.RIF.GSE207362 <- prcomp(t(exp.GSE207362[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE207362$PC1 <- as.data.frame(pc.RIF.GSE207362[5]$x)$PC1
mdata.GSE207362$PC2 <- as.data.frame(pc.RIF.GSE207362[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.5a <- ggplot(mdata.GSE207362, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE207362))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE207362))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE207362)
ps2.5b <- ggplot(mdata.GSE207362, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.301', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.5 <- wrap_plots(ps2.5a, ps2.5b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE207362', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE243550
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE243550))), 500))
pc.RIF.GSE243550 <- prcomp(t(exp.GSE243550[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE243550$PC1 <- as.data.frame(pc.RIF.GSE243550[5]$x)$PC1
mdata.GSE243550$PC2 <- as.data.frame(pc.RIF.GSE243550[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.6a <- ggplot(mdata.GSE243550, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE243550))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE243550))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE243550)
ps2.6b <- ggplot(mdata.GSE243550, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.142', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.6 <- wrap_plots(ps2.6a, ps2.6b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE243550', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE26787
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE26787))), 500))
pc.RIF.GSE26787 <- prcomp(t(exp.GSE26787[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE26787$PC1 <- as.data.frame(pc.RIF.GSE26787[5]$x)$PC1
mdata.GSE26787$PC2 <- as.data.frame(pc.RIF.GSE26787[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.7a <- ggplot(mdata.GSE26787, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE26787))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE26787))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE26787)
ps2.7b <- ggplot(mdata.GSE26787, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.015', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.7 <- wrap_plots(ps2.7a, ps2.7b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE26787', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE58144 - Cohort 1
table(mdata.GSE58144.subset$Batch)
mdata.GSE58144.subset.cohort1 <- subset(mdata.GSE58144.subset, Batch == 'Cohort 1')
exp.GSE58144.max.subset.cohort1 <- exp.GSE58144.max.subset[, rownames(mdata.GSE58144.subset.cohort1)]

hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE58144.max.subset.cohort1))), 500))
pc.RIF.GSE58144.subset.cohort1 <- prcomp(t(exp.GSE58144.max.subset.cohort1[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE58144.subset.cohort1$PC1 <- as.data.frame(pc.RIF.GSE58144.subset.cohort1[5]$x)$PC1
mdata.GSE58144.subset.cohort1$PC2 <- as.data.frame(pc.RIF.GSE58144.subset.cohort1[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.8a <- ggplot(mdata.GSE58144.subset.cohort1, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE58144.subset.cohort1))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE58144.subset.cohort1))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE58144.subset.cohort1)
ps2.8b <- ggplot(mdata.GSE58144.subset.cohort1, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'p = 0.750', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.8 <- wrap_plots(ps2.8a, ps2.8b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE58144 - Cohort 1', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE58144 - Cohort 2, batch 1
table(mdata.GSE58144.subset$Batch)
mdata.GSE58144.subset.cohort2.batch1 <- subset(mdata.GSE58144.subset, Batch == 'Cohort 2, batch 1')
exp.GSE58144.max.subset.cohort2.batch1 <- exp.GSE58144.max.subset[, rownames(mdata.GSE58144.subset.cohort2.batch1)]

hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE58144.max.subset.cohort2.batch1))), 500))
pc.RIF.GSE58144.subset.cohort2.batch1 <- prcomp(t(exp.GSE58144.max.subset.cohort2.batch1[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE58144.subset.cohort2.batch1$PC1 <- as.data.frame(pc.RIF.GSE58144.subset.cohort2.batch1[5]$x)$PC1
mdata.GSE58144.subset.cohort2.batch1$PC2 <- as.data.frame(pc.RIF.GSE58144.subset.cohort2.batch1[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.9a <- ggplot(mdata.GSE58144.subset.cohort2.batch1, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE58144.subset.cohort2.batch1))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE58144.subset.cohort2.batch1))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE58144.subset.cohort2.batch1)
ps2.9b <- ggplot(mdata.GSE58144.subset.cohort2.batch1, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'p = 0.950', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.9 <- wrap_plots(ps2.9a, ps2.9b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE58144 - Cohort 2, batch 1', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE71331
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE71331))), 500))
pc.RIF.GSE71331 <- prcomp(t(exp.GSE71331[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE71331$PC1 <- as.data.frame(pc.RIF.GSE71331[5]$x)$PC1
mdata.GSE71331$PC2 <- as.data.frame(pc.RIF.GSE71331[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.10a <- ggplot(mdata.GSE71331, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE71331))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE71331))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE71331)
ps2.10b <- ggplot(mdata.GSE71331, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.056', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.10 <- wrap_plots(ps2.10a, ps2.10b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE71331', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## GSE92324
hvg <- names(tail(sort(rowVars(as.matrix(exp.GSE92324))), 500))
pc.RIF.GSE92324 <- prcomp(t(exp.GSE92324[hvg, ]), center = TRUE, scale. = TRUE)
mdata.GSE92324$PC1 <- as.data.frame(pc.RIF.GSE92324[5]$x)$PC1
mdata.GSE92324$PC2 <- as.data.frame(pc.RIF.GSE92324[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.11a <- ggplot(mdata.GSE92324, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.GSE92324))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.GSE92324))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.GSE92324)
ps2.11b <- ggplot(mdata.GSE92324, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 4.55e-4', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.11 <- wrap_plots(ps2.11a, ps2.11b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'GSE92324', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## PRJNA314429
hvg <- names(tail(sort(rowVars(as.matrix(exp.PRJNA314429))), 500))
pc.RIF.PRJNA314429 <- prcomp(t(exp.PRJNA314429[hvg, ]), center = TRUE, scale. = TRUE)
mdata.PRJNA314429$PC1 <- as.data.frame(pc.RIF.PRJNA314429[5]$x)$PC1
mdata.PRJNA314429$PC2 <- as.data.frame(pc.RIF.PRJNA314429[5]$x)$PC2

# A - PCA by EndEst and Phenotype
ps2.12a <- ggplot(mdata.PRJNA314429, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.RIF.PRJNA314429))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.RIF.PRJNA314429))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24)) +
  guides(color = guide_legend(order = 1), shape = guide_legend(order = 2))

# B - Boxplot by EndEst and Phenotype
t.test(EndEst ~ Phenotype, mdata.PRJNA314429)
ps2.12b <- ggplot(mdata.PRJNA314429, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.442', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

ps2.12 <- wrap_plots(ps2.12a, ps2.12b, ncol = 2, widths = c(2, 1)) +
  plot_annotation(title = 'PRJNA314429', theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))


## Figure S2 panel
ps2 <- wrap_plots(ps2.1, ps2.2, ps2.3, ps2.4, ps2.5, ps2.6,
                  ps2.7, ps2.8, ps2.9, ps2.10, ps2.11, ps2.12, ncol = 3)
ggsave(plot = ps2, filename = 'visualization/Figure S2.png', width = 16, height = 12, dpi = 300)






### Figure 2 - Notable cases in analysis 

# Select five datasets with the largest sample size 
mdata.GSE103465 <- read.csv('processed.data/mdata.GSE103465.csv', sep = ',', header = T, row.names = 1)
mdata.GSE106602 <- read.csv('processed.data/mdata.GSE106602.csv', sep = ',', header = T, row.names = 1)[, -9]
mdata.GSE188409 <- read.csv('processed.data/mdata.GSE188409.csv', sep = ',', header = T, row.names = 1)
mdata.GSE205398 <- read.csv('processed.data/mdata.GSE205398.csv', sep = ',', header = T, row.names = 1)[, -9]
mdata.GSE207362 <- read.csv('processed.data/mdata.GSE207362.csv', sep = ',', header = T, row.names = 1)[, -9]
mdata.GSE243550 <- read.csv('processed.data/mdata.GSE243550.csv', sep = ',', header = T, row.names = 1)[, -9]
mdata.GSE26787 <- read.csv('processed.data/mdata.GSE26787.csv', sep = ',', header = T, row.names = 1)
mdata.GSE58144.subset <- read.csv('processed.data/mdata.GSE58144.subset.csv', sep = ',', header = T, row.names = 1)
mdata.GSE58144.subset.cohort1 <- subset(mdata.GSE58144.subset, Batch == 'Cohort 1')[, -9]
mdata.GSE58144.subset.cohort1$GEO.accession <- 'GSE58144 - Cohort 1'
mdata.GSE58144.subset.cohort2.batch1 <- subset(mdata.GSE58144.subset, Batch == 'Cohort 2, batch 1')[, -9]
mdata.GSE58144.subset.cohort2.batch1$GEO.accession <- 'GSE58144 - Cohort 2, batch 1'
mdata.GSE71331 <- read.csv('processed.data/mdata.GSE71331.csv', sep = ',', header = T, row.names = 1)
mdata.GSE92324 <- read.csv('processed.data/mdata.GSE92324.csv', sep = ',', header = T, row.names = 1)
mdata.PRJNA314429 <- read.csv('processed.data/mdata.PRJNA314429.csv', sep = ',', header = T, row.names = 1)[, -9]

mdata.merged <- rbind(mdata.GSE103465, mdata.GSE106602, mdata.GSE188409, mdata.GSE205398,
                      mdata.GSE207362, mdata.GSE243550, mdata.GSE26787, mdata.GSE58144.subset.cohort1,
                      mdata.GSE58144.subset.cohort2.batch1, mdata.GSE71331, mdata.GSE92324, mdata.PRJNA314429)
table(mdata.merged$GEO.accession, mdata.merged$Phenotype) # all samples

#                              Fertile RIF
# GSE103465                          3   3
# GSE106602                         16  19 *
# GSE188409                          5   5
# GSE205398                          3   3
# GSE207362                         26  12 *
# GSE243550                         20  20 *
# GSE26787                           5   5
# GSE58144 - Cohort 1               23  22 *
# GSE58144 - Cohort 2, batch 1      45  20 *
# GSE71331                           5   7
# GSE92324                           8  10 *
# PRJNA314429                        3   5

exp.GSE106602 <- read.csv('processed.data/vst.GSE106602.mtx.csv', sep = ',', header = T, row.names = 1)
exp.GSE207362 <- read.csv('processed.data/vst.GSE207362.mtx.csv', sep = ',', header = T, row.names = 1)
exp.GSE243550 <- read.csv('processed.data/vst.GSE243550.mtx.csv', sep = ',', header = T, row.names = 1)
exp.GSE58144.subset <- read.csv('processed.data/exp.GSE58144.max.subset.csv', sep = ',', header = T, row.names = 1)
exp.GSE92324 <- read.csv('processed.data/exp.GSE92324.max.csv', sep = ',', header = T, row.names = 1)


## GSE106602

# Figure 2A - Dotplot by EndEst and Phenotype GSE106602
p2a <- ggplot(mdata.GSE106602, aes(EndEst, Phenotype)) + 
  geom_point(position = position_dodge2(width = 0.3), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = 54, color = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 85, color = 'black', linetype = 'dashed') +
  labs(title = '', x = 'EndEst', y = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  scale_x_continuous(n.breaks = 6, limits = c(40, 85))

# Prepare subset of data to make groups comparable by EndEst
mdata.GSE106602.subset <- subset(mdata.GSE106602, EndEst > 54 & EndEst < 85)
exp.GSE106602.subset <- exp.GSE106602[, rownames(mdata.GSE106602.subset)]

# DEA using simple model using limma (9.2 Two Groups) on the subset of data
Group <- factor(mdata.GSE106602.subset$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(exp.GSE106602.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE106602.simple.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# DEA using mixed model using limma (9.6.2 Many time points) on the subset of data
X <- ns(mdata.GSE106602.subset$EndEst, df = 3)
Group <- factor(mdata.GSE106602.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(exp.GSE106602.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE106602.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# Identify representative genes
dea.GSE106602.simple.subset.up <- dea.GSE106602.simple.subset[
  dea.GSE106602.simple.subset$logFC > 0.667 & dea.GSE106602.simple.subset$P.Value < 0.05, ]
dea.GSE106602.simple.subset.down <- dea.GSE106602.simple.subset[
  dea.GSE106602.simple.subset$logFC < -0.667 & dea.GSE106602.simple.subset$P.Value < 0.05, ]
dea.GSE106602.complex.subset.up <- dea.GSE106602.complex.subset[
  dea.GSE106602.complex.subset$logFC > 0.667 & dea.GSE106602.complex.subset$P.Value < 0.05, ]
dea.GSE106602.complex.subset.down <- dea.GSE106602.complex.subset[
  dea.GSE106602.complex.subset$logFC < -0.667 & dea.GSE106602.complex.subset$P.Value < 0.05, ]

dea.GSE106602.simple.subset.up[!(rownames(dea.GSE106602.simple.subset.up) %in% rownames(dea.GSE106602.complex.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in simple mode
dea.GSE106602.complex.subset.up[!(rownames(dea.GSE106602.complex.subset.up) %in% rownames(dea.GSE106602.simple.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in complex mode
dea.GSE106602.simple.subset.up[rownames(dea.GSE106602.simple.subset.up) %in% rownames(dea.GSE106602.complex.subset.up), ] %>% 
  arrange(-logFC) %>% head() # intersection, take CLEC3B as an example

dea.GSE106602.simple.subset.down[!(rownames(dea.GSE106602.simple.subset.down) %in% rownames(dea.GSE106602.complex.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in simple mode
dea.GSE106602.complex.subset.down[!(rownames(dea.GSE106602.complex.subset.down) %in% rownames(dea.GSE106602.simple.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in complex mode
dea.GSE106602.simple.subset.down[rownames(dea.GSE106602.simple.subset.down) %in% rownames(dea.GSE106602.complex.subset.down), ] %>% 
  arrange(logFC) %>% head() # intersection, SRCAP is the top

mdata.GSE106602.subset$CLEC3B <- as.numeric(exp.GSE106602.subset['CLEC3B', ])

# Figure 2B - Venn diagrams for DEGs GSE106602
degs.GSE106602.simple.subset.up <- rownames(dea.GSE106602.simple.subset.up[
  !(rownames(dea.GSE106602.simple.subset.up) %in% rownames(dea.GSE106602.complex.subset.up)), ])
degs.GSE106602.complex.subset.up <- rownames(dea.GSE106602.complex.subset.up[
  !(rownames(dea.GSE106602.complex.subset.up) %in% rownames(dea.GSE106602.simple.subset.up)), ])
degs.GSE106602.intersection.subset.up <- rownames(dea.GSE106602.simple.subset.up[
  rownames(dea.GSE106602.simple.subset.up) %in% rownames(dea.GSE106602.complex.subset.up), ])
degs.GSE106602.simple.subset.down <- rownames(dea.GSE106602.simple.subset.down[
  !(rownames(dea.GSE106602.simple.subset.down) %in% rownames(dea.GSE106602.complex.subset.down)), ])
degs.GSE106602.complex.subset.down <- rownames(dea.GSE106602.complex.subset.down[
  !(rownames(dea.GSE106602.complex.subset.down) %in% rownames(dea.GSE106602.simple.subset.down)), ])
degs.GSE106602.intersection.subset.down <- rownames(dea.GSE106602.simple.subset.down[
  rownames(dea.GSE106602.simple.subset.down) %in% rownames(dea.GSE106602.complex.subset.down), ])

genesets.down <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                        c(degs.GSE106602.simple.subset.down, degs.GSE106602.intersection.subset.down),
                      'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                        c(degs.GSE106602.complex.subset.down, degs.GSE106602.intersection.subset.down))
genesets.up <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                      c(degs.GSE106602.simple.subset.up, degs.GSE106602.intersection.subset.up),
                    'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                      c(degs.GSE106602.complex.subset.up, degs.GSE106602.intersection.subset.up))

p2b1 <- ggVennDiagram(genesets.down,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.6, 0.01))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 850)) +
  scale_y_reverse() +
  annotate('text', label = c('DEA model:\nGene expression ~ 1 + Phenotype', 
                             'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst'), 
           x = c(-20, -20), y = c(-1.4, 5.4), size = 4.7)

p2b2 <- ggVennDiagram(genesets.up,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number\nof genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 850)) +
  scale_y_reverse()

p2b <- wrap_plots(p2b1, p2b2, ncol = 2)

# Figure 2 C - Top Up DEG GSE106602, boxplot and curves
dea.GSE106602.simple.subset['CLEC3B',]
p2c1 <- ggplot(mdata.GSE106602.subset, aes(Phenotype, CLEC3B)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'CLEC3B\np = 2.20e-07', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13))

dea.GSE106602.complex.subset['CLEC3B',]
p2c2 <- ggplot(mdata.GSE106602.subset, aes(EndEst, CLEC3B, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'CLEC3B\np = 7.14e-3', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

p2c <- wrap_plots(p2c1, p2c2, ncol = 2, widths = c(1, 3))



## GSE92324

# Figure 2D - Dotplot by EndEst and Phenotype GSE92324
p2d <- ggplot(mdata.GSE92324, aes(EndEst, Phenotype)) + 
  geom_point(position = position_dodge2(width = 0.3), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = 63, color = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 88, color = 'black', linetype = 'dashed') +
  labs(title = '', x = 'EndEst', y = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  scale_x_continuous(n.breaks = 6, limits = c(63, 88))

# Prepare subset of data to make groups comparable by EndEst
mdata.GSE92324.subset <- mdata.GSE92324 # no filtering for this dataset
exp.GSE92324.subset <- exp.GSE92324[, rownames(mdata.GSE92324.subset)]

# DEA using simple model using limma (9.2 Two Groups) on the subset of data
Group <- factor(mdata.GSE92324.subset$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(exp.GSE92324.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE92324.simple.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# DEA using mixed model using limma (9.6.2 Many time points) on the subset of data
X <- ns(mdata.GSE92324.subset$EndEst, df = 3)
Group <- factor(mdata.GSE92324.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(exp.GSE92324.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE92324.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# Identify representative genes
dea.GSE92324.simple.subset.up <- dea.GSE92324.simple.subset[
  dea.GSE92324.simple.subset$logFC > 0.667 & dea.GSE92324.simple.subset$P.Value < 0.05, ]
dea.GSE92324.simple.subset.down <- dea.GSE92324.simple.subset[
  dea.GSE92324.simple.subset$logFC < -0.667 & dea.GSE92324.simple.subset$P.Value < 0.05, ]
dea.GSE92324.complex.subset.up <- dea.GSE92324.complex.subset[
  dea.GSE92324.complex.subset$logFC > 0.667 & dea.GSE92324.complex.subset$P.Value < 0.05, ]
dea.GSE92324.complex.subset.down <- dea.GSE92324.complex.subset[
  dea.GSE92324.complex.subset$logFC < -0.667 & dea.GSE92324.complex.subset$P.Value < 0.05, ]

dea.GSE92324.simple.subset.up[!(rownames(dea.GSE92324.simple.subset.up) %in% rownames(dea.GSE92324.complex.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in simple mode
dea.GSE92324.complex.subset.up[!(rownames(dea.GSE92324.complex.subset.up) %in% rownames(dea.GSE92324.simple.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in complex mode
dea.GSE92324.simple.subset.up[rownames(dea.GSE92324.simple.subset.up) %in% rownames(dea.GSE92324.complex.subset.up), ] %>% 
  arrange(-logFC) %>% head() # intersection, take TMEM144 as an example

dea.GSE92324.simple.subset.down[!(rownames(dea.GSE92324.simple.subset.down) %in% rownames(dea.GSE92324.complex.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in simple mode
dea.GSE92324.complex.subset.down[!(rownames(dea.GSE92324.complex.subset.down) %in% rownames(dea.GSE92324.simple.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in complex mode
dea.GSE92324.simple.subset.down[rownames(dea.GSE92324.simple.subset.down) %in% rownames(dea.GSE92324.complex.subset.down), ] %>% 
  arrange(logFC) %>% head() # intersection

mdata.GSE92324.subset$TMEM144 <- as.numeric(exp.GSE92324.subset['TMEM144', ])

# Figure 2E - Venn diagrams for DEGs GSE92324
degs.GSE92324.simple.subset.up <- rownames(dea.GSE92324.simple.subset.up[
  !(rownames(dea.GSE92324.simple.subset.up) %in% rownames(dea.GSE92324.complex.subset.up)), ])
degs.GSE92324.complex.subset.up <- rownames(dea.GSE92324.complex.subset.up[
  !(rownames(dea.GSE92324.complex.subset.up) %in% rownames(dea.GSE92324.simple.subset.up)), ])
degs.GSE92324.intersection.subset.up <- rownames(dea.GSE92324.simple.subset.up[
  rownames(dea.GSE92324.simple.subset.up) %in% rownames(dea.GSE92324.complex.subset.up), ])
degs.GSE92324.simple.subset.down <- rownames(dea.GSE92324.simple.subset.down[
  !(rownames(dea.GSE92324.simple.subset.down) %in% rownames(dea.GSE92324.complex.subset.down)), ])
degs.GSE92324.complex.subset.down <- rownames(dea.GSE92324.complex.subset.down[
  !(rownames(dea.GSE92324.complex.subset.down) %in% rownames(dea.GSE92324.simple.subset.down)), ])
degs.GSE92324.intersection.subset.down <- rownames(dea.GSE92324.simple.subset.down[
  rownames(dea.GSE92324.simple.subset.down) %in% rownames(dea.GSE92324.complex.subset.down), ])

genesets.down <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                        c(degs.GSE92324.simple.subset.down, degs.GSE92324.intersection.subset.down),
                      'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                        c(degs.GSE92324.complex.subset.down, degs.GSE92324.intersection.subset.down))
genesets.up <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                      c(degs.GSE92324.simple.subset.up, degs.GSE92324.intersection.subset.up),
                    'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                      c(degs.GSE92324.complex.subset.up, degs.GSE92324.intersection.subset.up))

p2e1 <- ggVennDiagram(genesets.down,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.6, 0.01))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 1050)) +
  scale_y_reverse() +
  annotate('text', label = c('DEA model:\nGene expression ~ 1 + Phenotype', 
                             'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst'), 
           x = c(-20, -20), y = c(-1.4, 5.4), size = 4.7)

p2e2 <- ggVennDiagram(genesets.up,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number\nof genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 1050)) +
  scale_y_reverse()

p2e <- wrap_plots(p2e1, p2e2, ncol = 2)

# Figure 2 F - Top Up DEG GSE92324, boxplot and curves
dea.GSE92324.simple.subset['TMEM144',]
p2f1 <- ggplot(mdata.GSE92324.subset, aes(Phenotype, TMEM144)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'TMEM144\np = 1.25e-08', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13))

dea.GSE92324.complex.subset['TMEM144',]
p2f2 <- ggplot(mdata.GSE92324.subset, aes(EndEst, TMEM144, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'TMEM144\np = 0.016', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

p2f <- wrap_plots(p2f1, p2f2, ncol = 2, widths = c(1, 3))



## GSE243550

# Figure 2G - Dotplot by EndEst and Phenotype GSE243550
p2g1 <- ggplot(mdata.GSE243550, aes(EndEst, Phenotype)) + 
  geom_point(position = position_dodge2(width = 0.3), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = 66.5, color = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 80, color = 'black', linetype = 'dashed') +
  labs(title = '', x = 'EndEst', y = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  scale_x_continuous(n.breaks = 2, limits = c(0, 10))

p2g2 <- ggplot(mdata.GSE243550, aes(EndEst, Phenotype)) + 
  geom_point(position = position_dodge2(width = 0.3), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = 66.5, color = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 79.5, color = 'black', linetype = 'dashed') +
  labs(title = '', x = 'EndEst', y = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 13),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  xlim(c(65, 83))

p2g <- wrap_plots(p2g1, p2g2, widths = c(1, 4), ncol = 2)

# Prepare subset of data to make groups comparable by EndEst 
mdata.GSE243550.subset <- subset(mdata.GSE243550, EndEst > 66 & EndEst < 80) 
exp.GSE243550.subset <- exp.GSE243550[, rownames(mdata.GSE243550.subset)]
# More soft filtering by EndEst > 62 & EndEst < 84 enable identification of a larger number of DEGs,
# however all of this DEGs would be artefacts of unequal coverage of samples at both ends of the EndEst scores range

# DEA using simple model using limma (9.2 Two Groups) on the subset of data
Group <- factor(mdata.GSE243550.subset$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(exp.GSE243550.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE243550.simple.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# DEA using mixed model using limma (9.6.2 Many time points) on the subset of data
X <- ns(mdata.GSE243550.subset$EndEst, df = 3)
Group <- factor(mdata.GSE243550.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(exp.GSE243550.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE243550.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# Identify representative genes
dea.GSE243550.simple.subset.up <- dea.GSE243550.simple.subset[
  dea.GSE243550.simple.subset$logFC > 0.667 & dea.GSE243550.simple.subset$P.Value < 0.05, ]
dea.GSE243550.simple.subset.down <- dea.GSE243550.simple.subset[
  dea.GSE243550.simple.subset$logFC < -0.667 & dea.GSE243550.simple.subset$P.Value < 0.05, ]
dea.GSE243550.complex.subset.up <- dea.GSE243550.complex.subset[
  dea.GSE243550.complex.subset$logFC > 0.667 & dea.GSE243550.complex.subset$P.Value < 0.05, ]
dea.GSE243550.complex.subset.down <- dea.GSE243550.complex.subset[
  dea.GSE243550.complex.subset$logFC < -0.667 & dea.GSE243550.complex.subset$P.Value < 0.05, ]

dea.GSE243550.simple.subset.up[!(rownames(dea.GSE243550.simple.subset.up) %in% rownames(dea.GSE243550.complex.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in simple mode
dea.GSE243550.complex.subset.up[!(rownames(dea.GSE243550.complex.subset.up) %in% rownames(dea.GSE243550.simple.subset.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in complex mode
dea.GSE243550.simple.subset.up[rownames(dea.GSE243550.simple.subset.up) %in% rownames(dea.GSE243550.complex.subset.up), ] %>% 
  arrange(-logFC) %>% head() # intersection, TNC (!)

dea.GSE243550.simple.subset.down[!(rownames(dea.GSE243550.simple.subset.down) %in% rownames(dea.GSE243550.complex.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in simple mode
dea.GSE243550.complex.subset.down[!(rownames(dea.GSE243550.complex.subset.down) %in% rownames(dea.GSE243550.simple.subset.down)), ] %>% 
  arrange(logFC) %>% head() # down only in complex mode
dea.GSE243550.simple.subset.down[rownames(dea.GSE243550.simple.subset.down) %in% rownames(dea.GSE243550.complex.subset.down), ] %>% 
  arrange(logFC) %>% head() # intersection, take COL27A1 as an example

mdata.GSE243550.subset$COL27A1 <- as.numeric(exp.GSE243550.subset['COL27A1', ])

# Figure 2H - Venn diagrams for DEGs GSE243550
degs.GSE243550.simple.subset.up <- rownames(dea.GSE243550.simple.subset.up[
  !(rownames(dea.GSE243550.simple.subset.up) %in% rownames(dea.GSE243550.complex.subset.up)), ])
degs.GSE243550.complex.subset.up <- rownames(dea.GSE243550.complex.subset.up[
  !(rownames(dea.GSE243550.complex.subset.up) %in% rownames(dea.GSE243550.simple.subset.up)), ])
degs.GSE243550.intersection.subset.up <- rownames(dea.GSE243550.simple.subset.up[
  rownames(dea.GSE243550.simple.subset.up) %in% rownames(dea.GSE243550.complex.subset.up), ])
degs.GSE243550.simple.subset.down <- rownames(dea.GSE243550.simple.subset.down[
  !(rownames(dea.GSE243550.simple.subset.down) %in% rownames(dea.GSE243550.complex.subset.down)), ])
degs.GSE243550.complex.subset.down <- rownames(dea.GSE243550.complex.subset.down[
  !(rownames(dea.GSE243550.complex.subset.down) %in% rownames(dea.GSE243550.simple.subset.down)), ])
degs.GSE243550.intersection.subset.down <- rownames(dea.GSE243550.simple.subset.down[
  rownames(dea.GSE243550.simple.subset.down) %in% rownames(dea.GSE243550.complex.subset.down), ])

genesets.down <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                        c(degs.GSE243550.simple.subset.down, degs.GSE243550.intersection.subset.down),
                      'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                        c(degs.GSE243550.complex.subset.down, degs.GSE243550.intersection.subset.down))
genesets.up <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                      c(degs.GSE243550.simple.subset.up, degs.GSE243550.intersection.subset.up),
                    'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                      c(degs.GSE243550.complex.subset.up, degs.GSE243550.intersection.subset.up))

p2h1 <- ggVennDiagram(genesets.down,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.6, 0.01))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 103)) +
  scale_y_reverse() +
  annotate('text', label = c('DEA model:\nGene expression ~ 1 + Phenotype', 
                             'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst'), 
           x = c(-20, -20), y = c(-1.4, 5.4), size = 4.7)

p2h2 <- ggVennDiagram(genesets.up,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number\nof genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 103)) +
  scale_y_reverse()

p2h <- wrap_plots(p2h1, p2h2, ncol = 2)

# Figure 2 I - Top Up DEG GSE243550, boxplot and curves
dea.GSE243550.simple.subset['COL27A1', ]
p2i1 <- ggplot(mdata.GSE243550.subset, aes(Phenotype, COL27A1)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'COL27A1\np = 4.87e-05', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13))

dea.GSE243550.complex.subset['COL27A1',]
p2i2 <- ggplot(mdata.GSE243550.subset, aes(EndEst, COL27A1, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'COL27A1\np = 0.029', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

p2i <- wrap_plots(p2i1, p2i2, ncol = 2, widths = c(1, 3))



## GSE207362

# Figure 2J - Dotplot by EndEst and Phenotype GSE207362
p2j <- ggplot(mdata.GSE207362, aes(EndEst, Phenotype)) + 
  geom_point(position = position_dodge2(width = 0.3), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = 72.5, color = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 82.5, color = 'black', linetype = 'dashed') +
  labs(title = '', x = 'EndEst', y = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  scale_x_continuous(n.breaks = 6, limits = c(70, 91))

# Prepare subset of data to make groups comparable by EndEst
mdata.GSE207362.subset <- subset(mdata.GSE207362, EndEst > 72 & EndEst < 83)
exp.GSE207362.subset <- exp.GSE207362[, rownames(mdata.GSE207362.subset)]

# DEA using simple model using limma (9.2 Two Groups) on the subset of data
Group <- factor(mdata.GSE207362.subset$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(exp.GSE207362.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE207362.simple.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# DEA using mixed model using limma (9.6.2 Many time points) on the subset of data
X <- ns(mdata.GSE207362.subset$EndEst, df = 3)
Group <- factor(mdata.GSE207362.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(exp.GSE207362.subset, design)
fit <- eBayes(fit)
head(design)
dea.GSE207362.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# Identify representative genes
dea.GSE207362.simple.subset.up <- dea.GSE207362.simple.subset[
  dea.GSE207362.simple.subset$logFC > 0.667 & dea.GSE207362.simple.subset$P.Value < 0.05, ]
dea.GSE207362.simple.subset.down <- dea.GSE207362.simple.subset[
  dea.GSE207362.simple.subset$logFC < -0.667 & dea.GSE207362.simple.subset$P.Value < 0.05, ]
dea.GSE207362.complex.subset.up <- dea.GSE207362.complex.subset[
  dea.GSE207362.complex.subset$logFC > 0.667 & dea.GSE207362.complex.subset$P.Value < 0.05, ]
dea.GSE207362.complex.subset.down <- dea.GSE207362.complex.subset[
  dea.GSE207362.complex.subset$logFC < -0.667 & dea.GSE207362.complex.subset$P.Value < 0.05, ]

dea.GSE207362.simple.subset.up[!(rownames(dea.GSE207362.simple.subset.up) %in% rownames(dea.GSE207362.complex.subset.up)), ] %>%
  arrange(-logFC) %>% head() # up only in simple mode
dea.GSE207362.complex.subset.up[!(rownames(dea.GSE207362.complex.subset.up) %in% rownames(dea.GSE207362.simple.subset.up)), ] %>%
  arrange(-logFC) %>% head() # up only in complex mode
dea.GSE207362.simple.subset.up[rownames(dea.GSE207362.simple.subset.up) %in% rownames(dea.GSE207362.complex.subset.up), ] %>%
  arrange(-logFC) %>% head() # intersection - take CPZ 

dea.GSE207362.simple.subset.down[!(rownames(dea.GSE207362.simple.subset.down) %in% rownames(dea.GSE207362.complex.subset.down)), ] %>%
  arrange(logFC) %>% head() # down only in simple mode
dea.GSE207362.complex.subset.down[!(rownames(dea.GSE207362.complex.subset.down) %in% rownames(dea.GSE207362.simple.subset.down)), ] %>%
  arrange(logFC) %>% head() # down only in complex mode
dea.GSE207362.simple.subset.down[rownames(dea.GSE207362.simple.subset.down) %in% rownames(dea.GSE207362.complex.subset.down), ] %>%
  arrange(logFC) %>% head() # intersection

mdata.GSE207362.subset$CPZ <- as.numeric(exp.GSE207362.subset['CPZ', ])

# Figure 2K - Venn diagrams for DEGs GSE207362
degs.GSE207362.simple.subset.up <- rownames(dea.GSE207362.simple.subset.up[
  !(rownames(dea.GSE207362.simple.subset.up) %in% rownames(dea.GSE207362.complex.subset.up)), ])
degs.GSE207362.complex.subset.up <- rownames(dea.GSE207362.complex.subset.up[
  !(rownames(dea.GSE207362.complex.subset.up) %in% rownames(dea.GSE207362.simple.subset.up)), ])
degs.GSE207362.intersection.subset.up <- rownames(dea.GSE207362.simple.subset.up[
  rownames(dea.GSE207362.simple.subset.up) %in% rownames(dea.GSE207362.complex.subset.up), ])
degs.GSE207362.simple.subset.down <- rownames(dea.GSE207362.simple.subset.down[
  !(rownames(dea.GSE207362.simple.subset.down) %in% rownames(dea.GSE207362.complex.subset.down)), ])
degs.GSE207362.complex.subset.down <- rownames(dea.GSE207362.complex.subset.down[
  !(rownames(dea.GSE207362.complex.subset.down) %in% rownames(dea.GSE207362.simple.subset.down)), ])
degs.GSE207362.intersection.subset.down <- rownames(dea.GSE207362.simple.subset.down[
  rownames(dea.GSE207362.simple.subset.down) %in% rownames(dea.GSE207362.complex.subset.down), ])

genesets.down <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                        c(degs.GSE207362.simple.subset.down, degs.GSE207362.intersection.subset.down),
                      'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                        c(degs.GSE207362.complex.subset.down, degs.GSE207362.intersection.subset.down))
genesets.up <- list('DEA model:\nGene expression ~ 1 + Phenotype' = 
                      c(degs.GSE207362.simple.subset.up, degs.GSE207362.intersection.subset.up),
                    'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst' = 
                      c(degs.GSE207362.complex.subset.up, degs.GSE207362.intersection.subset.up))

p2k1 <- ggVennDiagram(genesets.down,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.6, 0.01))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 631)) +
  scale_y_reverse() +
  annotate('text', label = c('DEA model:\nGene expression ~ 1 + Phenotype', 
                             'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst'), 
           x = c(-20, -20), y = c(-1.4, 5.4), size = 4.7)

p2k2 <- ggVennDiagram(genesets.up,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number\nof genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 631)) +
  scale_y_reverse()

p2k <- wrap_plots(p2k1, p2k2, ncol = 2)


# Figure 2 L - Top Up DEG GSE243550, boxplot and curves
dea.GSE207362.simple.subset['CPZ', ]
p2l1 <- ggplot(mdata.GSE207362.subset, aes(Phenotype, CPZ)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'CPZ\np = 6.56e-03', x = 'Phenotype', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13))

dea.GSE207362.complex.subset['CPZ',]
p2l2 <- ggplot(mdata.GSE207362.subset, aes(EndEst, CPZ, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'CPZ\np = 0.046', x = 'EndEst', y = 'mRNA levels,\nlog-normalized values') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

p2l <- wrap_plots(p2l1, p2l2, ncol = 2, widths = c(1, 3))



## GSE58144

# Figure 2M - Dotplot by EndEst and Phenotype GSE58144

mdata.GSE58144.subset$Phenotype.Batch <- paste0(mdata.GSE58144.subset$Phenotype, 
                                                ' - ',
                                                mdata.GSE58144.subset$Batch)
table(mdata.GSE58144.subset$Phenotype.Batch)
mdata.GSE58144.subset$Phenotype.Batch <- factor(mdata.GSE58144.subset$Phenotype.Batch, levels = c(
  'Fertile - Cohort 2, batch 1', 'RIF - Cohort 2, batch 1', 'Fertile - Cohort 1', 'RIF - Cohort 1'))

p2m <- ggplot(mdata.GSE58144.subset, aes(EndEst, Phenotype.Batch)) + 
  geom_point(position = position_dodge2(width = 0.3), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  geom_vline(xintercept = 58.5, color = 'black', linetype = 'dashed') +
  geom_vline(xintercept = 75.5, color = 'black', linetype = 'dashed') +
  labs(title = '', x = 'EndEst', y = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  scale_x_continuous(n.breaks = 6, limits = c(6, 84))

# Prepare subset of data to make groups comparable by EndEst
mdata.GSE58144.subset2 <- subset(mdata.GSE58144.subset, EndEst > 58 & EndEst < 76)
exp.GSE58144.subset2 <- exp.GSE58144.subset[, rownames(mdata.GSE58144.subset2)]

# DEA using simple model using limma (9.2 Two Groups) on the subset of data
Group <- factor(mdata.GSE58144.subset2$Phenotype)
Batch <- factor(mdata.GSE58144.subset2$Batch)
design <- model.matrix(~ Group + Batch)
fit <- lmFit(exp.GSE58144.subset2, design)
fit <- eBayes(fit)
head(design)
dea.GSE58144.simple.subset2 <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# DEA using mixed model using limma (9.6.2 Many time points) on the subset of data
X <- ns(mdata.GSE58144.subset2$EndEst, df = 3)
design <- model.matrix(~ Group * X + Batch)
fit <- lmFit(exp.GSE58144.subset2, design)
fit <- eBayes(fit)
head(design)
dea.GSE58144.complex.subset2 <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)

# Identify representative genes
dea.GSE58144.simple.subset2.up <- dea.GSE58144.simple.subset2[
  dea.GSE58144.simple.subset2$logFC > 0.667 & dea.GSE58144.simple.subset2$P.Value < 0.05, ]
dea.GSE58144.simple.subset2.down <- dea.GSE58144.simple.subset2[
  dea.GSE58144.simple.subset2$logFC < -0.667 & dea.GSE58144.simple.subset2$P.Value < 0.05, ]
dea.GSE58144.complex.subset2.up <- dea.GSE58144.complex.subset2[
  dea.GSE58144.complex.subset2$logFC > 0.667 & dea.GSE58144.complex.subset2$P.Value < 0.05, ]
dea.GSE58144.complex.subset2.down <- dea.GSE58144.complex.subset2[
  dea.GSE58144.complex.subset2$logFC < -0.667 & dea.GSE58144.complex.subset2$P.Value < 0.05, ]

dea.GSE58144.simple.subset2.up[!(rownames(dea.GSE58144.simple.subset2.up) %in% rownames(dea.GSE58144.complex.subset2.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in simple mode
dea.GSE58144.complex.subset2.up[!(rownames(dea.GSE58144.complex.subset2.up) %in% rownames(dea.GSE58144.simple.subset2.up)), ] %>% 
  arrange(-logFC) %>% head() # up only in complex mode
dea.GSE58144.simple.subset2.up[rownames(dea.GSE58144.simple.subset2.up) %in% rownames(dea.GSE58144.complex.subset2.up), ] %>% 
  arrange(-logFC) %>% head() # intersection

dea.GSE58144.simple.subset2.down[!(rownames(dea.GSE58144.simple.subset2.down) %in% rownames(dea.GSE58144.complex.subset2.down)), ] %>% 
  arrange(logFC) %>% head() # down only in simple mode
dea.GSE58144.complex.subset2.down[!(rownames(dea.GSE58144.complex.subset2.down) %in% rownames(dea.GSE58144.simple.subset2.down)), ] %>% 
  arrange(logFC) %>% head() # down only in complex mode
dea.GSE58144.simple.subset2.down[rownames(dea.GSE58144.simple.subset2.down) %in% rownames(dea.GSE58144.complex.subset2.down), ] %>% 
  arrange(logFC) %>% head() # intersection

# There is no DEGs detected in simple model and intersections 

# Figure 2N - Venn diagrams for DEGs GSE58144
degs.GSE58144.simple.subset2.up <- rownames(dea.GSE58144.simple.subset2.up[
  !(rownames(dea.GSE58144.simple.subset2.up) %in% rownames(dea.GSE58144.complex.subset2.up)), ])
degs.GSE58144.complex.subset2.up <- rownames(dea.GSE58144.complex.subset2.up[
  !(rownames(dea.GSE58144.complex.subset2.up) %in% rownames(dea.GSE58144.simple.subset2.up)), ])
degs.GSE58144.intersection.subset2.up <- rownames(dea.GSE58144.simple.subset2.up[
  rownames(dea.GSE58144.simple.subset2.up) %in% rownames(dea.GSE58144.complex.subset2.up), ])
degs.GSE58144.simple.subset2.down <- rownames(dea.GSE58144.simple.subset2.down[
  !(rownames(dea.GSE58144.simple.subset2.down) %in% rownames(dea.GSE58144.complex.subset2.down)), ])
degs.GSE58144.complex.subset2.down <- rownames(dea.GSE58144.complex.subset2.down[
  !(rownames(dea.GSE58144.complex.subset2.down) %in% rownames(dea.GSE58144.simple.subset2.down)), ])
degs.GSE58144.intersection.subset2.down <- rownames(dea.GSE58144.simple.subset2.down[
  rownames(dea.GSE58144.simple.subset2.down) %in% rownames(dea.GSE58144.complex.subset2.down), ])

genesets.down <- list('DEA model:\nGene expression ~ 1 + Phenotype + Batch' = 
                        c(degs.GSE58144.simple.subset2.down, degs.GSE58144.intersection.subset2.down),
                      'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst + Batch' = 
                        c(degs.GSE58144.complex.subset2.down, degs.GSE58144.intersection.subset2.down))
genesets.up <- list('DEA model:\nGene expression ~ 1 + Phenotype + Batch' = 
                      c(degs.GSE58144.simple.subset2.up, degs.GSE58144.intersection.subset2.up),
                    'DEA model:\nGene expression ~ 1 + Phenotype + EndEst\n+ Phenotype:EndEst + Batch' = 
                      c(degs.GSE58144.complex.subset2.up, degs.GSE58144.intersection.subset2.up))

p2n1 <- ggVennDiagram(genesets.down,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.6, 0.01))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 28)) +
  scale_y_reverse() +
  annotate('text', label = c('DEA model:\nGene expression ~ 1 + Phenotype + Cohort', 
                             'DEA model:\nGene expression ~ 1 + Phenotype\n+ EndEst + Phenotype:EndEst + Cohort'), 
           x = c(-25, -25), y = c(-1.4, 5.4), size = 4.7)

p2n2 <- ggVennDiagram(genesets.up,
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number\nof genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_fill_gradient(low = '#F4FAFE', high = '#878787', limits = c(0, 28)) +
  scale_y_reverse()

p2n <- wrap_plots(p2n1, p2n2, ncol = 2)

# Figure 2 panel
p2.line1 <- wrap_plots(p2a, p2b, p2c, ncol = 3, widths = c(1.3, 3.7, 3))
p2.line2 <- wrap_plots(p2d, p2e, p2f, ncol = 3, widths = c(1.3, 3.7, 3))
p2.line3 <- wrap_plots(p2g, p2h, p2i, ncol = 3, widths = c(1.3, 3.7, 3))
p2.line4 <- wrap_plots(p2j, p2k, p2l, ncol = 3, widths = c(1.3, 3.7, 3))
p2lines1234 <- wrap_plots(p2.line1, p2.line2, p2.line3, p2.line4, nrow = 4)
ggsave(plot = p2lines1234, filename = 'visualization/Figure 2 a-l.png', width = 14, height = 14, dpi = 300)

p2.line5 <- wrap_plots(p2m, p2n, ncol = 2, widths = c(2, 4.7))
ggsave(plot = p2.line5, filename = 'visualization/Figure 2 mn.png', width = 12, height = 2.5, dpi = 300)






### Figure 3 based on PRJNA1141235 scRNAseq dataset

# Load PRJNA1141235 transcript quantification results and create a Seurat object
samples.files.names <- list.files(path = '/home/pd/datasets/hs_endometrium/sc_rna_seq/PRJNA1141235.RIF/quants/', full.names = TRUE)
samples.files.names

samples.seu.list <- lapply(samples.files.names, function(sample.name) {
  sample.quant.mtx <- tximport(files = paste0(sample.name, '/alevin/quants_mat.gz'), type = 'alevin')
  sample.seu <- CreateSeuratObject(sample.quant.mtx$counts, min.features = 300)
  sample.seu$orig.ident <- gsub('^([^_]*).*', '\\1', gsub('\\..*', '', gsub('.*//', '', sample.name)))
  return(sample.seu)
})

seu.PRJNA1141235 <- merge(x = samples.seu.list[[1]], y = samples.seu.list[-1])
table(seu.PRJNA1141235$orig.ident)

mdata.PRJNA1141235 <- read.csv('/home/pd/datasets/hs_endometrium/sc_rna_seq/PRJNA1141235.RIF/SraRunTable.csv', 
                               header = T, sep = ',', row.names = 1)
head(mdata.PRJNA1141235)

mdata.PRJNA1141235$names <- rownames(mdata.PRJNA1141235)
mdata.PRJNA1141235$GEO.accession <- 'PRJNA1141235'
mdata.PRJNA1141235$GSM.accession <- rownames(mdata.PRJNA1141235)
mdata.PRJNA1141235$Sample.name <- mdata.PRJNA1141235$Sample.Name
mdata.PRJNA1141235$Phenotype <- mdata.PRJNA1141235$disease
mdata.PRJNA1141235$Menstrual.cycle.time <- NA
mdata.PRJNA1141235$Menstrual.cycle.phase <- sub('.*-', '', mdata.PRJNA1141235$Sample.Name)
mdata.PRJNA1141235$Age.years <- mdata.PRJNA1141235$AGE
mdata.PRJNA1141235 <- mdata.PRJNA1141235[, c(39:46)]

mdata.PRJNA1141235 <- mdata.PRJNA1141235 %>%
  mutate(Phenotype = case_when(
    Phenotype == 'unavailable' ~ 'Fertile',
    Phenotype == 'recurrent implantation failure with normal endometrium' ~ 'RIF')) %>% 
  mutate(Menstrual.cycle.phase = case_when(
    Menstrual.cycle.phase == 'LPP' ~ 'Proliferative late',
    Menstrual.cycle.phase == 'MSP' ~ 'Secretory mid'))
table(mdata.PRJNA1141235$Phenotype, mdata.PRJNA1141235$Menstrual.cycle.phase)

seu.PRJNA1141235.metadata <- merge(x = seu.PRJNA1141235@meta.data, y = mdata.PRJNA1141235, 
                                by.x = 'orig.ident', by.y = 'names', all.x = TRUE)
rownames(seu.PRJNA1141235.metadata) <- rownames(seu.PRJNA1141235@meta.data)
seu.PRJNA1141235@meta.data <- seu.PRJNA1141235.metadata

# Filter cells by MT genes expression and doublet annotation
seu.PRJNA1141235 <- JoinLayers(seu.PRJNA1141235)
# create and annotate gse.PRJNA314429 again (line 892) if necessary
features.ids <- data.frame(ensg = rownames(GetAssayData(seu.PRJNA1141235, assay = 'RNA', layer = 'counts')),
                           symbols = mapIds(org.Hs.eg.db,
                                            keys = substr(rownames(GetAssayData(seu.PRJNA1141235, assay = 'RNA', layer = 'counts')), 1, 15),
                                            column = 'SYMBOL', 
                                            keytype = 'ENSEMBL', 
                                            multiVals = 'first'),
                           chr = (as.data.frame(gse.PRJNA314429@rowRanges)[
                             rownames(GetAssayData(seu.PRJNA1141235, assay = 'RNA', layer = 'counts')), ])$seqnames)
seu.PRJNA1141235$percent.mt <- PercentageFeatureSet(seu.PRJNA1141235, features = subset(features.ids, chr == 'chrM')$ensg)
Idents(seu.PRJNA1141235) <- 'orig.ident'
VlnPlot(seu.PRJNA1141235, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, raster = FALSE)

seu.PRJNA1141235 <- subset(seu.PRJNA1141235, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)

samples.seu.list <- SplitObject(seu.PRJNA1141235, split.by = 'orig.ident')
for (i in 1:length(samples.seu.list)) {
  sample.seu <- NormalizeData(samples.seu.list[[i]])
  sample.seu <- FindVariableFeatures(sample.seu, selection.method = 'vst', nfeatures = 2000)
  sample.seu <- ScaleData(sample.seu)
  sample.seu <- RunPCA(sample.seu, dims = 1:20)
  sample.seu <- RunUMAP(sample.seu, dims = 1:20)
  sample.seu <- FindNeighbors(object = sample.seu, dims = 1:20)              
  sample.seu <- FindClusters(object = sample.seu, resolution = 0.2)
  sweep.res.list <- paramSweep(sample.seu, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  annotations <- sample.seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.05 * nrow(sample.seu@meta.data))  ## Assuming 5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  sample.seu <- doubletFinder(sample.seu, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(sample.seu@meta.data)[ncol(sample.seu@meta.data)] <- 'DF_prediction'
  samples.seu.list[[i]] <- sample.seu
  remove(sample.seu)
}

seu.PRJNA1141235.DF <- merge(x = samples.seu.list[[1]], y = samples.seu.list[-1])
seu.PRJNA1141235$DF_prediction <- seu.PRJNA1141235.DF$DF_prediction
seu.PRJNA1141235 <- JoinLayers(seu.PRJNA1141235)
table(seu.PRJNA1141235$DF_prediction, seu.PRJNA1141235$orig.ident)
seu.PRJNA1141235 <- subset(seu.PRJNA1141235, DF_prediction == 'Singlet')
saveRDS(seu.PRJNA1141235, file.path('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.Rds'))


# Prepare separate object for cell type deconvolution in Figure 3 - both phases
seu.PRJNA1141235.fig3 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.Rds')
seu.PRJNA1141235.fig3[['RNA']] <- split(seu.PRJNA1141235.fig3[['RNA']], f = seu.PRJNA1141235.fig3$orig.ident)
seu.PRJNA1141235.fig3 <- NormalizeData(seu.PRJNA1141235.fig3)
seu.PRJNA1141235.fig3 <- FindVariableFeatures(seu.PRJNA1141235.fig3, selection.method = 'vst', nfeatures = 2000)
# Regress out cell cycle variability
seu.PRJNA1141235.fig3 <- CellCycleScoring(seu.PRJNA1141235.fig3,
                                          g2m.features = features.ids[features.ids$symbols %in% cc.genes$g2m.genes, ]$ensg,
                                          s.features = features.ids[features.ids$symbols %in% cc.genes$s.genes, ]$ensg)
seu.PRJNA1141235.fig3 <- ScaleData(seu.PRJNA1141235.fig3, vars.to.regress = c('S.Score', 'G2M.Score'))
seu.PRJNA1141235.fig3 <- RunPCA(seu.PRJNA1141235.fig3, npcs = 20)
# Indicate Menstrual.cycle.phase as an additional variability to harmonise
seu.PRJNA1141235.fig3 <- RunHarmony(seu.PRJNA1141235.fig3, group.by.vars = c('orig.ident', 'Menstrual.cycle.phase'), reduction.save = 'harmony')
seu.PRJNA1141235.fig3 <- RunUMAP(seu.PRJNA1141235.fig3, reduction = 'harmony', dims = 1:20, reduction.name = 'umap.harmony')
seu.PRJNA1141235.fig3 <- FindNeighbors(seu.PRJNA1141235.fig3, reduction = 'harmony', dims = 1:20)
seu.PRJNA1141235.fig3 <- FindClusters(seu.PRJNA1141235.fig3, resolution = 0.2) # less sensitive to get major groups

cells.order.random <- sample(Cells(seu.PRJNA1141235.fig3))
DimPlot(seu.PRJNA1141235.fig3, reduction = 'umap.harmony', group.by = 'orig.ident', 
        label = FALSE, raster=FALSE,cells = cells.order.random)

# Figure 3A - integrated object by phenotype 
p3a <- DimPlot(seu.PRJNA1141235.fig3, reduction = 'umap.harmony', group.by = 'Phenotype', 
        label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom', 
        legend.direction = 'vertical') +
  scale_colour_brewer(palette = 'Set2')

# Figure 3B - integrated object by phase 
p3b <- DimPlot(seu.PRJNA1141235.fig3, reduction = 'umap.harmony', group.by = 'Menstrual.cycle.phase', 
               label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Menstrual cycle phase') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom', 
        legend.direction = 'vertical') +
  scale_color_manual(values = c('Proliferative late' = viridis::viridis(6)[2], 
                                'Secretory mid' = viridis::viridis(6)[4]))

# Get marker genes for cell clusters
seu.PRJNA1141235.fig3 <- JoinLayers(seu.PRJNA1141235.fig3)
all.markers <- FindAllMarkers(seu.PRJNA1141235.fig3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
all.markers$symbol <- mapIds(org.Hs.eg.db, keys = substr(all.markers$gene, 1, 15), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
all.markers.filtered <- all.markers %>%
  filter(!(is.na(symbol))) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = abs(avg_log2FC)) %>%
  arrange(cluster, -abs(avg_log2FC))
summary(as.factor(all.markers$cluster))
summary(as.factor(all.markers.filtered$cluster))

# Annotate marker genes by MSigDB C8 ontology
head(as.data.frame(msigdbr_collections()), 30)
m_t2g <- msigdbr(species = 'Homo sapiens', category = 'C8') %>% dplyr::select(gs_name, gene_symbol)
ck <- compareCluster(symbol ~ cluster,
                     data = all.markers.filtered,
                     fun = 'enricher',
                     TERM2GENE = m_t2g,
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.01)
View(summary(ck))
dotplot(ck, showCategory = 3) + scale_y_discrete(labels=function(x) str_wrap(x, width=80))
DimPlot(seu.PRJNA1141235.fig3, reduction = 'umap.harmony', group.by = 'seurat_clusters', label = TRUE)

# Explore clusters by classic cell type marker genes 
# 'EMCN', 'VWF'     Endothelial
# 'EPCAM', 'CLDN3'  Epithelial
# 'CCL5', 'GZMA'    Lymphoid
# 'AIF1', 'MNDA'    Myeliod
# 'COL1A1', 'DCN'   Stromal
DimPlot(seu.PRJNA1141235.fig3, reduction = 'umap.harmony', group.by = 'seurat_clusters', label = TRUE) +
  FeaturePlot(seu.PRJNA1141235.fig3, feature = subset(features.ids, symbols %in% c('DCN'))$ensg)

# Set cell clusters annotation 
seu.PRJNA1141235.fig3@meta.data <- seu.PRJNA1141235.fig3@meta.data %>% mutate(cell_type = case_when(
  seurat_clusters %in% c('9') ~ 'Myeloid',
  seurat_clusters %in% c('2', '5', '11', '12') ~ 'Lymphoid',
  seurat_clusters %in% c('6') ~ 'Endothelial',
  seurat_clusters %in% c('1', '4', '13') ~ 'Epithelial',
  seurat_clusters %in% c('0', '3', '7', '8', '10') ~ 'Stromal'
))

# Figure 3C - integrated object by cell type
p3c <- DimPlot(seu.PRJNA1141235.fig3, reduction = 'umap.harmony', group.by = 'cell_type', 
               label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom', 
        legend.direction = 'vertical') +
  scale_color_manual(values = c('Myeloid' = brewer.pal(n = 7, name = 'Set2')[3],
                                'Lymphoid' = brewer.pal(n = 7, name = 'Set2')[4],
                                'Endothelial' = brewer.pal(n = 7, name = 'Set2')[5],
                                'Epithelial' = brewer.pal(n = 7, name = 'Set2')[6],
                                'Stromal' = brewer.pal(n = 7, name = 'Set2')[7]))
  
Idents(seu.PRJNA1141235.fig3) <- 'cell_type'
all.markers <- FindAllMarkers(seu.PRJNA1141235.fig3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
all.markers$symbol <- mapIds(org.Hs.eg.db, keys = substr(all.markers$gene, 1, 15), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
all.markers.filtered <- all.markers %>%
  filter(!(is.na(symbol))) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = abs(avg_log2FC)) %>%
  arrange(cluster, -abs(avg_log2FC))

features.ids <- data.frame(ensg = rownames(GetAssayData(seu.PRJNA1141235.fig3, assay = 'RNA', layer = 'counts')),
                           symbols = mapIds(org.Hs.eg.db,
                                            keys = substr(rownames(GetAssayData(seu.PRJNA1141235.fig3, assay = 'RNA', layer = 'counts')), 1, 15),
                                            column = 'SYMBOL', 
                                            keytype = 'ENSEMBL', 
                                            multiVals = 'first'),
                           chr = (as.data.frame(gse.PRJNA314429@rowRanges)[
                             rownames(GetAssayData(seu.PRJNA1141235.fig3, assay = 'RNA', layer = 'counts')), ])$seqnames)

marker.genes <- c(
  'EMCN', 'VWF', # Endothelial
  'EPCAM', 'CLDN3', # Epithelial
  'CCL5', 'GZMA', # Lymphoid
  'AIF1', 'MNDA', # Myeliod
  'COL1A1', 'DCN' # Stromal
)
marker.genes.df <- subset(features.ids, symbols %in% marker.genes)
rownames(marker.genes.df) <- marker.genes.df$symbols
marker.genes.df <- marker.genes.df[marker.genes, ]

p3d <- DotPlot(seu.PRJNA1141235.fig3, assay = 'RNA', features = marker.genes.df$ensg, group.by = 'cell_type', 
        dot.scale = 5, cluster.idents = FALSE) +
  labs(title = '', x = 'Marker genes', y = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'bottom',
        legend.box = 'vertical', 
        legend.margin = margin()) +
  scale_x_discrete(labels = marker.genes) +
  scale_y_discrete(limits = rev(levels(as.factor(seu.PRJNA1141235.fig3$cell_type)))) + 
  scale_color_distiller(type = 'seq', direction = 1, palette = 'Greys', breaks = c(0, 1))
# p3d[['guides']][['guides']][['colour']][['params']][['title']] <- 'Average\nexpression'
# p3d[['guides']][['guides']][['size']][['params']][['title']] <- 'Percent\nexpressed'

gc()
saveRDS(seu.PRJNA1141235.fig3, file.path('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.fig3.Rds'))
write.csv(features.ids, 'processed.data/features.ids.csv')
write.csv(marker.genes.df, 'processed.data/marker.genes.df.csv')
# Use this files to replicate Figure 3 plots
seu.PRJNA1141235.fig3 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.fig3.Rds')
features.ids <- read.csv('processed.data/features.ids.csv')
marker.genes.df <- read.csv('processed.data/marker.genes.df.csv')


# Prepare reference for deconvolution
seu.PRJNA1141235.fig3.reference <- AverageExpression(seu.PRJNA1141235.fig3, group.by = c('cell_type'), 
                                                     assays = 'RNA', slot = 'counts', return.seurat = FALSE)$RNA
seu.PRJNA1141235.fig3.reference <- as.data.frame(seu.PRJNA1141235.fig3.reference)
features.ids.dcnv.reference <- features.ids
rownames(features.ids.dcnv.reference) <- features.ids.dcnv.reference$ensg
features.ids.dcnv.reference <- features.ids.dcnv.reference[rownames(features.ids.dcnv.reference) %in% rownames(seu.PRJNA1141235.fig3.reference), ]
features.ids.dcnv.reference <- features.ids.dcnv.reference[rownames(seu.PRJNA1141235.fig3.reference), ]
seu.PRJNA1141235.fig3.reference$symbols <- features.ids.dcnv.reference$symbols
seu.PRJNA1141235.fig3.reference <- subset(seu.PRJNA1141235.fig3.reference, !(is.na(symbols)) & !(duplicated(symbols)))
rownames(seu.PRJNA1141235.fig3.reference) <- seu.PRJNA1141235.fig3.reference$symbols
head(seu.PRJNA1141235.fig3.reference)
# slice reference mtx to remove low expressed genes with high variance
seu.PRJNA1141235.fig3.reference$sum.exp <- rowSums(seu.PRJNA1141235.fig3.reference[, c(1:5)])
seu.PRJNA1141235.fig3.reference.sliced <- seu.PRJNA1141235.fig3.reference %>% slice_max(sum.exp, n = 10000)
seu.PRJNA1141235.fig3.reference.sliced <- seu.PRJNA1141235.fig3.reference.sliced[, c(1:5)]
write.csv(seu.PRJNA1141235.fig3.reference.sliced, 'processed.data/seu.PRJNA1141235.fig3.reference.sliced.csv')



# Deconvolute bulk datasets for Figure 3

# clone GitHub repository https://github.com/BNadel/GEDITExpanded
# use python3 and run GEDITExpanded-main/GEDITv3.0/GEDIT3.py script in terminal as follows
# python3 GEDIT3.py -mix /home/pd/projects/article_2025_RIF/project_RIF/processed.data/exp.GSE111974.max.csv -ref /home/pd/projects/article_2025_RIF/project_RIF/processed.data/seu.PRJNA1141235.fig3.reference.sliced.csv -outFile /home/pd/projects/article_2025_RIF/project_RIF/processed.data/GSE111974.dcnv.results
# python3 GEDIT3.py -mix /home/pd/projects/article_2025_RIF/project_RIF/processed.data/exp.GSE58144.max.subset.csv -ref /home/pd/projects/article_2025_RIF/project_RIF/processed.data/seu.PRJNA1141235.fig3.reference.sliced.csv -outFile /home/pd/projects/article_2025_RIF/project_RIF/processed.data/GSE58144.subset.dcnv.results
# python3 GEDIT3.py -mix /home/pd/projects/article_2025_RIF/project_RIF/processed.data/vst.GSE207362.mtx.csv -ref /home/pd/projects/article_2025_RIF/project_RIF/processed.data/seu.PRJNA1141235.fig3.reference.sliced.csv -outFile /home/pd/projects/article_2025_RIF/project_RIF/processed.data/GSE207362.dcnv.results
# python3 GEDIT3.py -mix /home/pd/projects/article_2025_RIF/project_RIF/processed.data/vst.GSE243550.mtx.csv -ref /home/pd/projects/article_2025_RIF/project_RIF/processed.data/seu.PRJNA1141235.fig3.reference.sliced.csv -outFile /home/pd/projects/article_2025_RIF/project_RIF/processed.data/GSE243550.dcnv.results

gedit_result.GSE111974 <- read.csv('processed.data/GSE111974.dcnv.results_CTPredictions.tsv', sep = '\t', header = T, row.names = 1)
gedit_result.GSE207362 <- read.csv('processed.data/GSE207362.dcnv.results_CTPredictions.tsv', sep = '\t', header = T, row.names = 1)
gedit_result.GSE243550 <- read.csv('processed.data/GSE243550.dcnv.results_CTPredictions.tsv', sep = '\t', header = T, row.names = 1)
gedit_result.GSE58144.subset <- read.csv('processed.data/GSE58144.subset.dcnv.results_CTPredictions.tsv', sep = '\t', header = T, row.names = 1)

mdata.GSE111974 <- read.csv('processed.data/mdata.GSE111974.csv', sep = ',', header = T, row.names = 1)
mdata.GSE207362 <- read.csv('processed.data/mdata.GSE207362.csv', sep = ',', header = T, row.names = 1)
mdata.GSE243550 <- read.csv('processed.data/mdata.GSE243550.csv', sep = ',', header = T, row.names = 1)
mdata.GSE58144.subset <- read.csv('processed.data/mdata.GSE58144.subset.csv', sep = ',', header = T, row.names = 1)

mdata.GSE111974 <- cbind(mdata.GSE111974, gedit_result.GSE111974)
mdata.GSE207362 <- cbind(mdata.GSE207362, gedit_result.GSE207362)
mdata.GSE243550 <- cbind(mdata.GSE243550, gedit_result.GSE243550)
mdata.GSE58144.subset <- cbind(mdata.GSE58144.subset, gedit_result.GSE58144.subset)

mdata.GSE111974.subset <- subset(mdata.GSE111974, EndEst > 31 & EndEst < 85)
mdata.GSE207362.subset <- subset(mdata.GSE207362, EndEst > 72 & EndEst < 83)
mdata.GSE243550.subset <- subset(mdata.GSE243550, EndEst > 66 & EndEst < 80)
mdata.GSE58144.subset2 <- subset(mdata.GSE58144.subset, EndEst > 58 & EndEst < 76)

mdata.merged.dcnv.plot.GSE111974.subset <- mdata.GSE111974.subset %>% 
  pivot_longer(cols = colnames(mdata.GSE111974.subset)[c(12:16)], names_to = 'cell_type', values_to = 'Proportion')

p3e1 <- ggplot(mdata.merged.dcnv.plot.GSE111974.subset, aes(EndEst, Proportion * 100, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  facet_rep_grid(GEO.accession ~ cell_type) +
  # facet_wrap(GEO.accession ~ cell_type, scales = 'free', ncol = 5) +
  labs(title = '', x = 'EndEst', y = 'Cell type percent') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 13),
        strip.text.x = element_text(size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'none',
        panel.border = element_blank(), 
        axis.line = element_line()) +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

mdata.merged.dcnv.plot.GSE207362.subset <- mdata.GSE207362.subset %>% 
  pivot_longer(cols = colnames(mdata.GSE207362.subset)[c(13:17)], names_to = 'cell_type', values_to = 'Proportion')

p3e2 <- ggplot(mdata.merged.dcnv.plot.GSE207362.subset, aes(EndEst, Proportion * 100, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  facet_rep_grid(GEO.accession ~ cell_type) +
  # facet_wrap(GEO.accession ~ cell_type, scales = 'free', ncol = 5) +
  labs(title = '', x = 'EndEst', y = 'Cell type percent') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 13),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'none',
        panel.border = element_blank(), 
        axis.line = element_line()) +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

mdata.merged.dcnv.plot.GSE243550.subset <- mdata.GSE243550.subset %>% 
  pivot_longer(cols = colnames(mdata.GSE243550.subset)[c(13:17)], names_to = 'cell_type', values_to = 'Proportion')

p3e3 <- ggplot(mdata.merged.dcnv.plot.GSE243550.subset, aes(EndEst, Proportion * 100, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  facet_rep_grid(GEO.accession ~ cell_type) +
  # facet_wrap(GEO.accession ~ cell_type, scales = 'free', ncol = 5) +
  labs(title = '', x = 'EndEst', y = 'Cell type percent') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 13),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'none',
        panel.border = element_blank(), 
        axis.line = element_line()) +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

mdata.merged.dcnv.plot.GSE58144.subset2 <- mdata.GSE58144.subset2 %>% 
  pivot_longer(cols = colnames(mdata.GSE58144.subset)[c(13:17)], names_to = 'cell_type', values_to = 'Proportion')

head(mdata.merged.dcnv.plot.GSE58144.subset2)

mdata.merged.dcnv.plot.GSE58144.subset2 <- mdata.merged.dcnv.plot.GSE58144.subset2 %>% 
  mutate(GEO.accession = case_when(
    Batch == 'Cohort 1' ~ 'GSE58144\nCohort 1',
    Batch == 'Cohort 2, batch 1' ~ 'GSE58144\nCohort 2, batch 1'))

p3e4 <- ggplot(subset(mdata.merged.dcnv.plot.GSE58144.subset2, Batch == 'Cohort 1'), 
               aes(EndEst, Proportion * 100, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  facet_rep_grid(GEO.accession ~ cell_type) +
  # facet_wrap(GEO.accession ~ cell_type, scales = 'free', ncol = 5) +
  labs(title = '', x = 'EndEst', y = 'Cell type percent') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 13),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'none',
        panel.border = element_blank(), 
        axis.line = element_line()) +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

p3e5 <- ggplot(subset(mdata.merged.dcnv.plot.GSE58144.subset2, Batch == 'Cohort 2, batch 1'), 
               aes(EndEst, Proportion * 100, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  facet_rep_grid(GEO.accession ~ cell_type) +
  # facet_wrap(GEO.accession ~ cell_type, scales = 'free', ncol = 5) +
  labs(title = '', x = 'EndEst', y = 'Cell type percent') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        strip.background = element_blank(),
        strip.text.y = element_text(size = 13),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom',
        panel.border = element_blank(), 
        axis.line = element_line()) +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

# Figure 3 panel
X <- ns(mdata.GSE111974.subset$EndEst, df = 3)
Group <- factor(mdata.GSE111974.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(t(mdata.GSE111974.subset[, c(12:16)]), design)
fit <- eBayes(fit)
dea.dcnv.GSE111974.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.dcnv.GSE111974.complex.subset

X <- ns(mdata.GSE207362.subset$EndEst, df = 3)
Group <- factor(mdata.GSE207362.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(t(mdata.GSE207362.subset[, c(13:17)]), design)
fit <- eBayes(fit)
dea.dcnv.GSE207362.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.dcnv.GSE207362.complex.subset

X <- ns(mdata.GSE243550.subset$EndEst, df = 3)
Group <- factor(mdata.GSE243550.subset$Phenotype)
design <- model.matrix(~ Group * X)
fit <- lmFit(t(mdata.GSE243550.subset[, c(13:17)]), design)
fit <- eBayes(fit)
dea.dcnv.GSE243550.complex.subset <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.dcnv.GSE243550.complex.subset

X <- ns(mdata.GSE58144.subset2$EndEst, df = 3)
Group <- factor(mdata.GSE58144.subset2$Phenotype)
Batch <- factor(mdata.GSE58144.subset2$Batch)
design <- model.matrix(~ Group * X + Batch)
fit <- lmFit(t(mdata.GSE58144.subset2[, c(13:17)]), design)
fit <- eBayes(fit)
dea.dcnv.GSE58144.complex.subset2 <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.dcnv.GSE58144.complex.subset2

p3.line1 <- wrap_plots(p3a, p3b, p3c, p3d, ncol = 4, widths = c(1, 1, 1, 1.5))
p3.line2 <- wrap_plots(p3e1, p3e2, p3e3, p3e4, p3e5, nrow = 5, heights = c(1, 1, 1, 1, 1))
p3 <- wrap_plots(p3.line1, p3.line2, nrow = 2, heights = c(1, 6))
ggsave(plot = p3, filename = 'visualization/Figure 3.png', width = 12, height = 15, dpi = 300)







### Figure 4 & S3 - scRNAseq multi-dataset analysis

features.ids <- read.csv('processed.data/features.ids.csv', sep = ',', header = T, row.names = 1)
marker.genes.df <- read.csv('processed.data/marker.genes.df.csv', sep = ',', header = T, row.names = 1)
marker.genes <- c(
  'EMCN', 'VWF', # Endothelial
  'EPCAM', 'CLDN3', # Epithelial
  'CCL5', 'GZMA', # Lymphoid
  'AIF1', 'MNDA', # Myeliod
  'COL1A1', 'DCN' # Stromal
)

## PRJNA1141235 subset of secretory phase samples
seu.PRJNA1141235.fig3 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.fig3.Rds')
table(seu.PRJNA1141235.fig3$Menstrual.cycle.phase, seu.PRJNA1141235.fig3$Phenotype)
seu.PRJNA1141235.fig4 <- subset(seu.PRJNA1141235.fig3, Menstrual.cycle.phase == 'Secretory mid')
# Rerun projection on a reduced number of objects, harmonisation only by donor
seu.PRJNA1141235.fig4 <- RunPCA(seu.PRJNA1141235.fig4, npcs = 20)
seu.PRJNA1141235.fig4 <- RunHarmony(seu.PRJNA1141235.fig4, group.by.vars = c('orig.ident'), reduction.save = 'harmony')
seu.PRJNA1141235.fig4 <- RunUMAP(seu.PRJNA1141235.fig4, reduction = 'harmony', dims = 1:20, reduction.name = 'umap.harmony')
seu.PRJNA1141235.fig4 <- FindNeighbors(seu.PRJNA1141235.fig4, reduction = 'harmony', dims = 1:20)
seu.PRJNA1141235.fig4 <- FindClusters(seu.PRJNA1141235.fig4, resolution = 0.2)

cells.order.random <- sample(Cells(seu.PRJNA1141235.fig4))
ps3a <- DimPlot(seu.PRJNA1141235.fig4, reduction = 'umap.harmony', group.by = 'Phenotype', 
               label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right') +
  scale_colour_brewer(palette = 'Set2')

ps3b <- DimPlot(seu.PRJNA1141235.fig4, reduction = 'umap.harmony', group.by = 'cell_type', 
               label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right') +
  scale_color_manual(values = c('Myeloid' = brewer.pal(n = 7, name = 'Set2')[3],
                                'Lymphoid' = brewer.pal(n = 7, name = 'Set2')[4],
                                'Endothelial' = brewer.pal(n = 7, name = 'Set2')[5],
                                'Epithelial' = brewer.pal(n = 7, name = 'Set2')[6],
                                'Stromal' = brewer.pal(n = 7, name = 'Set2')[7]))

seu.PRJNA1141235.fig4 <- JoinLayers(seu.PRJNA1141235.fig4)
Idents(seu.PRJNA1141235.fig4) <- 'cell_type'
ps3c <- DotPlot(seu.PRJNA1141235.fig4, assay = 'RNA', features = marker.genes.df$ensg, group.by = 'cell_type', 
               dot.scale = 5, cluster.idents = FALSE) +
  labs(title = '', x = 'Marker genes', y = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right',
        legend.box = 'horizontal',
        legend.margin = margin()) +
  scale_x_discrete(labels = marker.genes) +
  scale_y_discrete(limits = rev(levels(as.factor(seu.PRJNA1141235.fig4$cell_type)))) + 
  scale_color_distiller(type = 'seq', direction = 1, palette = 'Greys')
ps3c[['guides']][['guides']][['colour']][['params']][['title']] <- 'Average\nexpression'
ps3c[['guides']][['guides']][['size']][['params']][['title']] <- 'Percent\nexpressed'

gc()
saveRDS(seu.PRJNA1141235.fig4, file.path('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.fig4.Rds'))
seu.PRJNA1141235.fig4 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.fig4.Rds')



## GSE183837

samples.files.names <- list.files(path = '/home/pd/datasets/hs_endometrium/sc_rna_seq/GSE183837.RIF/quants/', full.names = TRUE)
samples.files.names

samples.seu.list <- lapply(samples.files.names, function(sample.name) {
  sample.quant.mtx <- tximport(files = paste0(sample.name, '/alevin/quants_mat.gz'), type = 'alevin')
  sample.seu <- CreateSeuratObject(sample.quant.mtx$counts, min.features = 300)
  sample.seu$orig.ident <- gsub('^([^_]*).*', '\\1', gsub('\\..*', '', gsub('.*//', '', sample.name)))
  return(sample.seu)
})

seu.GSE183837 <- merge(x = samples.seu.list[[1]], y = samples.seu.list[-1])
table(seu.GSE183837$orig.ident)

mdata.GSE183837 <- read.csv('/home/pd/datasets/hs_endometrium/sc_rna_seq/GSE183837.RIF/SraRunTable.csv', 
                            header = T, sep = ',', row.names = 1)
head(mdata.GSE183837)
table(mdata.GSE183837$Group)

mdata.GSE183837$names <- rownames(mdata.GSE183837)
mdata.GSE183837$GEO.accession <- 'GSE183837'
mdata.GSE183837$GSM.accession <- mdata.GSE183837$GEO_Accession..exp.
mdata.GSE183837$Sample.name <- mdata.GSE183837$Sample.Name
mdata.GSE183837$Phenotype <- mdata.GSE183837$Group
mdata.GSE183837$Menstrual.cycle.time <- mdata.GSE183837$STAGE
mdata.GSE183837$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE183837$Age.years <- NA
mdata.GSE183837 <- mdata.GSE183837[, c(29:36)]

mdata.GSE183837 <- mdata.GSE183837 %>% 
  mutate(Phenotype = case_when(
    Phenotype == 'control' ~ 'Fertile',
    Phenotype == 'RIF' ~ 'RIF'))
table(mdata.GSE183837$Phenotype)

head(seu.GSE183837@meta.data)
head(mdata.GSE183837)
seu.GSE183837.metadata <- merge(x = seu.GSE183837@meta.data, y = mdata.GSE183837, by.x = 'orig.ident', by.y = 'names', all.x = TRUE)
rownames(seu.GSE183837.metadata) <- rownames(seu.GSE183837@meta.data)
seu.GSE183837@meta.data <- seu.GSE183837.metadata

seu.GSE183837 <- JoinLayers(seu.GSE183837)
seu.GSE183837$percent.mt <- PercentageFeatureSet(seu.GSE183837, features = subset(features.ids, chr == 'chrM')$ensg)
Idents(seu.GSE183837) <- 'orig.ident'
VlnPlot(seu.GSE183837, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, raster = FALSE)

seu.GSE183837 <- subset(seu.GSE183837, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)

samples.seu.list <- SplitObject(seu.GSE183837, split.by = 'orig.ident')
for (i in 1:length(samples.seu.list)) {
  sample.seu <- NormalizeData(samples.seu.list[[i]])
  sample.seu <- FindVariableFeatures(sample.seu, selection.method = 'vst', nfeatures = 2000)
  sample.seu <- ScaleData(sample.seu)
  sample.seu <- RunPCA(sample.seu, dims = 1:20)
  sample.seu <- RunUMAP(sample.seu, dims = 1:20)
  sample.seu <- FindNeighbors(object = sample.seu, dims = 1:20)              
  sample.seu <- FindClusters(object = sample.seu, resolution = 0.2)
  sweep.res.list <- paramSweep(sample.seu, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  annotations <- sample.seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.05 * nrow(sample.seu@meta.data))  ## Assuming 5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  sample.seu <- doubletFinder(sample.seu, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(sample.seu@meta.data)[ncol(sample.seu@meta.data)] <- 'DF_prediction'
  samples.seu.list[[i]] <- sample.seu
  remove(sample.seu)
}

seu.GSE183837.DF <- merge(x = samples.seu.list[[1]], y = samples.seu.list[-1])
seu.GSE183837$DF_prediction <- seu.GSE183837.DF$DF_prediction
seu.GSE183837 <- JoinLayers(seu.GSE183837)
table(seu.GSE183837$DF_prediction, seu.GSE183837$orig.ident)
seu.GSE183837 <- subset(seu.GSE183837, DF_prediction == 'Singlet')

rm(list = setdiff(ls(), c('seu.GSE183837', 'features.ids')))
gc()

seu.GSE183837[['RNA']] <- split(seu.GSE183837[['RNA']], f = seu.GSE183837$orig.ident)
seu.GSE183837 <- NormalizeData(seu.GSE183837)
seu.GSE183837 <- FindVariableFeatures(seu.GSE183837, selection.method = 'vst', nfeatures = 2000)
seu.GSE183837 <- CellCycleScoring(seu.GSE183837,
                                  g2m.features = features.ids[features.ids$symbols %in% cc.genes$g2m.genes, ]$ensg,
                                  s.features = features.ids[features.ids$symbols %in% cc.genes$s.genes, ]$ensg)
seu.GSE183837 <- ScaleData(seu.GSE183837, vars.to.regress = c('S.Score', 'G2M.Score'))
seu.GSE183837 <- RunPCA(seu.GSE183837, npcs = 20)
seu.GSE183837 <- RunHarmony(seu.GSE183837, group.by.vars = c('orig.ident'), reduction.save = 'harmony')
seu.GSE183837 <- RunUMAP(seu.GSE183837, reduction = 'harmony', dims = 1:20, reduction.name = 'umap.harmony')
seu.GSE183837 <- FindNeighbors(seu.GSE183837, reduction = 'harmony', dims = 1:20)
seu.GSE183837 <- FindClusters(seu.GSE183837, resolution = 0.2)

cells.order.random <- sample(Cells(seu.GSE183837))
ps3d <- DimPlot(seu.GSE183837, reduction = 'umap.harmony', group.by = 'Phenotype', 
                label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right') +
  scale_colour_brewer(palette = 'Set2')

seu.GSE183837 <- JoinLayers(seu.GSE183837)
all.markers <- FindAllMarkers(seu.GSE183837, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
all.markers$symbol <- mapIds(org.Hs.eg.db, keys = substr(all.markers$gene, 1, 15), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
all.markers.filtered <- all.markers %>%
  filter(!(is.na(symbol))) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = abs(avg_log2FC)) %>%
  arrange(cluster, -abs(avg_log2FC))
summary(as.factor(all.markers$cluster))
summary(as.factor(all.markers.filtered$cluster))

# head(as.data.frame(msigdbr_collections()), 30)
m_t2g <- msigdbr(species = 'Homo sapiens', category = 'C8') %>% dplyr::select(gs_name, gene_symbol)
ck <- compareCluster(symbol ~ cluster,
                     data = all.markers.filtered,
                     fun = 'enricher',
                     TERM2GENE = m_t2g,
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.01)
View(summary(ck))
dotplot(ck, showCategory = 3) + scale_y_discrete(labels=function(x) str_wrap(x, width=80))
DimPlot(seu.GSE183837, reduction = 'umap.harmony', group.by = 'seurat_clusters', label = TRUE)
marker.genes
FeaturePlot(seu.GSE183837, feature = subset(features.ids, symbols %in% c('EPCAM'))$ensg)

seu.GSE183837@meta.data <- seu.GSE183837@meta.data %>% mutate(cell_type = case_when(
  seurat_clusters %in% c('6') ~ 'Myeloid',
  seurat_clusters %in% c('3', '4') ~ 'Lymphoid',
  seurat_clusters %in% c('9') ~ 'Endothelial',
  seurat_clusters %in% c('7') ~ 'Epithelial',
  seurat_clusters %in% c('0', '1', '2', '5', '8') ~ 'Stromal'
))

ps3e <- DimPlot(seu.GSE183837, reduction = 'umap.harmony', group.by = 'cell_type', 
                label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right') +
  scale_color_manual(values = c('Myeloid' = brewer.pal(n = 7, name = 'Set2')[3],
                                'Lymphoid' = brewer.pal(n = 7, name = 'Set2')[4],
                                'Endothelial' = brewer.pal(n = 7, name = 'Set2')[5],
                                'Epithelial' = brewer.pal(n = 7, name = 'Set2')[6],
                                'Stromal' = brewer.pal(n = 7, name = 'Set2')[7]))

Idents(seu.GSE183837) <- 'cell_type'

ps3f <- DotPlot(seu.GSE183837, assay = 'RNA', features = marker.genes.df$ensg, group.by = 'cell_type', 
                dot.scale = 5, cluster.idents = FALSE) +
  labs(title = '', x = 'Marker genes', y = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right',
        legend.box = 'horizontal',
        legend.margin = margin()) +
  scale_x_discrete(labels = marker.genes) +
  scale_y_discrete(limits = rev(levels(as.factor(seu.GSE183837$cell_type)))) + 
  scale_color_distiller(type = 'seq', direction = 1, palette = 'Greys')
ps3f[['guides']][['guides']][['colour']][['params']][['title']] <- 'Average\nexpression'
ps3f[['guides']][['guides']][['size']][['params']][['title']] <- 'Percent\nexpressed'

gc()
saveRDS(seu.GSE183837, file.path('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.GSE183837.Rds'))
seu.GSE183837 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.GSE183837.Rds')



## GSE250130

samples.files.names <- list.files(path = '/home/pd/datasets/hs_endometrium/sc_rna_seq/GSE250130.RIF/quants/', full.names = TRUE)
samples.files.names

samples.seu.list <- lapply(samples.files.names, function(sample.name) {
  sample.quant.mtx <- tximport(files = paste0(sample.name, '/alevin/quants_mat.gz'), type = 'alevin')
  sample.seu <- CreateSeuratObject(sample.quant.mtx$counts, min.features = 300)
  sample.seu$orig.ident <- gsub('^([^_]*).*', '\\1', gsub('\\..*', '', gsub('.*//', '', sample.name)))
  return(sample.seu)
})

seu.GSE250130 <- merge(x = samples.seu.list[[1]], y = samples.seu.list[-1])
table(seu.GSE250130$orig.ident)

mdata.GSE250130 <- read.csv('/home/pd/datasets/hs_endometrium/sc_rna_seq/GSE250130.RIF/SraRunTable.csv', 
                            header = T, sep = ',', row.names = 1)
head(mdata.GSE250130)

mdata.GSE250130$names <- rownames(mdata.GSE250130)
mdata.GSE250130$GEO.accession <- 'GSE250130'
mdata.GSE250130$GSM.accession <- rownames(mdata.GSE250130)
mdata.GSE250130$Sample.name <- mdata.GSE250130$Sample.Name
mdata.GSE250130$Phenotype <- sub('.* ', '', sub('_.*', '', mdata.GSE250130$Library.Name))
mdata.GSE250130$Menstrual.cycle.time <- 'LH+7'
mdata.GSE250130$Menstrual.cycle.phase <- 'Secretory mid'
mdata.GSE250130$Age.years <- mdata.GSE250130$AGE
mdata.GSE250130 <- mdata.GSE250130[, c(35:42)]

mdata.GSE250130 <- mdata.GSE250130 %>% 
  mutate(Phenotype = case_when(
    Phenotype == 'LH7' ~ 'Fertile',
    Phenotype == 'RIF' ~ 'RIF'))
table(mdata.GSE250130$Phenotype)

head(seu.GSE250130@meta.data)
head(mdata.GSE250130)
seu.GSE250130.metadata <- merge(x = seu.GSE250130@meta.data, y = mdata.GSE250130, by.x = 'orig.ident', by.y = 'names', all.x = TRUE)
rownames(seu.GSE250130.metadata) <- rownames(seu.GSE250130@meta.data)
seu.GSE250130@meta.data <- seu.GSE250130.metadata

seu.GSE250130 <- JoinLayers(seu.GSE250130)
seu.GSE250130$percent.mt <- PercentageFeatureSet(seu.GSE250130, features = subset(features.ids, chr == 'chrM')$ensg)
Idents(seu.GSE250130) <- 'orig.ident'
VlnPlot(seu.GSE250130, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3, raster = FALSE)

seu.GSE250130 <- subset(seu.GSE250130, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 20)

samples.seu.list <- SplitObject(seu.GSE250130, split.by = 'orig.ident')
for (i in 1:length(samples.seu.list)) {
  sample.seu <- NormalizeData(samples.seu.list[[i]])
  sample.seu <- FindVariableFeatures(sample.seu, selection.method = 'vst', nfeatures = 2000)
  sample.seu <- ScaleData(sample.seu)
  sample.seu <- RunPCA(sample.seu, dims = 1:20)
  sample.seu <- RunUMAP(sample.seu, dims = 1:20)
  sample.seu <- FindNeighbors(object = sample.seu, dims = 1:20)              
  sample.seu <- FindClusters(object = sample.seu, resolution = 0.2)
  sweep.res.list <- paramSweep(sample.seu, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn %>% filter(BCmetric == max(BCmetric)) %>% dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  annotations <- sample.seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.05 * nrow(sample.seu@meta.data))  ## Assuming 5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  sample.seu <- doubletFinder(sample.seu, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  colnames(sample.seu@meta.data)[ncol(sample.seu@meta.data)] <- 'DF_prediction'
  samples.seu.list[[i]] <- sample.seu
  remove(sample.seu)
}

seu.GSE250130.DF <- merge(x = samples.seu.list[[1]], y = samples.seu.list[-1])
seu.GSE250130$DF_prediction <- seu.GSE250130.DF$DF_prediction
seu.GSE250130 <- JoinLayers(seu.GSE250130)
table(seu.GSE250130$DF_prediction, seu.GSE250130$orig.ident)
seu.GSE250130 <- subset(seu.GSE250130, DF_prediction == 'Singlet')

seu.GSE250130[['RNA']] <- split(seu.GSE250130[['RNA']], f = seu.GSE250130$orig.ident)
seu.GSE250130 <- NormalizeData(seu.GSE250130)
seu.GSE250130 <- FindVariableFeatures(seu.GSE250130, selection.method = 'vst', nfeatures = 2000)
seu.GSE250130 <- CellCycleScoring(seu.GSE250130,
                                  g2m.features = features.ids[features.ids$symbols %in% cc.genes$g2m.genes, ]$ensg,
                                  s.features = features.ids[features.ids$symbols %in% cc.genes$s.genes, ]$ensg)
seu.GSE250130 <- ScaleData(seu.GSE250130, vars.to.regress = c('S.Score', 'G2M.Score'))
seu.GSE250130 <- RunPCA(seu.GSE250130, npcs = 20)
seu.GSE250130 <- RunHarmony(seu.GSE250130, group.by.vars = c('orig.ident'), reduction.save = 'harmony')
seu.GSE250130 <- RunUMAP(seu.GSE250130, reduction = 'harmony', dims = 1:20, reduction.name = 'umap.harmony')
seu.GSE250130 <- FindNeighbors(seu.GSE250130, reduction = 'harmony', dims = 1:20)
seu.GSE250130 <- FindClusters(seu.GSE250130, resolution = 0.2)

cells.order.random <- sample(Cells(seu.GSE250130))
ps3g <- DimPlot(seu.GSE250130, reduction = 'umap.harmony', group.by = 'Phenotype', 
                 label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right') +
  scale_colour_brewer(palette = 'Set2')

seu.GSE250130 <- JoinLayers(seu.GSE250130)
all.markers <- FindAllMarkers(seu.GSE250130, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, return.thresh = 0.05)
all.markers$symbol <- mapIds(org.Hs.eg.db, keys = substr(all.markers$gene, 1, 15), column = 'SYMBOL', keytype = 'ENSEMBL', multiVals = 'first')
all.markers.filtered <- all.markers %>%
  filter(!(is.na(symbol))) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = abs(avg_log2FC)) %>%
  arrange(cluster, -abs(avg_log2FC))
summary(as.factor(all.markers$cluster))
summary(as.factor(all.markers.filtered$cluster))

# head(as.data.frame(msigdbr_collections()), 30)
m_t2g <- msigdbr(species = 'Homo sapiens', category = 'C8') %>% dplyr::select(gs_name, gene_symbol)
ck <- compareCluster(symbol ~ cluster,
                     data = all.markers.filtered,
                     fun = 'enricher',
                     TERM2GENE = m_t2g,
                     pAdjustMethod = 'BH',
                     pvalueCutoff = 0.01)
View(summary(ck))
dotplot(ck, showCategory = 3) + scale_y_discrete(labels=function(x) str_wrap(x, width=80))
DimPlot(seu.GSE250130, reduction = 'umap.harmony', group.by = 'seurat_clusters', label = TRUE)
marker.genes
FeaturePlot(seu.GSE250130, feature = subset(features.ids, symbols %in% c('GZMA'))$ensg)

seu.GSE250130@meta.data <- seu.GSE250130@meta.data %>% mutate(cell_type = case_when(
  seurat_clusters %in% c('4') ~ 'Myeloid',
  seurat_clusters %in% c('1', '3', '5', '6', '8', '9') ~ 'Lymphoid',
  seurat_clusters %in% c('12') ~ 'Endothelial',
  seurat_clusters %in% c('2', '7', '11', '10') ~ 'Epithelial',
  seurat_clusters %in% c('0') ~ 'Stromal'
))

ps3h <- DimPlot(seu.GSE250130, reduction = 'umap.harmony', group.by = 'cell_type', 
                label = FALSE, raster = FALSE, cells = cells.order.random) +
  labs(title = '', x = 'UMAP Harmony 2', y = 'UMAP Harmony 1', color = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right') +
  scale_color_manual(values = c('Myeloid' = brewer.pal(n = 7, name = 'Set2')[3],
                                'Lymphoid' = brewer.pal(n = 7, name = 'Set2')[4],
                                'Endothelial' = brewer.pal(n = 7, name = 'Set2')[5],
                                'Epithelial' = brewer.pal(n = 7, name = 'Set2')[6],
                                'Stromal' = brewer.pal(n = 7, name = 'Set2')[7]))

Idents(seu.GSE250130) <- 'cell_type'

ps3i <- DotPlot(seu.GSE250130, assay = 'RNA', features = marker.genes.df$ensg, group.by = 'cell_type', 
                dot.scale = 5, cluster.idents = FALSE) +
  labs(title = '', x = 'Marker genes', y = 'Cell type') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = 'right', 
        legend.box = 'horizontal',
        legend.margin = margin()) +
  scale_x_discrete(labels = marker.genes) +
  scale_y_discrete(limits = rev(levels(as.factor(seu.GSE250130$cell_type)))) + 
  scale_color_distiller(type = 'seq', direction = 1, palette = 'Greys')
ps3i[['guides']][['guides']][['colour']][['params']][['title']] <- 'Average\nexpression'
ps3i[['guides']][['guides']][['size']][['params']][['title']] <- 'Percent\nexpressed'

gc()
saveRDS(seu.GSE250130, file.path('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.GSE250130.Rds'))
seu.GSE250130 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.GSE250130.Rds')


## Figure S3
ps3 <- wrap_plots(ps3a, ps3b, ps3c, ps3d, ps3e, ps3f, ps3g, ps3h, ps3i, ncol = 3,
                  widths = c(0.9, 0.9, 1, 0.9, 0.9, 1, 0.9, 0.9, 1))
ggsave(plot = ps3, filename = 'visualization/Figure S3.png', width = 14, height = 9, dpi = 300)





## Figure 4A - Whole-sample pseudobulks EndEst annotation
seu.PRJNA1141235.fig4 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.PRJNA1141235.fig4.Rds')
seu.GSE250130 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.GSE250130.Rds')
seu.GSE183837 <- readRDS('/home/pd/projects/article_2025_RIF/project_RIF_Zenodo/seu.GSE183837.Rds')
features.ids <- read.csv('processed.data/features.ids.csv', sep = ',', header = T, row.names = 1)
marker.genes.df <- read.csv('processed.data/marker.genes.df.csv', sep = ',', header = T, row.names = 1)
rownames(features.ids) <- features.ids$ensg

# PRJNA1141235 WOI subset
head(seu.PRJNA1141235.fig4@meta.data)
pseudo.PRJNA1141235.mdata <- seu.PRJNA1141235.fig4@meta.data[ , c(4:10)] %>% distinct()
rownames(pseudo.PRJNA1141235.mdata) <- pseudo.PRJNA1141235.mdata$GSM.accession
pseudo.PRJNA1141235.exp <- AggregateExpression(seu.PRJNA1141235.fig4, group.by = c('orig.ident'), return.seurat = FALSE)$RNA
pseudo.PRJNA1141235.exp <- as.data.frame(pseudo.PRJNA1141235.exp)
features.ids.PRJNA1141235 <- features.ids[rownames(features.ids) %in% rownames(pseudo.PRJNA1141235.exp), ]
features.ids.PRJNA1141235 <- features.ids.PRJNA1141235[rownames(pseudo.PRJNA1141235.exp), ]
pseudo.PRJNA1141235.exp$symbols <- features.ids.PRJNA1141235$symbols
pseudo.PRJNA1141235.exp <- subset(pseudo.PRJNA1141235.exp, !(is.na(symbols)) & !(duplicated(symbols)))
rownames(pseudo.PRJNA1141235.exp) <- pseudo.PRJNA1141235.exp$symbols
pseudo.PRJNA1141235.exp <- pseudo.PRJNA1141235.exp[, c(1:4)]
ddsTC.pseudo.PRJNA1141235 <- DESeqDataSetFromMatrix(countData = round(as.matrix(pseudo.PRJNA1141235.exp), 0), 
                                                    colData = pseudo.PRJNA1141235.mdata, 
                                                    design = ~ 1)
keep.pseudo.PRJNA1141235 <- rowSums(counts(ddsTC.pseudo.PRJNA1141235) >= 1) >= ncol(counts(ddsTC.pseudo.PRJNA1141235))
table(keep.pseudo.PRJNA1141235)
ddsTC.pseudo.PRJNA1141235 <- ddsTC.pseudo.PRJNA1141235[keep.pseudo.PRJNA1141235, ]
vst.pseudo.PRJNA1141235 <- DESeq2::vst(ddsTC.pseudo.PRJNA1141235, blind = TRUE)
vst.pseudo.PRJNA1141235.mtx <- assay(vst.pseudo.PRJNA1141235)
boxplot(vst.pseudo.PRJNA1141235.mtx)

endest.results.pseudo.PRJNA1141235 <- estimate_cycle_time(exprs = vst.pseudo.PRJNA1141235.mtx,
                                                          entrez_ids = mapIds(org.Hs.eg.db,
                                                          keys = rownames(vst.pseudo.PRJNA1141235.mtx),
                                                          column = 'ENTREZID',
                                                          keytype = 'SYMBOL',
                                                          multiVals = first))
pseudo.PRJNA1141235.mdata$EndEst <- endest.results.pseudo.PRJNA1141235$estimated_time

hvg <- names(tail(sort(rowVars(vst.pseudo.PRJNA1141235.mtx)), 500))
pc.pseudo.PRJNA1141235 <- prcomp(t(vst.pseudo.PRJNA1141235.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pseudo.PRJNA1141235.mdata$PC1 <- as.data.frame(pc.pseudo.PRJNA1141235[5]$x)$PC1
pseudo.PRJNA1141235.mdata$PC2 <- as.data.frame(pc.pseudo.PRJNA1141235[5]$x)$PC2

p4a11 <- ggplot(pseudo.PRJNA1141235.mdata, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.pseudo.PRJNA1141235))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.pseudo.PRJNA1141235))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24))

t.test(EndEst ~ Phenotype, pseudo.PRJNA1141235.mdata)
p4a12 <- ggplot(pseudo.PRJNA1141235.mdata, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'p = 0.492', x = '', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

p4a1 <- wrap_plots(p4a11, p4a12, ncol = 2, widths = c(2, 1))


# GSE183837
head(seu.GSE183837@meta.data)
pseudo.GSE183837.mdata <- seu.GSE183837@meta.data[ , c(4:10)] %>% distinct()
rownames(pseudo.GSE183837.mdata) <- pseudo.GSE183837.mdata$GSM.accession
pseudo.GSE183837.exp <- AggregateExpression(seu.GSE183837, group.by = c('Sample.name'), return.seurat = FALSE)$RNA
pseudo.GSE183837.exp <- as.data.frame(pseudo.GSE183837.exp)
features.ids.GSE183837 <- features.ids[rownames(features.ids) %in% rownames(pseudo.GSE183837.exp), ]
features.ids.GSE183837 <- features.ids.GSE183837[rownames(pseudo.GSE183837.exp), ]
pseudo.GSE183837.exp$symbols <- features.ids.GSE183837$symbols
pseudo.GSE183837.exp <- subset(pseudo.GSE183837.exp, !(is.na(symbols)) & !(duplicated(symbols)))
rownames(pseudo.GSE183837.exp) <- pseudo.GSE183837.exp$symbols
pseudo.GSE183837.exp <- pseudo.GSE183837.exp[, c(1:9)]
ddsTC.pseudo.GSE183837 <- DESeqDataSetFromMatrix(countData = round(as.matrix(pseudo.GSE183837.exp), 0), 
                                                 colData = pseudo.GSE183837.mdata, 
                                                 design = ~ 1)
keep.pseudo.GSE183837 <- rowSums(counts(ddsTC.pseudo.GSE183837) >= 1) >= ncol(counts(ddsTC.pseudo.GSE183837))
table(keep.pseudo.GSE183837)
ddsTC.pseudo.GSE183837 <- ddsTC.pseudo.GSE183837[keep.pseudo.GSE183837, ]
vst.pseudo.GSE183837 <- DESeq2::vst(ddsTC.pseudo.GSE183837, blind = TRUE)
vst.pseudo.GSE183837.mtx <- assay(vst.pseudo.GSE183837)
boxplot(vst.pseudo.GSE183837.mtx)

endest.results.pseudo.GSE183837 <- estimate_cycle_time(exprs = vst.pseudo.GSE183837.mtx,
                                                       entrez_ids = mapIds(org.Hs.eg.db,
                                                       keys = rownames(vst.pseudo.GSE183837.mtx),
                                                       column = 'ENTREZID',
                                                       keytype = 'SYMBOL',
                                                       multiVals = first))
pseudo.GSE183837.mdata$EndEst <- endest.results.pseudo.GSE183837$estimated_time

hvg <- names(tail(sort(rowVars(vst.pseudo.GSE183837.mtx)), 500))
pc.pseudo.GSE183837 <- prcomp(t(vst.pseudo.GSE183837.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pseudo.GSE183837.mdata$PC1 <- as.data.frame(pc.pseudo.GSE183837[5]$x)$PC1
pseudo.GSE183837.mdata$PC2 <- as.data.frame(pc.pseudo.GSE183837[5]$x)$PC2

p4a21 <- ggplot(pseudo.GSE183837.mdata, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.pseudo.GSE183837))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.pseudo.GSE183837))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  scale_shape_manual(values = c(21, 24))

t.test(EndEst ~ Phenotype, pseudo.GSE183837.mdata)
p4a22 <- ggplot(pseudo.GSE183837.mdata, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  # geom_point(size = 3, alpha = 0.4, position = position_dodge2(width = 0.4)) +
  labs(title = 'p = 0.655', x = 'Phenotype', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

p4a2 <- wrap_plots(p4a21, p4a22, ncol = 2, widths = c(2, 1))


# GSE250130
head(seu.GSE250130@meta.data)
pseudo.GSE250130.mdata <- seu.GSE250130@meta.data[ , c(4:10)] %>% distinct()
rownames(pseudo.GSE250130.mdata) <- pseudo.GSE250130.mdata$GSM.accession
pseudo.GSE250130.exp <- AggregateExpression(seu.GSE250130, group.by = c('orig.ident'), return.seurat = FALSE)$RNA
pseudo.GSE250130.exp <- as.data.frame(pseudo.GSE250130.exp)
features.ids.GSE250130 <- features.ids[rownames(features.ids) %in% rownames(pseudo.GSE250130.exp), ]
features.ids.GSE250130 <- features.ids.GSE250130[rownames(pseudo.GSE250130.exp), ]
pseudo.GSE250130.exp$symbols <- features.ids.GSE250130$symbols
pseudo.GSE250130.exp <- subset(pseudo.GSE250130.exp, !(is.na(symbols)) & !(duplicated(symbols)))
rownames(pseudo.GSE250130.exp) <- pseudo.GSE250130.exp$symbols
pseudo.GSE250130.exp <- pseudo.GSE250130.exp[, c(1:16)]
ddsTC.pseudo.GSE250130 <- DESeqDataSetFromMatrix(countData = round(as.matrix(pseudo.GSE250130.exp), 0), 
                                                 colData = pseudo.GSE250130.mdata, 
                                                 design = ~ 1)
keep.pseudo.GSE250130 <- rowSums(counts(ddsTC.pseudo.GSE250130) >= 1) >= ncol(counts(ddsTC.pseudo.GSE250130))
table(keep.pseudo.GSE250130)
ddsTC.pseudo.GSE250130 <- ddsTC.pseudo.GSE250130[keep.pseudo.GSE250130, ]
vst.pseudo.GSE250130 <- DESeq2::vst(ddsTC.pseudo.GSE250130, blind = TRUE)
vst.pseudo.GSE250130.mtx <- assay(vst.pseudo.GSE250130)
boxplot(vst.pseudo.GSE250130.mtx)

endest.results.pseudo.GSE250130 <- estimate_cycle_time(exprs = vst.pseudo.GSE250130.mtx,
                                                       entrez_ids = mapIds(org.Hs.eg.db,
                                                                           keys = rownames(vst.pseudo.GSE250130.mtx),
                                                                           column = 'ENTREZID',
                                                                           keytype = 'SYMBOL',
                                                                           multiVals = first))
pseudo.GSE250130.mdata$EndEst <- endest.results.pseudo.GSE250130$estimated_time

hvg <- names(tail(sort(rowVars(vst.pseudo.GSE250130.mtx)), 500))
pc.pseudo.GSE250130 <- prcomp(t(vst.pseudo.GSE250130.mtx[hvg, ]), center = TRUE, scale. = TRUE)
pseudo.GSE250130.mdata$PC1 <- as.data.frame(pc.pseudo.GSE250130[5]$x)$PC1
pseudo.GSE250130.mdata$PC2 <- as.data.frame(pc.pseudo.GSE250130[5]$x)$PC2

p4a31 <- ggplot(pseudo.GSE250130.mdata, aes(x = PC1, y = PC2, fill = EndEst, shape = Phenotype)) +
  geom_point(color = 'black', size = 3, alpha = 0.9) +
  labs(title = '', 
       x = paste0('PC1: ', round((summary(pc.pseudo.GSE250130))$importance[2, 1], 4) * 100, '% variance'),
       y = paste0('PC2: ', round((summary(pc.pseudo.GSE250130))$importance[2, 2], 4) * 100, '% variance'),
       color = 'EndEst', shape = 'Phenotype') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13)) +
  scale_fill_viridis_c() +
  # scale_fill_gradientn(colors = as.vector(pals::kovesi.cyclic_mrybm_35_75_c68(100))) + 
  # cyclic color scale would be hard to perceive in other cases, retained directed viridis color pallete
  scale_shape_manual(values = c(21, 24))

t.test(EndEst ~ Phenotype, pseudo.GSE250130.mdata)
p4a32 <- ggplot(pseudo.GSE250130.mdata, aes(Phenotype, EndEst)) +
  geom_boxplot(color = 'black', fill = '#f0f0f0') +
  geom_point(position = position_dodge2(width = 0.4), color = 'black', fill = '#bdbdbd', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'p = 0.968', x = '', y = 'EndEst') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13)) +
  coord_cartesian(ylim = c(0, 100))

p4a3 <- wrap_plots(p4a31, p4a32, ncol = 2, widths = c(2, 1))




# DEA for whole cell types pseudobulks
Group <- factor(pseudo.PRJNA1141235.mdata$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.PRJNA1141235.mtx, design)
fit <- eBayes(fit)
head(design)
dea.pseudo.PRJNA1141235.whole <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.PRJNA1141235.whole.up <- subset(dea.pseudo.PRJNA1141235.whole, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.PRJNA1141235.whole.down <- subset(dea.pseudo.PRJNA1141235.whole, logFC < 0.667 & P.Value < 0.05)

Group <- factor(pseudo.GSE183837.mdata$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE183837.mtx, design)
fit <- eBayes(fit)
head(design)
dea.pseudo.GSE183837.whole <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE183837.whole.up <- subset(dea.pseudo.GSE183837.whole, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE183837.whole.down <- subset(dea.pseudo.GSE183837.whole, logFC < 0.667 & P.Value < 0.05)

Group <- factor(pseudo.GSE250130.mdata$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE250130.mtx, design)
fit <- eBayes(fit)
head(design)
dea.pseudo.GSE250130.whole <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE250130.whole.up <- subset(dea.pseudo.GSE250130.whole, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE250130.whole.down <- subset(dea.pseudo.GSE250130.whole, logFC < 0.667 & P.Value < 0.05)

genesets.whole.up <- list(GSE250130 = rownames(dea.pseudo.GSE250130.whole.up),
                          GSE183837 = rownames(dea.pseudo.GSE183837.whole.up),
                          PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.whole.up))
genesets.whole.down <- list(GSE250130 = rownames(dea.pseudo.GSE250130.whole.down),
                            GSE183837 = rownames(dea.pseudo.GSE183837.whole.down),
                            PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.whole.down))

p4b11 <- ggVennDiagram(genesets.whole.up, 
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = 'Down-regulated genes:\nlog2 fold change < 0.667\n& p < 0.05', x = 'Whole\nsample', y = '', fill = '') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b12 <- ggVennDiagram(genesets.whole.down, 
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = 'Up-regulated genes:\nlog2 fold change > 0.667\n& p < 0.05', x = '', y = '', fill = 'Number\nof genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = 'right',
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13)) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

wrap_plots(p4b11, p4b12, ncol = 2)


# Separate cell type pseudobulks generation and DEA 
View(seu.PRJNA1141235.fig4@meta.data)
pseudo.PRJNA1141235.sep.mdata <- seu.PRJNA1141235.fig4@meta.data[ , c(4:10, 18)] %>% distinct()
rownames(pseudo.PRJNA1141235.sep.mdata) <- paste0(pseudo.PRJNA1141235.sep.mdata$GSM.accession, '_', 
                                                  pseudo.PRJNA1141235.sep.mdata$cell_type)
pseudo.PRJNA1141235.sep.exp <- AggregateExpression(seu.PRJNA1141235.fig4, 
                                                   group.by = c('orig.ident', 'cell_type'), return.seurat = FALSE)$RNA
pseudo.PRJNA1141235.sep.exp <- as.data.frame(pseudo.PRJNA1141235.sep.exp)
pseudo.PRJNA1141235.sep.exp <- pseudo.PRJNA1141235.sep.exp[, rownames(pseudo.PRJNA1141235.sep.mdata)]
features.ids.PRJNA1141235.sep <- features.ids[rownames(features.ids) %in% rownames(pseudo.PRJNA1141235.sep.exp), ]
features.ids.PRJNA1141235.sep <- features.ids.PRJNA1141235.sep[rownames(pseudo.PRJNA1141235.sep.exp), ]
pseudo.PRJNA1141235.sep.exp$symbols <- features.ids.PRJNA1141235$symbols
pseudo.PRJNA1141235.sep.exp <- subset(pseudo.PRJNA1141235.sep.exp, !(is.na(symbols)) & !(duplicated(symbols)))
rownames(pseudo.PRJNA1141235.sep.exp) <- pseudo.PRJNA1141235.sep.exp$symbols
pseudo.PRJNA1141235.sep.exp <- pseudo.PRJNA1141235.sep.exp[, c(1:20)]

levels(as.factor(pseudo.PRJNA1141235.sep.mdata$cell_type))
pseudo.PRJNA1141235.sep.mdata.Endothelial <- subset(pseudo.PRJNA1141235.sep.mdata, cell_type == 'Endothelial')
pseudo.PRJNA1141235.sep.exp.Endothelial <- pseudo.PRJNA1141235.sep.exp[ , rownames(pseudo.PRJNA1141235.sep.mdata.Endothelial)]
ddsTC.pseudo.PRJNA1141235.sep.Endothelial <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.PRJNA1141235.sep.exp.Endothelial), 0), 
  colData = pseudo.PRJNA1141235.sep.mdata.Endothelial, design = ~ 1)
keep.pseudo.PRJNA1141235.sep.Endothelial <- rowSums(counts(ddsTC.pseudo.PRJNA1141235.sep.Endothelial) >= 1) >= 4
table(keep.pseudo.PRJNA1141235.sep.Endothelial)
ddsTC.pseudo.PRJNA1141235.sep.Endothelial <- ddsTC.pseudo.PRJNA1141235.sep.Endothelial[keep.pseudo.PRJNA1141235.sep.Endothelial, ]
vst.pseudo.PRJNA1141235.sep.Endothelial <- DESeq2::vst(ddsTC.pseudo.PRJNA1141235.sep.Endothelial, blind = TRUE)
vst.pseudo.PRJNA1141235.sep.Endothelial.mtx <- assay(vst.pseudo.PRJNA1141235.sep.Endothelial)
Group <- factor(pseudo.PRJNA1141235.sep.mdata.Endothelial$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.PRJNA1141235.sep.Endothelial.mtx, design)
fit <- eBayes(fit)
dea.pseudo.PRJNA1141235.sep.Endothelial <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.PRJNA1141235.sep.Endothelial.up <- subset(dea.pseudo.PRJNA1141235.sep.Endothelial, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.PRJNA1141235.sep.Endothelial.down <- subset(dea.pseudo.PRJNA1141235.sep.Endothelial, logFC < 0.667 & P.Value < 0.05)

pseudo.PRJNA1141235.sep.mdata.Epithelial <- subset(pseudo.PRJNA1141235.sep.mdata, cell_type == 'Epithelial')
pseudo.PRJNA1141235.sep.exp.Epithelial <- pseudo.PRJNA1141235.sep.exp[ , rownames(pseudo.PRJNA1141235.sep.mdata.Epithelial)]
ddsTC.pseudo.PRJNA1141235.sep.Epithelial <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.PRJNA1141235.sep.exp.Epithelial), 0), 
  colData = pseudo.PRJNA1141235.sep.mdata.Epithelial, design = ~ 1)
keep.pseudo.PRJNA1141235.sep.Epithelial <- rowSums(counts(ddsTC.pseudo.PRJNA1141235.sep.Epithelial) >= 1) >= 4
table(keep.pseudo.PRJNA1141235.sep.Epithelial)
ddsTC.pseudo.PRJNA1141235.sep.Epithelial <- ddsTC.pseudo.PRJNA1141235.sep.Epithelial[keep.pseudo.PRJNA1141235.sep.Epithelial, ]
vst.pseudo.PRJNA1141235.sep.Epithelial <- DESeq2::vst(ddsTC.pseudo.PRJNA1141235.sep.Epithelial, blind = TRUE)
vst.pseudo.PRJNA1141235.sep.Epithelial.mtx <- assay(vst.pseudo.PRJNA1141235.sep.Epithelial)
Group <- factor(pseudo.PRJNA1141235.sep.mdata.Epithelial$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.PRJNA1141235.sep.Epithelial.mtx, design)
fit <- eBayes(fit)
dea.pseudo.PRJNA1141235.sep.Epithelial <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.PRJNA1141235.sep.Epithelial.up <- subset(dea.pseudo.PRJNA1141235.sep.Epithelial, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.PRJNA1141235.sep.Epithelial.down <- subset(dea.pseudo.PRJNA1141235.sep.Epithelial, logFC < 0.667 & P.Value < 0.05)

pseudo.PRJNA1141235.sep.mdata.Lymphoid <- subset(pseudo.PRJNA1141235.sep.mdata, cell_type == 'Lymphoid')
pseudo.PRJNA1141235.sep.exp.Lymphoid <- pseudo.PRJNA1141235.sep.exp[ , rownames(pseudo.PRJNA1141235.sep.mdata.Lymphoid)]
ddsTC.pseudo.PRJNA1141235.sep.Lymphoid <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.PRJNA1141235.sep.exp.Lymphoid), 0), 
  colData = pseudo.PRJNA1141235.sep.mdata.Lymphoid, design = ~ 1)
keep.pseudo.PRJNA1141235.sep.Lymphoid <- rowSums(counts(ddsTC.pseudo.PRJNA1141235.sep.Lymphoid) >= 1) >= 4
table(keep.pseudo.PRJNA1141235.sep.Lymphoid)
ddsTC.pseudo.PRJNA1141235.sep.Lymphoid <- ddsTC.pseudo.PRJNA1141235.sep.Lymphoid[keep.pseudo.PRJNA1141235.sep.Lymphoid, ]
vst.pseudo.PRJNA1141235.sep.Lymphoid <- DESeq2::vst(ddsTC.pseudo.PRJNA1141235.sep.Lymphoid, blind = TRUE)
vst.pseudo.PRJNA1141235.sep.Lymphoid.mtx <- assay(vst.pseudo.PRJNA1141235.sep.Lymphoid)
Group <- factor(pseudo.PRJNA1141235.sep.mdata.Lymphoid$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.PRJNA1141235.sep.Lymphoid.mtx, design)
fit <- eBayes(fit)
dea.pseudo.PRJNA1141235.sep.Lymphoid <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.PRJNA1141235.sep.Lymphoid.up <- subset(dea.pseudo.PRJNA1141235.sep.Lymphoid, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.PRJNA1141235.sep.Lymphoid.down <- subset(dea.pseudo.PRJNA1141235.sep.Lymphoid, logFC < 0.667 & P.Value < 0.05)

pseudo.PRJNA1141235.sep.mdata.Myeloid <- subset(pseudo.PRJNA1141235.sep.mdata, cell_type == 'Myeloid')
pseudo.PRJNA1141235.sep.exp.Myeloid <- pseudo.PRJNA1141235.sep.exp[ , rownames(pseudo.PRJNA1141235.sep.mdata.Myeloid)]
ddsTC.pseudo.PRJNA1141235.sep.Myeloid <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.PRJNA1141235.sep.exp.Myeloid), 0), 
  colData = pseudo.PRJNA1141235.sep.mdata.Myeloid, design = ~ 1)
keep.pseudo.PRJNA1141235.sep.Myeloid <- rowSums(counts(ddsTC.pseudo.PRJNA1141235.sep.Myeloid) >= 1) >= 4
table(keep.pseudo.PRJNA1141235.sep.Myeloid)
ddsTC.pseudo.PRJNA1141235.sep.Myeloid <- ddsTC.pseudo.PRJNA1141235.sep.Myeloid[keep.pseudo.PRJNA1141235.sep.Myeloid, ]
vst.pseudo.PRJNA1141235.sep.Myeloid <- DESeq2::vst(ddsTC.pseudo.PRJNA1141235.sep.Myeloid, blind = TRUE)
vst.pseudo.PRJNA1141235.sep.Myeloid.mtx <- assay(vst.pseudo.PRJNA1141235.sep.Myeloid)
Group <- factor(pseudo.PRJNA1141235.sep.mdata.Myeloid$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.PRJNA1141235.sep.Myeloid.mtx, design)
fit <- eBayes(fit)
dea.pseudo.PRJNA1141235.sep.Myeloid <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.PRJNA1141235.sep.Myeloid.up <- subset(dea.pseudo.PRJNA1141235.sep.Myeloid, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.PRJNA1141235.sep.Myeloid.down <- subset(dea.pseudo.PRJNA1141235.sep.Myeloid, logFC < 0.667 & P.Value < 0.05)

pseudo.PRJNA1141235.sep.mdata.Stromal <- subset(pseudo.PRJNA1141235.sep.mdata, cell_type == 'Stromal')
pseudo.PRJNA1141235.sep.exp.Stromal <- pseudo.PRJNA1141235.sep.exp[ , rownames(pseudo.PRJNA1141235.sep.mdata.Stromal)]
ddsTC.pseudo.PRJNA1141235.sep.Stromal <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.PRJNA1141235.sep.exp.Stromal), 0), 
  colData = pseudo.PRJNA1141235.sep.mdata.Stromal, design = ~ 1)
keep.pseudo.PRJNA1141235.sep.Stromal <- rowSums(counts(ddsTC.pseudo.PRJNA1141235.sep.Stromal) >= 1) >= 4
table(keep.pseudo.PRJNA1141235.sep.Stromal)
ddsTC.pseudo.PRJNA1141235.sep.Stromal <- ddsTC.pseudo.PRJNA1141235.sep.Stromal[keep.pseudo.PRJNA1141235.sep.Stromal, ]
vst.pseudo.PRJNA1141235.sep.Stromal <- DESeq2::vst(ddsTC.pseudo.PRJNA1141235.sep.Stromal, blind = TRUE)
vst.pseudo.PRJNA1141235.sep.Stromal.mtx <- assay(vst.pseudo.PRJNA1141235.sep.Stromal)
Group <- factor(pseudo.PRJNA1141235.sep.mdata.Stromal$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.PRJNA1141235.sep.Stromal.mtx, design)
fit <- eBayes(fit)
dea.pseudo.PRJNA1141235.sep.Stromal <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.PRJNA1141235.sep.Stromal.up <- subset(dea.pseudo.PRJNA1141235.sep.Stromal, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.PRJNA1141235.sep.Stromal.down <- subset(dea.pseudo.PRJNA1141235.sep.Stromal, logFC < 0.667 & P.Value < 0.05)




View(seu.GSE183837@meta.data)
pseudo.GSE183837.sep.mdata <- seu.GSE183837@meta.data[ , c(4:10, 18)] %>% distinct()
rownames(pseudo.GSE183837.sep.mdata) <- paste0(pseudo.GSE183837.sep.mdata$GSM.accession, '_', 
                                                  pseudo.GSE183837.sep.mdata$cell_type)
pseudo.GSE183837.sep.exp <- AggregateExpression(seu.GSE183837, 
                                                   group.by = c('Sample.name', 'cell_type'), return.seurat = FALSE)$RNA
pseudo.GSE183837.sep.exp <- as.data.frame(pseudo.GSE183837.sep.exp)
pseudo.GSE183837.sep.exp <- pseudo.GSE183837.sep.exp[, rownames(pseudo.GSE183837.sep.mdata)]
features.ids.GSE183837.sep <- features.ids[rownames(features.ids) %in% rownames(pseudo.GSE183837.sep.exp), ]
features.ids.GSE183837.sep <- features.ids.GSE183837.sep[rownames(pseudo.GSE183837.sep.exp), ]
pseudo.GSE183837.sep.exp$symbols <- features.ids.GSE183837$symbols
pseudo.GSE183837.sep.exp <- subset(pseudo.GSE183837.sep.exp, !(is.na(symbols)) & !(duplicated(symbols)))
rownames(pseudo.GSE183837.sep.exp) <- pseudo.GSE183837.sep.exp$symbols
pseudo.GSE183837.sep.exp <- pseudo.GSE183837.sep.exp[, c(1:45)]

levels(as.factor(pseudo.GSE183837.sep.mdata$cell_type))
pseudo.GSE183837.sep.mdata.Endothelial <- subset(pseudo.GSE183837.sep.mdata, cell_type == 'Endothelial')
pseudo.GSE183837.sep.exp.Endothelial <- pseudo.GSE183837.sep.exp[ , rownames(pseudo.GSE183837.sep.mdata.Endothelial)]
ddsTC.pseudo.GSE183837.sep.Endothelial <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE183837.sep.exp.Endothelial), 0), 
  colData = pseudo.GSE183837.sep.mdata.Endothelial, design = ~ 1)
keep.pseudo.GSE183837.sep.Endothelial <- rowSums(counts(ddsTC.pseudo.GSE183837.sep.Endothelial) >= 1) >= 4
table(keep.pseudo.GSE183837.sep.Endothelial)
ddsTC.pseudo.GSE183837.sep.Endothelial <- ddsTC.pseudo.GSE183837.sep.Endothelial[keep.pseudo.GSE183837.sep.Endothelial, ]
vst.pseudo.GSE183837.sep.Endothelial <- DESeq2::vst(ddsTC.pseudo.GSE183837.sep.Endothelial, blind = TRUE)
vst.pseudo.GSE183837.sep.Endothelial.mtx <- assay(vst.pseudo.GSE183837.sep.Endothelial)
Group <- factor(pseudo.GSE183837.sep.mdata.Endothelial$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE183837.sep.Endothelial.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE183837.sep.Endothelial <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE183837.sep.Endothelial.up <- subset(dea.pseudo.GSE183837.sep.Endothelial, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE183837.sep.Endothelial.down <- subset(dea.pseudo.GSE183837.sep.Endothelial, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE183837.sep.mdata.Epithelial <- subset(pseudo.GSE183837.sep.mdata, cell_type == 'Epithelial')
pseudo.GSE183837.sep.exp.Epithelial <- pseudo.GSE183837.sep.exp[ , rownames(pseudo.GSE183837.sep.mdata.Epithelial)]
ddsTC.pseudo.GSE183837.sep.Epithelial <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE183837.sep.exp.Epithelial), 0), 
  colData = pseudo.GSE183837.sep.mdata.Epithelial, design = ~ 1)
keep.pseudo.GSE183837.sep.Epithelial <- rowSums(counts(ddsTC.pseudo.GSE183837.sep.Epithelial) >= 1) >= 4
table(keep.pseudo.GSE183837.sep.Epithelial)
ddsTC.pseudo.GSE183837.sep.Epithelial <- ddsTC.pseudo.GSE183837.sep.Epithelial[keep.pseudo.GSE183837.sep.Epithelial, ]
vst.pseudo.GSE183837.sep.Epithelial <- DESeq2::vst(ddsTC.pseudo.GSE183837.sep.Epithelial, blind = TRUE)
vst.pseudo.GSE183837.sep.Epithelial.mtx <- assay(vst.pseudo.GSE183837.sep.Epithelial)
Group <- factor(pseudo.GSE183837.sep.mdata.Epithelial$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE183837.sep.Epithelial.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE183837.sep.Epithelial <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE183837.sep.Epithelial.up <- subset(dea.pseudo.GSE183837.sep.Epithelial, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE183837.sep.Epithelial.down <- subset(dea.pseudo.GSE183837.sep.Epithelial, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE183837.sep.mdata.Lymphoid <- subset(pseudo.GSE183837.sep.mdata, cell_type == 'Lymphoid')
pseudo.GSE183837.sep.exp.Lymphoid <- pseudo.GSE183837.sep.exp[ , rownames(pseudo.GSE183837.sep.mdata.Lymphoid)]
ddsTC.pseudo.GSE183837.sep.Lymphoid <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE183837.sep.exp.Lymphoid), 0), 
  colData = pseudo.GSE183837.sep.mdata.Lymphoid, design = ~ 1)
keep.pseudo.GSE183837.sep.Lymphoid <- rowSums(counts(ddsTC.pseudo.GSE183837.sep.Lymphoid) >= 1) >= 4
table(keep.pseudo.GSE183837.sep.Lymphoid)
ddsTC.pseudo.GSE183837.sep.Lymphoid <- ddsTC.pseudo.GSE183837.sep.Lymphoid[keep.pseudo.GSE183837.sep.Lymphoid, ]
vst.pseudo.GSE183837.sep.Lymphoid <- DESeq2::vst(ddsTC.pseudo.GSE183837.sep.Lymphoid, blind = TRUE)
pseudo.GSE183837.sep.Lymphoid.mtx <- assay(vst.pseudo.GSE183837.sep.Lymphoid)
Group <- factor(pseudo.GSE183837.sep.mdata.Lymphoid$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(pseudo.GSE183837.sep.Lymphoid.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE183837.sep.Lymphoid <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE183837.sep.Lymphoid.up <- subset(dea.pseudo.GSE183837.sep.Lymphoid, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE183837.sep.Lymphoid.down <- subset(dea.pseudo.GSE183837.sep.Lymphoid, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE183837.sep.mdata.Myeloid <- subset(pseudo.GSE183837.sep.mdata, cell_type == 'Myeloid')
pseudo.GSE183837.sep.exp.Myeloid <- pseudo.GSE183837.sep.exp[ , rownames(pseudo.GSE183837.sep.mdata.Myeloid)]
ddsTC.pseudo.GSE183837.sep.Myeloid <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE183837.sep.exp.Myeloid), 0), 
  colData = pseudo.GSE183837.sep.mdata.Myeloid, design = ~ 1)
keep.pseudo.GSE183837.sep.Myeloid <- rowSums(counts(ddsTC.pseudo.GSE183837.sep.Myeloid) >= 1) >= 4
table(keep.pseudo.GSE183837.sep.Myeloid)
ddsTC.pseudo.GSE183837.sep.Myeloid <- ddsTC.pseudo.GSE183837.sep.Myeloid[keep.pseudo.GSE183837.sep.Myeloid, ]
vst.pseudo.GSE183837.sep.Myeloid <- DESeq2::vst(ddsTC.pseudo.GSE183837.sep.Myeloid, blind = TRUE)
pseudo.GSE183837.sep.Myeloid.mtx <- assay(vst.pseudo.GSE183837.sep.Myeloid)
Group <- factor(pseudo.GSE183837.sep.mdata.Myeloid$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(pseudo.GSE183837.sep.Myeloid.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE183837.sep.Myeloid <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE183837.sep.Myeloid.up <- subset(dea.pseudo.GSE183837.sep.Myeloid, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE183837.sep.Myeloid.down <- subset(dea.pseudo.GSE183837.sep.Myeloid, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE183837.sep.mdata.Stromal <- subset(pseudo.GSE183837.sep.mdata, cell_type == 'Stromal')
pseudo.GSE183837.sep.exp.Stromal <- pseudo.GSE183837.sep.exp[ , rownames(pseudo.GSE183837.sep.mdata.Stromal)]
ddsTC.pseudo.GSE183837.sep.Stromal <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE183837.sep.exp.Stromal), 0), 
  colData = pseudo.GSE183837.sep.mdata.Stromal, design = ~ 1)
keep.pseudo.GSE183837.sep.Stromal <- rowSums(counts(ddsTC.pseudo.GSE183837.sep.Stromal) >= 1) >= 4
table(keep.pseudo.GSE183837.sep.Stromal)
ddsTC.pseudo.GSE183837.sep.Stromal <- ddsTC.pseudo.GSE183837.sep.Stromal[keep.pseudo.GSE183837.sep.Stromal, ]
vst.pseudo.GSE183837.sep.Stromal <- DESeq2::vst(ddsTC.pseudo.GSE183837.sep.Stromal, blind = TRUE)
pseudo.GSE183837.sep.Stromal.mtx <- assay(vst.pseudo.GSE183837.sep.Stromal)
Group <- factor(pseudo.GSE183837.sep.mdata.Stromal$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(pseudo.GSE183837.sep.Stromal.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE183837.sep.Stromal <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE183837.sep.Stromal.up <- subset(dea.pseudo.GSE183837.sep.Stromal, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE183837.sep.Stromal.down <- subset(dea.pseudo.GSE183837.sep.Stromal, logFC < 0.667 & P.Value < 0.05)





View(seu.GSE250130@meta.data)
pseudo.GSE250130.sep.mdata <- seu.GSE250130@meta.data[ , c(4:10, 19)] %>% distinct()
rownames(pseudo.GSE250130.sep.mdata) <- paste0(pseudo.GSE250130.sep.mdata$GSM.accession, '_', 
                                                  pseudo.GSE250130.sep.mdata$cell_type)
pseudo.GSE250130.sep.exp <- AggregateExpression(seu.GSE250130, 
                                                   group.by = c('orig.ident', 'cell_type'), return.seurat = FALSE)$RNA
pseudo.GSE250130.sep.exp <- as.data.frame(pseudo.GSE250130.sep.exp)
pseudo.GSE250130.sep.exp <- pseudo.GSE250130.sep.exp[, rownames(pseudo.GSE250130.sep.mdata)]
features.ids.GSE250130.sep <- features.ids[rownames(features.ids) %in% rownames(pseudo.GSE250130.sep.exp), ]
features.ids.GSE250130.sep <- features.ids.GSE250130.sep[rownames(pseudo.GSE250130.sep.exp), ]
pseudo.GSE250130.sep.exp$symbols <- features.ids.GSE250130$symbols
pseudo.GSE250130.sep.exp <- subset(pseudo.GSE250130.sep.exp, !(is.na(symbols)) & !(duplicated(symbols)))
rownames(pseudo.GSE250130.sep.exp) <- pseudo.GSE250130.sep.exp$symbols
pseudo.GSE250130.sep.exp <- pseudo.GSE250130.sep.exp[, c(1:79)]

levels(as.factor(pseudo.GSE250130.sep.mdata$cell_type))
pseudo.GSE250130.sep.mdata.Endothelial <- subset(pseudo.GSE250130.sep.mdata, cell_type == 'Endothelial')
pseudo.GSE250130.sep.exp.Endothelial <- pseudo.GSE250130.sep.exp[ , rownames(pseudo.GSE250130.sep.mdata.Endothelial)]
ddsTC.pseudo.GSE250130.sep.Endothelial <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE250130.sep.exp.Endothelial), 0), 
  colData = pseudo.GSE250130.sep.mdata.Endothelial, design = ~ 1)
keep.pseudo.GSE250130.sep.Endothelial <- rowSums(counts(ddsTC.pseudo.GSE250130.sep.Endothelial) >= 1) >= 4
table(keep.pseudo.GSE250130.sep.Endothelial)
ddsTC.pseudo.GSE250130.sep.Endothelial <- ddsTC.pseudo.GSE250130.sep.Endothelial[keep.pseudo.GSE250130.sep.Endothelial, ]
vst.pseudo.GSE250130.sep.Endothelial <- DESeq2::vst(ddsTC.pseudo.GSE250130.sep.Endothelial, blind = TRUE)
vst.pseudo.GSE250130.sep.Endothelial.mtx <- assay(vst.pseudo.GSE250130.sep.Endothelial)
Group <- factor(pseudo.GSE250130.sep.mdata.Endothelial$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE250130.sep.Endothelial.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE250130.sep.Endothelial <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE250130.sep.Endothelial.up <- subset(dea.pseudo.GSE250130.sep.Endothelial, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE250130.sep.Endothelial.down <- subset(dea.pseudo.GSE250130.sep.Endothelial, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE250130.sep.mdata.Epithelial <- subset(pseudo.GSE250130.sep.mdata, cell_type == 'Epithelial')
pseudo.GSE250130.sep.exp.Epithelial <- pseudo.GSE250130.sep.exp[ , rownames(pseudo.GSE250130.sep.mdata.Epithelial)]
ddsTC.pseudo.GSE250130.sep.Epithelial <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE250130.sep.exp.Epithelial), 0), 
  colData = pseudo.GSE250130.sep.mdata.Epithelial, design = ~ 1)
keep.pseudo.GSE250130.sep.Epithelial <- rowSums(counts(ddsTC.pseudo.GSE250130.sep.Epithelial) >= 1) >= 4
table(keep.pseudo.GSE250130.sep.Epithelial)
ddsTC.pseudo.GSE250130.sep.Epithelial <- ddsTC.pseudo.GSE250130.sep.Epithelial[keep.pseudo.GSE250130.sep.Epithelial, ]
vst.pseudo.GSE250130.sep.Epithelial <- DESeq2::vst(ddsTC.pseudo.GSE250130.sep.Epithelial, blind = TRUE)
vst.pseudo.GSE250130.sep.Epithelial.mtx <- assay(vst.pseudo.GSE250130.sep.Epithelial)
Group <- factor(pseudo.GSE250130.sep.mdata.Epithelial$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE250130.sep.Epithelial.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE250130.sep.Epithelial <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE250130.sep.Epithelial.up <- subset(dea.pseudo.GSE250130.sep.Epithelial, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE250130.sep.Epithelial.down <- subset(dea.pseudo.GSE250130.sep.Epithelial, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE250130.sep.mdata.Lymphoid <- subset(pseudo.GSE250130.sep.mdata, cell_type == 'Lymphoid')
pseudo.GSE250130.sep.exp.Lymphoid <- pseudo.GSE250130.sep.exp[ , rownames(pseudo.GSE250130.sep.mdata.Lymphoid)]
ddsTC.pseudo.GSE250130.sep.Lymphoid <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE250130.sep.exp.Lymphoid), 0), 
  colData = pseudo.GSE250130.sep.mdata.Lymphoid, design = ~ 1)
keep.pseudo.GSE250130.sep.Lymphoid <- rowSums(counts(ddsTC.pseudo.GSE250130.sep.Lymphoid) >= 1) >= 4
table(keep.pseudo.GSE250130.sep.Lymphoid)
ddsTC.pseudo.GSE250130.sep.Lymphoid <- ddsTC.pseudo.GSE250130.sep.Lymphoid[keep.pseudo.GSE250130.sep.Lymphoid, ]
vst.pseudo.GSE250130.sep.Lymphoid <- DESeq2::vst(ddsTC.pseudo.GSE250130.sep.Lymphoid, blind = TRUE)
vst.pseudo.GSE250130.sep.Lymphoid.mtx <- assay(vst.pseudo.GSE250130.sep.Lymphoid)
Group <- factor(pseudo.GSE250130.sep.mdata.Lymphoid$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE250130.sep.Lymphoid.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE250130.sep.Lymphoid <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE250130.sep.Lymphoid.up <- subset(dea.pseudo.GSE250130.sep.Lymphoid, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE250130.sep.Lymphoid.down <- subset(dea.pseudo.GSE250130.sep.Lymphoid, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE250130.sep.mdata.Myeloid <- subset(pseudo.GSE250130.sep.mdata, cell_type == 'Myeloid')
pseudo.GSE250130.sep.exp.Myeloid <- pseudo.GSE250130.sep.exp[ , rownames(pseudo.GSE250130.sep.mdata.Myeloid)]
ddsTC.pseudo.GSE250130.sep.Myeloid <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE250130.sep.exp.Myeloid), 0), 
  colData = pseudo.GSE250130.sep.mdata.Myeloid, design = ~ 1)
keep.pseudo.GSE250130.sep.Myeloid <- rowSums(counts(ddsTC.pseudo.GSE250130.sep.Myeloid) >= 1) >= 4
table(keep.pseudo.GSE250130.sep.Myeloid)
ddsTC.pseudo.GSE250130.sep.Myeloid <- ddsTC.pseudo.GSE250130.sep.Myeloid[keep.pseudo.GSE250130.sep.Myeloid, ]
vst.pseudo.GSE250130.sep.Myeloid <- DESeq2::vst(ddsTC.pseudo.GSE250130.sep.Myeloid, blind = TRUE)
vst.pseudo.GSE250130.sep.Myeloid.mtx <- assay(vst.pseudo.GSE250130.sep.Myeloid)
Group <- factor(pseudo.GSE250130.sep.mdata.Myeloid$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE250130.sep.Myeloid.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE250130.sep.Myeloid <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE250130.sep.Myeloid.up <- subset(dea.pseudo.GSE250130.sep.Myeloid, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE250130.sep.Myeloid.down <- subset(dea.pseudo.GSE250130.sep.Myeloid, logFC < 0.667 & P.Value < 0.05)

pseudo.GSE250130.sep.mdata.Stromal <- subset(pseudo.GSE250130.sep.mdata, cell_type == 'Stromal')
pseudo.GSE250130.sep.exp.Stromal <- pseudo.GSE250130.sep.exp[ , rownames(pseudo.GSE250130.sep.mdata.Stromal)]
ddsTC.pseudo.GSE250130.sep.Stromal <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(pseudo.GSE250130.sep.exp.Stromal), 0), 
  colData = pseudo.GSE250130.sep.mdata.Stromal, design = ~ 1)
keep.pseudo.GSE250130.sep.Stromal <- rowSums(counts(ddsTC.pseudo.GSE250130.sep.Stromal) >= 1) >= 4
table(keep.pseudo.GSE250130.sep.Stromal)
ddsTC.pseudo.GSE250130.sep.Stromal <- ddsTC.pseudo.GSE250130.sep.Stromal[keep.pseudo.GSE250130.sep.Stromal, ]
vst.pseudo.GSE250130.sep.Stromal <- DESeq2::vst(ddsTC.pseudo.GSE250130.sep.Stromal, blind = TRUE)
vst.pseudo.GSE250130.sep.Stromal.mtx <- assay(vst.pseudo.GSE250130.sep.Stromal)
Group <- factor(pseudo.GSE250130.sep.mdata.Stromal$Phenotype)
design <- model.matrix(~ Group)
fit <- lmFit(vst.pseudo.GSE250130.sep.Stromal.mtx, design)
fit <- eBayes(fit)
dea.pseudo.GSE250130.sep.Stromal <- topTable(fit, coef = 2, adjust = 'BH', number = Inf)
dea.pseudo.GSE250130.sep.Stromal.up <- subset(dea.pseudo.GSE250130.sep.Stromal, logFC > 0.667 & P.Value < 0.05)
dea.pseudo.GSE250130.sep.Stromal.down <- subset(dea.pseudo.GSE250130.sep.Stromal, logFC < 0.667 & P.Value < 0.05)



genesets.Endothelial.up <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Endothelial.up),
                          GSE183837 = rownames(dea.pseudo.GSE183837.sep.Endothelial.up),
                          PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Endothelial.up))
genesets.Endothelial.down <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Endothelial.down),
                            GSE183837 = rownames(dea.pseudo.GSE183837.sep.Endothelial.down),
                            PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Endothelial.down))

genesets.Epithelial.up <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Epithelial.up),
                          GSE183837 = rownames(dea.pseudo.GSE183837.sep.Epithelial.up),
                          PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Epithelial.up))
genesets.Epithelial.down <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Epithelial.down),
                            GSE183837 = rownames(dea.pseudo.GSE183837.sep.Epithelial.down),
                            PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Epithelial.down))

genesets.Lymphoid.up <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Lymphoid.up),
                          GSE183837 = rownames(dea.pseudo.GSE183837.sep.Lymphoid.up),
                          PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Lymphoid.up))
genesets.Lymphoid.down <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Lymphoid.down),
                            GSE183837 = rownames(dea.pseudo.GSE183837.sep.Lymphoid.down),
                            PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Lymphoid.down))

genesets.Myeloid.up <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Myeloid.up),
                          GSE183837 = rownames(dea.pseudo.GSE183837.sep.Myeloid.up),
                          PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Myeloid.up))
genesets.Myeloid.down <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Myeloid.down),
                            GSE183837 = rownames(dea.pseudo.GSE183837.sep.Myeloid.down),
                            PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Myeloid.down))

genesets.Stromal.up <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Stromal.up),
                          GSE183837 = rownames(dea.pseudo.GSE183837.sep.Stromal.up),
                          PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Stromal.up))
genesets.Stromal.down <- list(GSE250130 = rownames(dea.pseudo.GSE250130.sep.Stromal.down),
                            GSE183837 = rownames(dea.pseudo.GSE183837.sep.Stromal.down),
                            PRJNA1141235 = rownames(dea.pseudo.PRJNA1141235.sep.Stromal.down))

p4b21 <- ggVennDiagram(genesets.Endothelial.up, 
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b22 <- ggVennDiagram(genesets.Endothelial.down, 
                      label = 'count', 
                      label_alpha = 0, 
                      label_size = 4.7, 
                      set_size = 0, 
                      color = 1,
                      edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b31 <- ggVennDiagram(genesets.Epithelial.up, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b32 <- ggVennDiagram(genesets.Epithelial.down, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b41 <- ggVennDiagram(genesets.Lymphoid.up, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) + 
  scale_y_reverse()

p4b42 <- ggVennDiagram(genesets.Lymphoid.down, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b51 <- ggVennDiagram(genesets.Myeloid.up, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) + 
  scale_y_reverse()

p4b52 <- ggVennDiagram(genesets.Myeloid.down, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b61 <- ggVennDiagram(genesets.Stromal.up, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()

p4b62 <- ggVennDiagram(genesets.Stromal.down, 
                       label = 'count', 
                       label_alpha = 0, 
                       label_size = 4.7, 
                       set_size = 0, 
                       color = 1,
                       edge_size = 0.5) + 
  labs(title = '', x = '', y = '', fill = 'Number of genes') + 
  theme_void(base_size = 13) +
  theme(plot.title = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_gradient(low = '#F4FAFE', high = '#525252', limits = c(0, 1088)) +
  scale_y_reverse()



## Figure 4 
p4a <- wrap_plots(p4a1, p4a2, p4a3, nrow = 3)
p4b <- wrap_plots(p4b11, p4b12,
                  p4b21, p4b22,
                  p4b31, p4b32,
                  p4b41, p4b42,
                  p4b51, p4b52,
                  p4b61, p4b62, nrow = 6)
p4 <- wrap_plots(p4a, p4b, ncol = 2, widths = c(1, 1.5))
ggsave(plot = p4, filename = 'visualization/Figure 4 (1).png', width = 13, height = 9, dpi = 300)
ggsave(plot = p4, filename = 'visualization/Figure 4 (2).png', width = 13, height = 13, dpi = 300)






## Graphical abstract
mdata.GSE111974 <- read.csv('processed.data/mdata.GSE111974.csv', sep = ',', header = T, row.names = 1)
exp.GSE111974 <- read.csv('processed.data/exp.GSE111974.max.csv', sep = ',', header = T, row.names = 1)
mdata.GSE111974$PDZK1 <- as.numeric(exp.GSE111974['PDZK1', ])

pgabs1 <- ggplot(mdata.GSE111974, aes(Phenotype, PDZK1, fill = Phenotype)) +
  geom_boxplot(color = 'black', alpha = 0.5) + 
  geom_point(position = position_dodge2(width = 0.4), color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'Significant\ndifference', x = 'Phenotype', y = 'mRNA levels') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 13),
        legend.position = 'none') +
  scale_fill_brewer(palette = 'Set2')

pgabs2 <- ggplot(mdata.GSE111974, aes(EndEst, PDZK1, color = Phenotype, fill = Phenotype)) +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  geom_point(color = 'black', size = 3, shape = 21, alpha = 0.8) +
  labs(title = 'Non-significant\ndifference', x = 'Menstrual cycle progression', y = 'mRNA levels') +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = 'italic'),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13),
        legend.text = element_text(size = 13), 
        legend.title = element_text(size = 13),
        legend.position = 'bottom') +
  scale_colour_brewer(palette = 'Set2') +
  scale_fill_brewer(palette = 'Set2')

wrap_plots(pgabs1, pgabs2, ncol = 2, widths = c(1, 3)) +
  plot_annotation(title = 'Differential expression analysis for an example gene\nbefore and after considering the menstrual cycle effect',
                  theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5)))
ggsave('visualization/GA.png', width = 6, height = 4.5, dpi = 600, limitsize = FALSE) 
  
