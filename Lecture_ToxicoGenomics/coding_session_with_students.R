#analysis
library(DESeq2)
library(PCAtools)

#plotting 
library(devtools)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(grid)
library(ggpubr)

sig_lvl <- .05  #you can adjust this

#set you working directory
#path <- choose.dir()
path <- "~/Desktop/Essen/theoretical_ecotox_course/"
setwd(path)
list.files()

#load the gene counts 
rawCounts <- read.table('lep_counts_filtered.csv', header = T, sep = ',', row.names= 1)
#inspect the data
head(rawCounts)
dim(rawCounts)


#change sample names
samples <- colnames(rawCounts)
sample_names <- gsub('X', '', samples)
sample_names

colnames(rawCounts) <- sample_names 
head(rawCounts)


#load the experimental data
sampleData <- read.csv('sample_treatments.csv', header = T, sep = ';', row.names = 1)
colnames(sampleData) <- c('Pesticide', 'Biotic_interaction', 'Treatment')

#sets reference levels for factors; very important
sampleData$Pesticide <- factor(sampleData$Pesticide, levels = c('control', 'low', 'medium', 'high'))
sampleData$Biotic_interaction <- factor(sampleData$Biotic_interaction, levels = c('L', 'G', 'GL'))


#if we want to model expression as a function of insecticide concentration with numerical predictors 
#(log fold change for each unit change in concentration), we must
#add insecticide conc. from day of sampling (exposure concencentration, not nominal conc.)
sampleData$Concentration[sampleData$Pesticide == 'control'] <- 0
sampleData$Concentration[sampleData$Pesticide == 'low'] <- 0.3
sampleData$Concentration[sampleData$Pesticide == 'medium'] <- 4.93
sampleData$Concentration[sampleData$Pesticide == 'high'] <- 19.5


#order the data
sampleData <- sampleData[sample_names, ]
colnames(rawCounts)
all(colnames(rawCounts) == rownames(sampleData)) #only if this is TRUE, proceed

#create DESeq2 object; rawCounts must be integers for DESeq (if they are not, use the round() function)
dds1 <- DESeqDataSetFromMatrix(countData=round(rawCounts), colData=sampleData, design= ~ Biotic_interaction*Concentration)

#Deseq2 object = summarized experiment object
#to access the count data:
head(assay(dds1))#is the same as
assays(dds1)$counts #note 'assay' vs. 'assays'; is the same as
head(counts(dds1))

#to access the metadata: 
colData(dds1)
rowData(dds1)

#before checking the actual gene expression changes, we explore the properties of count data and 
#check for outliers with a PCA. For this, we must normalize the counts with DESeq2 

#DESeq2 normalization
dds1 <- estimateSizeFactors(dds1)

#the heteroscedasticity of counts
#raw counts
ct <- counts(dds1, normalized = T) #this is a matrix that must be transformed to df
ct <- data.frame(ct)

#log transformed data
ct.log <- log2(ct)

#vst data
ct.vst <- data.frame(assay(vst(dds1)))

#explore mean-variance relationship, depending on the applied transformation
ct$gene <- row.names(ct)
ct$mean <- rowMeans(ct[,1:24])
ct$sd <- rowSds(as.matrix(ct[,1:24]))

ct.log$gene <- row.names(ct.log)
ct.log$mean <- rowMeans(ct.log[,1:24])
ct.log$sd <- rowSds(as.matrix(ct.log[,1:24]))


ct.vst$gene <- row.names(ct.vst)
ct.vst$mean <- rowMeans(ct.vst[,1:24])
ct.vst$sd <- rowSds(as.matrix(ct.vst[,1:24]))


#check
head(ct)
head(ct.log)
head(ct.vst)


# Reshape data from wide to long format for plotting
ct_long <- melt(ct, id.vars = c("gene", 'mean', 'sd'), variable.name = "sample", 
                value.name = "expression")

ct.log_long <- melt(ct.log, id.vars = c("gene", 'mean', 'sd'), variable.name = "sample", 
                    value.name = "expression")

ct.vst_long <- melt(ct.vst, id.vars = c("gene", 'mean', 'sd'), variable.name = "sample", 
                    value.name = "expression")

#check
head(ct_long)
head(ct.log_long)
head(ct.vst_long)


#subset the data to make sure we are looking at the same genes
ct_sub <- ct_long[ct_long$mean <= 10000,]
ct.log_sub <- ct.log_long[ct.log_long$gene %in% ct_sub$gene, ]
ct.vst_sub <- ct.vst_long[ct.vst_long$gene %in% ct_sub$gene, ]


#reduced data set
p1 <- ggplot(ct_sub, aes(x=mean, y=sd)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  #scale_x_continuous(limits = c(0,6000)) +
  #scale_y_continuous(limits = c(0,2000))+
  ggtitle("raw counts") +
  stat_smooth ()+
  theme_bw()
p1


##reduced data logged
p2 <- ggplot(ct.log_sub, aes(x=mean, y=sd)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  #scale_x_continuous(limits = c(1,15)) +
  #scale_y_continuous(limits = c(1,10))+
  ggtitle("logged counts") +
  stat_smooth()+
  theme_bw()
p2

##reduced data vst
p3 <- ggplot(ct.vst_sub, aes(x=mean, y=sd)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  #scale_x_continuous(limits = c(1,15)) +
  #scale_y_continuous(limits = c(1,10))+
  ggtitle("vst counts") +
  stat_smooth()+
  theme_bw()

p3


#plot all together
ggarrange(p1, p2, p3, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)


#now, we will run a PCA on normalized & transformed counts. The reason for the transformation is 
#the heteroscedasticity ib count data. If one performs PCA directly on raw RNA-seq counts, 
#the resulting plot typically depends mostly on the genes with highest counts (because they show 
#the largest absolute differences between samples) 
#log transformation requires the addition of a small pseudocount BUT depending on the choice 
#of this pseudocount, now the genes with the very lowest counts will contribute a considerable amount of 
#noise to the resulting plot, because taking the logarithm of small counts actually inflates their variance; 
#VST is designed to mitigate this problem a little bit, therefore we use the VST for PCA


#PCA on vst counts
vsd <- vst(dds1)

#define colors for the PCA plot
cols2 = c('grey', "#35B779FF", "#31688EFF","#440154FF")

#run the PCA
p <- pca(assay(vsd), metadata = sampleData)

#create a biplot (PC1 vs PC2)
biplot(p, showLoadings = F, sizeLoadingsNames = 2, boxedLoadingsNames = F, legendPosition = 'right', colby = 'Pesticide', colkey = cols2) 

#sample 35 might be an outlier, but we will leave it for now in the data set



#back the actual differential expression analysis with DESeq2
#run DESeq(); this includes normalization (done alrady), dispersion estimation and testing for DE
dds1 <- DESeq(dds1) 

#what we want to have in the end are expression changes between conditions; 
#for this we must specify contrasts
resultsNames(dds1)

#LFC estimate for 'Concentration' (was specified as numerical predictor) is given as 
#function of insecticide conc. = gene expression change [as LFC] per unit increase in insecticide [ug/L]
#(this is the effect of insecticide with biotic interaciton at reference level i.e, only lepidostoma present)
summary(results(dds1, name=c("Concentration"), alpha = sig_lvl))
conc <- results(dds1, name=c("Concentration"), alpha = sig_lvl)


#genes which are differentially expressed due to biotic interaction (with Concentration at reference level i.e., 0 ug/L)
summary(results(dds1, name=c("Biotic_interaction_GL_vs_L"), alpha = sig_lvl))

#genes showing significant interaction between biotic interaction & insecticide exposure 
summary(results(dds1, name=c("Biotic_interactionGL.Concentration"), alpha = sig_lvl))

#genes changing their expression as function of insecticide concentration under biotic interaction
#note, that this is not the combined effect
summary(results(dds1, contrast = list(c("Concentration", "Biotic_interactionGL.Concentration")), lfcThreshold = 0, alpha = sig_lvl))
conc_bi <- results(dds1, contrast = list(c("Concentration", "Biotic_interactionGL.Concentration")), lfcThreshold = 0, alpha = sig_lvl)


#despite the summary function, we can plot the results in different ways, e.g. with an MAplot 
plotMA(conc, ylim=c(-1,1))
plotMA(conc_bi, ylim=c(-1,1))



####functional enrichment to identidy the underlying molecular pathways
library(tidyr)
library(GO.db)
library(topGO)

setwd('~/Desktop/Doktorarbeit/Indoor_Genexpression/functional_annotation_eggnog/lepidostoma_trinity')

pcutoff = 0.05 ##you can adjust this

#######if gene2go is already prepared
GO_ids = read.csv('Trinity_Gene2GO.csv', sep=';', header = F)
##check that no invisible line breaks are in the file ;open the file one time a just save it
#transform to longformat

long_GO <- gather(GO_ids, Gen, IDs, V2:V673)

# take out genes without GO terms
long_GO <- long_GO[which(long_GO$`IDs` != ""),] 

##remove variable column
long_GO <- long_GO[, c(1, 3)]

#sort by transcript/gene
gene.go <- long_GO[order(long_GO$V1), ]

# Create list with element for each gene, containing vectors with all terms for each gene
gene2GO <- tapply(gene.go$`IDs`, gene.go$V1, function(x)x)

head(gene2GO) #go IDs as strings


#Define vector that is 1 if gene is significantly DE (adj.P.Val < chosen pcutoff) and 0 otherwise
##differentiate between up and downregulated genes
DE <- conc

DE$X <- row.names(DE)
DE <- DE[, c(7,2,6)]
DE$padj[is.na(DE$padj)] <- 1

DE_up <- DE
DE_up$up <-ifelse(DE$log2FoldChange < 0, 0, 1)
DE_up$padj <- ifelse(DE_up$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_up$padj == 1 & DE_up$up == 1, 1, 0)
geneList_up <- tmp


DE_down <- DE
DE_down$down <-ifelse(DE_down$log2FoldChange > 0, 0, 1)
DE_down$padj <- ifelse(DE_down$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_down$padj == 1 & DE_down$down == 1, 1, 0)
geneList_down <- tmp


##geneList need the same names as in the match for GO terms (gene2GO)
names(geneList_up) <- unlist(lapply(DE_up$X, function(x)x[1]))
names(geneList_down) <- unlist(lapply(DE_down$X, function(x)x[1]))

##Create topGOdata object:

GOdata_down <- new("topGOdata",
                   ontology = "BP",  #ontology criteria = biological process
                   allGenes = geneList_down, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                   geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

##run enrichment test: here Fishers Exact Test
resultFisher.weight.down <- runTest(GOdata_down, algorithm = "weight01", statistic = "fisher")
#summarize in table
down <- GenTable(GOdata_down, weight01 = resultFisher.weight.down, topNodes = length(resultFisher.weight.down@score), numChar = 500)


GOdata_up <- new("topGOdata",
                 ontology = "BP",  #ontology criteria = biological process
                 allGenes = geneList_up, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                 geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                 annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

resultFisher.weight.up <- runTest(GOdata_up, algorithm = "weight01", statistic = "fisher")
up <- GenTable(GOdata_up, weight01 = resultFisher.weight.up, topNodes = length(resultFisher.weight.up@score), numChar = 500)