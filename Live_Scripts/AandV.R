#!/usr/bin/Rscript
# Description: Code for session 3 of RNa-seq course - annotation and visualisation

# Load libraries ----
library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(patchwork)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)

# set default theme for ggplot
theme_set(theme_bw())

# Load results from previous session ----
ddsObj.interaction <- readRDS('RObjects/DESeqDataSet.interaction.rds')
results.interaction.11 <- readRDS('RObjects/DESeqResults.interaction_d11.rds')
results.interaction.33 <- readRDS('RObjects/DESeqResults.interaction_d33.rds')

# Annotation ----
ah <- AnnotationHub()
# download the mouse annotation from ensembl
MouseEnsDb <- AnnotationHub::query(ah, c('EnsDb', 'Mus musculus', '102'))[[1]]
# extract gene information from the mouse annotation
annotations <- genes(MouseEnsDb, return.type='data.frame')

# simplify annottion table
annot <- annotations %>%
  select(gene_id, gene_name, entrezid) %>%
  filter(gene_id %in% rownames(results.interaction.11))

# load premade annotation table
ensemblAnnot <- readRDS('RObjects/Ensembl_annotations.rds')

# annotate my DESeq2 result ----
annot.interaction.11 <- as.data.frame(results.interaction.11) %>%
  rownames_to_column('GeneID') %>%          # move rownames as a column
  left_join(ensemblAnnot, 'GeneID') %>%     # merge annotation table to my results
  rename(logFC=log2FoldChange, FDR=padj)    # rename some columns

# save results
write_tsv(annot.interaction.11, 'results/Interaction.11_results_Annotated.txt')

# Visualisation ----
# p-value plot for sanity check
hist(annot.interaction.11$pvalue)
abline(v=0.5)

# Shrink my log2FC values to account for lowly expressed genes

ddsSchrink.11 <- lfcShrink(dds = ddsObj.interaction,
                           res = results.interaction.11,
                           type = 'ashr')


schrinkTab.11 <- as.data.frame(ddsSchrink.11) %>%
  rownames_to_column('GeneID') %>%          # move rownames as a column
  left_join(ensemblAnnot, 'GeneID') %>%     # merge annotation table to my results
  rename(logFC=log2FoldChange, FDR=padj)    # rename some columns


# Plotting, finally!
par(mfrow=c(2,1)) # plot formatting options
plotMA(results.interaction.11, alpha=0.05)
plotMA(ddsSchrink.11, alpha=0.05)

## ggplot2
# p-value histogram ggplot version
ggplot(annot.interaction.11, aes(pvalue)) + 
  geom_histogram(col='black', fill='grey40') + 
  theme_bw()

volcanoTab.11 <- schrinkTab.11 %>%
  mutate(`-log10(pvalue)`=-log10(pvalue))

ggplot(volcanoTab.11, 
       aes(x=logFC, y=`-log10(pvalue)`)) +
  geom_point(aes(col=FDR<0.05), size=1) +
  geom_text(data=~top_n(.x, 1,wt=-FDR), 
            aes(label=Symbol), size=15)



# exercise 1 ----
# repeat analysis with d33
## shrink
ddsSchrink.33 <- lfcShrink(dds = ddsObj.interaction,
                           res = results.interaction.33,
                           type = 'ashr')

## create -log10(pvalue)
# first create my table with the results
schrinkTab.33 <- as.data.frame(ddsSchrink.33) %>%
  rownames_to_column('GeneID') %>%          # move rownames as a column
  left_join(ensemblAnnot, 'GeneID') %>%     # merge annotation table to my results
  rename(logFC=log2FoldChange, FDR=padj)    # rename some columns

volcanoTab.33 <- schrinkTab.33 %>%
  mutate(`-log10(pvalue)`=-log10(pvalue))


## Volcano plot
v33 <- ggplot(volcanoTab.33, 
       aes(x=logFC, y=`-log10(pvalue)`)) +
  geom_point(aes(col=FDR<0.05), size=1) +
  geom_text_repel(data=~top_n(.x, 10,wt=-FDR), 
            aes(label=Symbol), size=5) +
  geom_text_repel(data=~top_n(.x, 1,wt=-logFC), 
            aes(label=Symbol), size=5)


## Compare
v11 <- ggplot(volcanoTab.11, 
       aes(x=logFC, y=`-log10(pvalue)`)) +
  geom_point(aes(col=FDR<0.05), size=1) +
  geom_text_repel(data=~top_n(.x, 10,wt=-FDR), 
                  aes(label=Symbol), size=5) +
  geom_text_repel(data=~top_n(.x, 3,wt=-logFC), 
                  aes(label=Symbol), size=5)


v11 + v33
v11 / v33

# exercise 2 ----
## MA plot
par(mfrow=c(1,2))
plotMA(results.interaction.33, alpha=0.05)
plotMA(ddsSchrink.33, alpha=0.05)

# repeat plot but in ggplot format
maTab.33 <- schrinkTab.33 %>%
  mutate(`M`=log2(baseMean))

ggplot(maTab.33, aes(x=M, y= logFC)) +
  geom_point(aes(col=FDR < 0.05), size=1) +
  scale_y_continuous(limit=c(-6,6), oob = scales::squish) + 
  theme_bw()

# Strip charts ----


geneID <- filter(schrinkTab.11, Symbol=='Il10ra') %>%
  pull(GeneID)


plotCounts(ddsObj.interaction,
           gene = geneID,
           intgroup = c('TimePoint', 'Status', 'Replicate'),
           returnData = T) %>%
  ggplot(aes(x=Status, y=log2(count))) +
  geom_point(aes(fill=Replicate), shape=21, size=2) +
  facet_wrap(~TimePoint) +
  expand_limits(y=0) +
  labs(title = 'Normalised counts - Il10ra') + 
  theme_bw()


# Venn diagram ----

vennDat <- tibble(Geneid=rownames(results.interaction.11)) %>%
  
  mutate(Upregulated_11 = results.interaction.11$padj < 0.05 &
           !is.na(results.interaction.11$padj) &
           results.interaction.11$log2FoldChange > 0) %>% 
  mutate(Downregulated_11 = results.interaction.11$padj < 0.05 &
           !is.na(results.interaction.11$padj) &
           results.interaction.11$log2FoldChange < 0) %>%
  
  mutate(Upregulated_33 = results.interaction.33$padj < 0.05 & 
           !is.na(results.interaction.33$padj) & 
           results.interaction.33$log2FoldChange > 0) %>%
  
  mutate(Downregulated_33 = results.interaction.33$padj < 0.05 & 
           !is.na(results.interaction.33$padj) & 
           results.interaction.33$log2FoldChange < 0) 

ggvenn(vennDat, set_name_size = 4)


# Heatmap ----
# get top significant genes
sigGenes <- schrinkTab.11 %>%
  top_n(300, wt=-FDR) %>%
  pull('GeneID')

# filter the data for these genes
plotDat <- vst(ddsObj.interaction)[sigGenes,] %>%
  assay()

z.mat <- t(scale(t(plotDat), center = T, scale = T))

Heatmap(plotDat)

# colour palette
myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette)

Heatmap(z.mat, 
        name='z-scores', 
        col=myRamp, 
        show_row_names = F)


ha <- HeatmapAnnotation(df=colData(ddsObj.interaction)[,c('Status', 'TimePoint')])


Heatmap(z.mat, 
        name='z-scores', 
        col=myRamp, 
        show_row_names = F,
        top_annotation = ha)

# specify colors for annotation
ha1 = HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")], 
                        col = list(Status = c("Uninfected" = "darkgreen", 
                                              "Infected" = "palegreen"), 
                                   TimePoint = c("d11" = "lightblue", 
                                                 "d33" = "darkblue")))

Heatmap(z.mat, 
        name='z-scores', 
        col=myRamp, 
        show_row_names = F,
        top_annotation = ha1, split = 3)



