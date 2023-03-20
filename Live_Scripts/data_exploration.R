# load packages
library(tximport)
library(DESeq2)
library(tidyverse)

# read the sample information into a data frame
sampleinfo <- read_tsv("data/samplesheet.tsv", col_types=c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)

# import count data
# construct path to Salmon quantification file
files <- file.path("salmon", sampleinfo$SampleName, "quant.sf")
# match sample names in the sampleinfo table
files <- set_names(files, sampleinfo$SampleName)

# relate transcript ID to gene ID
tx2gene <- read_tsv("references/tx2gene.tsv")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

head(txi$counts)
# save the txi object for use in later sessions
saveRDS(txi, file="salmon_outputs/txi.rds")

# Exercise 1
tpm <- tximport(files, type = "salmon", tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM")

# prepare count matrix for data exploration
rawCounts <- round(txi$counts, 0)
# check dimension of count matrix
dim(rawCounts)

# filter the genes whose total number of reads across all samples is no greater than 5
keep <- rowSums(rawCounts) > 5
# summary of the "logical"
table(keep, useNA="always")
# subset genes where keep was TRUE
filtCounts <- rawCounts[keep,]
# check dimension of new count matrix
dim(filtCounts)

summary(filtCounts)
boxplot(filtCounts, main='Raw counts', las=2)
plot(rowMeans(filtCounts), rowSds(filtCounts),
     main='Raw counts: sd VS mean',
     xlim=c(0,10000), ylim=c(0,5000))

# simple log2 transformation
logcounts <- log2(filtCounts + 1)
# make a colour vector (vectorised if-else)
statusCols <- case_when(sampleinfo$Status=="Infected" ~ "red",
                        sampleinfo$Status=="Uninfected" ~ "orange")

# check distributions of samples using boxplots
boxplot(logcounts, xlab='', ylab='Log2(Counts)', las=2, col=statusCols, main="Log2(Counts)")
# add a blue horizontal line that corresponds to the median
abline(h=median(logcounts), col="blue")

# log2 counts standard deviation vs mean expression
plot(rowMeans(logcounts), rowSds(logcounts),
     main="Log2 Counts: sd VS mean")

# VST: variance stabilizing transformation
vst_counts <- vst(filtCounts)
# check distribution using boxplot
boxplot(vst_counts, xlab="", ylab="VST coutns",
        las=2, col=statusCols)
# add a blue horizontal line corresponding to the median
abline(h=median(vst_counts), col="blue")
# VST counts sd VS mean expression
plot(rowMeans(vst_counts), rowSds(vst_counts),
     main = "VST counts: sd vs mean")

# Exercise 2
# rlog: regularized log transformation
rlog_counts <- rlog(filtCounts)
# check distribution using boxplot
boxplot(rlog_counts, xlab="", ylab="rlog counts",las=2, col=StatusCols)
abline(h=median(rlog_counts),col="blue")
# rlog counts sd VS mean expression
plot(rowMeans(rlog_counts), rowSds(rlog_counts),
     main="rlog counts: sd VS mean")

# Principal Component Analysis
library(ggfortify)
# run PCA
pcDat <- prcomp(t(rlog_counts))
# plot PCA
autoplot(pcDat)
# colour- and shape- code Cell Type and Status respectively
autoplot(pcDat, data=sampleinfo, colour="Status",
         shape="TimePoint", size=5)

# Exercise 3
autoplot(pcDat, data=sampleinfo, x=2, y=3,
         colour="Status", shape="TimePoint", size=5)

library(ggrepel)
autoplot(pcDat, data=sampleinfo, colour="Status", shape="TimePoint", size=5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=SampleName), box.padding = 0.8)

# use dplyr::mutate to fix the sample sheet
sampleinfo <- mutate(sampleinfo, Status=case_when(SampleName=="SRR7657882" ~ "Uninfected", SampleName=="SRR7657873" ~ "Infected", TRUE ~ Status))

# export the correct version for later use
write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")

# look at the PCA now
autoplot(pcDat, data=sampleinfo, colour="Status", shape="TimePoint", size=5)

# hierachical clustering
library(ggdendro)
hclDat <- t(rlog_counts) %>%
  dist(method="euclidean") %>%
  hclust()
ggdendrogram(hclDat, rotate=TRUE)

# replace labels for easier reading
hclDat2 <- hclDat
hclDat2$labels <- str_c(sampleinfo$Status, ":", sampleinfo$TimePoint)
ggdendrogram(hclDat2, rotate=TRUE)

