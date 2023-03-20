# Last week we loaded int he data and resolved a sample swap by looking 
# at the PCA and checked some basic stats of the data. 

# and this  morning we have been through the statistical analysis of DESEq2 
# This afternoon we will put this into practice 

# we will start off with a simple model and build up through the afternoon 
# to the full interaction model and how to test  which model is the best
# the idea of this session is two familarise yourself witht eh functions 
# of deseq2 you are not a stastician and neither am I so do not panic. 


# For this course we use DEseq2, but other packages for RNASeq analysis are 
# availableincluding edgeR - there is a nice comparison here
# https://mikelove.wordpress.com/2016/09/28/deseq2-or-edger/



# Further to this there is a great bioconductor forum here:
# https://support.bioconductor.org/
# and MIkeLove Regulary responds 
# but also your local bioifnroamtician is good port of call (with biscuits)




###
# Completely restart environment this avoids any conflicts
# set fontsize correctly 
###

# first we are going to do the very simple model just as an example 
# we will do this together and get some results just to get the feel of 
# the deseq2 model design and object creation.


################################################################################
# load the libarries 
library(DESeq2)
library(tidyverse)
################################################################################

################################################################################
## Load in the data a setup ----

# You need 3 things to set a DESeq2 object
# you need the count data
# you need the meta-data 
# you need the model - what we would like to test is changing the gene expressoin


# we reload the data in from yesterdays session
txi <- readRDS("RObjects/txi.rds")
# along with the sample table
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", 
                       col_types = "cccc")

#first thing to do is check that the order of the rows of the object have the same 
# order as the sample table
all(colnames(txi$counts)==sampleinfo$SampleName)

################################################################################

# We have two variables in our experiment: “Status” and “Time Point”.

# We will fit two models under two assumptions: no interaction and interaction of 
# these two factors, however, to demonstrate the how DESeq2 is used we will start 
# with a simple model which considers Status but ignores Time Point.
################################################################################
## set up the model ( simple ) ----

# the formula uses a column from the sample table
sampleinfo
simple.model <- as.formula(~ Status)
#check the matrix model
# the intercept is chosen automatically via alphabetical
model.matrix(simple.model, data = sampleinfo)
# this is quite hard to see
cbind(sampleinfo, model.matrix(simple.model, data = sampleinfo))


# but we want the uninfected as the comparison level
sampleinfo <- mutate(sampleinfo, 
                     Status = fct_relevel(Status, "Uninfected"))
# now when we make the model matrix shows that the uninfected is 
# now the reference group
model.matrix(simple.model, data = sampleinfo)
cbind(sampleinfo, model.matrix(simple.model, data = sampleinfo))


# Build the deseq2 object - there are other options here for building the
# deseq2 object 
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)
# here explain about the deseq2 object and what is in it 
# it is a container that contains alot of things 




# have a look at the object
ddsObj.raw


# take the genes that have a count of more than 5 in ALL samples 
# this sets up a logcal vector the same size as the number of genes 
keep <- rowSums(counts(ddsObj.raw)) > 5
# keep these genes that have the high counts 
ddsObj.filt <- ddsObj.raw[keep,]
# these genes are the ones with low expression that would not be significant anyway 




################################################################################
# DESeq2 Workflow manual  ----

# Lets do the DESEq2 workflow step by step to try and understand what DESeq2 is 
# doing under the hood. in reality we won't do this but it is good to see 
# how DESeq2 works 

##### size factors ----

# we are not going assign a new variable 
# just so that we have the raw and filtered object as static variables

# how what this has done (added normalisation factor)
ddsObj <- estimateSizeFactors(ddsObj.filt)

# Not in the original data 
normalizationFactors(ddsObj.filt)
# it is in the new object we have created 
normalizationFactors(ddsObj)
# if you do not have transcript counts (from salmon) you will just 
# one normalisation factor per sample but we have one per gene per sample
# lets see the sum of the normalisation factors for each sample. 
colSums(normalizationFactors(ddsObj))

# get the log counts with no normalization
logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)

# array =5 chooses the sample with the worst normalization factors
# all other samples have lower expression!
# explain here about the :: when choosing a function from a specific package
limma::plotMA(logcounts, array = 5, ylim = c(-5,5))
abline(h=0, col = "red")

#but with normalization the counts now look good! 
# although we should be aware here that DESEq2 does not normalise the
# counts and then test for significance, it actually has these normalisation
# factors for each gene modeled into the GLM
logNormalisedCounts <- log2(counts(ddsObj, 
                                   normalized = TRUE) + 1)
limma::plotMA(logNormalisedCounts, array = 5, ylim = c(-5,5))
abline(h=0, col = "red")



##### dispersion ----
# never seen a bad dispertion
# we can do step 2 with the following command
# this estimates the dispersion for each gene using the neighbouring genes 
# as mentioned earlier 
ddsObj <- estimateDispersions(ddsObj)

# we can check them like this
plotDispEsts(ddsObj)
# remember, you are not expected to be a statistician and neither am I, 
# it is important to know what a normal one looks like
# dispertion decreases as mean increases (inversly proprtinal to the variance)
# most points are moved towards the fitted line (blue dots)
# few points as outliers 
# if it is not like this then you may have a serious problem and it is important 
# to go and speak to somebody


##### GLM and  Wald statistic ----

ddsObj <- nbinomWaldTest(ddsObj)



#  and that is it! the whole deseq2 workflow in 3 commands, but it is 
# even easier than that. and this is why we always have to be careful what 
# we put into DESeq2 


################################################################################
##### DESeq2 automatically  ----

ddsObj <- DESeq(ddsObj.filt)

# results

results.simple <- results(ddsObj, alpha = 0.05)
results.simple
# explain here the results table 


# now do exercise 1 10 minutes 
################################################################################
################################################################################
#### Exercise 1 ----


# how many genes are upregulated and down regulated at the 0.05 FDR cutoff level
sum(results.simple$padj < 0.05)

sum(is.na(results.simple$padj))
# independant filtering 


# these are genes that are too lowly expressed 

# up-regulated?

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm = TRUE)
#1879
# down-regulated?

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm = TRUE)
#1005

# explain here that it maybe silly to do this

################################################################################
# COFFE Break!
################################################################################

# Now the more complicated additive model 
# lets re-do the steps for loading in the data as this will help us remember 
# what is going on. Repetition is important for learning

# we need 3 things -
# the sampleInfo
# the counts
# the model 

#the counts
txi <- readRDS("RObjects/txi.rds")

# so let us load in the sample sheet and make sure we re-level
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", 
                       col_types = "cccc") %>%
    mutate(Status = fct_relevel(Status, "Uninfected"))


# now make a new model with the status also 
additive.model <- as.formula(~ TimePoint + Status)

# create the object 
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additive.model)

# keep the highly expressed genes 
keep <- rowSums(counts(ddsObj.raw)) > 5 
ddsObj.filt <- ddsObj.raw[keep,]


################################################################################
# Now do excersise 2 - 10 mins! 
# Run the size factor estimation, dispersion estimation and modelling steps 
#using the DESeq command as above.

# 1. run DESeq2
ddsObj.filt <- DESeq(ddsObj.filt)


# 2  Extract the default contrast using the results command into a new object
# called results.additive
results.additive <- results(ddsObj.filt, alpha = 0.05)


# 2. a which results?

results.additive
# you can see that the table is Infected vs Uninfected - why is this?
# it is because of the default results ordering 


# 2.b How many genes have an adjusted p-value of less than 0.05

sum(results.additive$padj < 0.05, na.rm = TRUE)
# 2766
# Due to the default contrasts we have the same contrasts as the simple model 

################################################################################
##### Default contrasts in DESEq 2 ----

# by checking the matrix we can see that the last column in the 
# model matrix is the results we receive and this 
# is the default
model.matrix(additive.model, data = sampleinfo)
cbind(sampleinfo, model.matrix(additive.model, data = sampleinfo))
resultsNames(ddsObj.filt)


# so lets just rename the infected / uninfected results
results.InfectedvUninfected <- results.additive
rm(results.additive)



################################################################################




################################################################################
##### Exercise 3 the other Beta of the model -  10 minutes -----

# check the manual for ?results to see how to use the contrast argument
?results

# 1. Retrieve the results for the contrast of d33 versus d11.
results.d33vd11 <- results(ddsObj.filt, 
                           name = "TimePoint_d33_vs_d11",
                           alpha = 0.05)

# 2. ow many differentially expressed genes are there at FDR < 0.05?
sum(results.d33vd11$padj < 0.05, na.rm = TRUE)


# so now we now how many genes change through the infection and through 
# the timePoints! 

################################################################################
# but should we use an interaction model - this will tell us those genes that are 
# changing through time and with a further unexpected difference in infection. 
# Interaction Model?

#  Go through the diagram

# What kind of genes would be interacting those that appear late or early stage 
# infection (interestign ones)


# lets have another look at the PCA that we looked at before 
# the first step is to transform the counts 
vstcounts <- vst(ddsObj.raw, blind = TRUE)
# then we plot the PCA
plotPCA(vstcounts, intgroup = c("Status", "TimePoint"))
# so what is happening here? we have spearation of infected but not uninfected


################################################################################
# LRT test for simple vs additive 
# the LRT takes two models and then reports how many genes 
# are better explained by this model that the more simple model

ddsObj.LRT <- DESeq(ddsObj.filt, test = "LRT", 
                    reduced = simple.model)

results.Additive_v_Simple <- results(ddsObj.LRT)
results.Additive_v_Simple

# for 66 genes - the more complex model is better 
sum(results.Additive_v_Simple$padj < 0.05, na.rm = TRUE)


# if a gene is significant here we can say that the more complicated model is 
# better at explaining the variance in gene expression

# here we only have 66 is this surprising given the PCA?









################################################################################
# Excercise 4 - 20 mins
# test the LRT with the interaction model vs the additive model
# the LRT takes two models and then reports how many genes 
# are better explained by this model that the more simple model 

# remeber what we need to sep up the DESE2 object
# sampleinfo
# gene counts (txi)
# model (interaction model)


# lets make a new interation model
interaction.model <- as.formula(~ TimePoint * Status)

# now re-create a new DESEq2 object
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = interaction.model)
# remove the low count genes
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

# any run DESeq2
ddsObj.interaction <- DESeq(ddsObj.filt)



# Now test 
ddsObj.LRT <- DESeq(ddsObj.interaction, test = "LRT", 
                    reduced = additive.model)

# get the results 
results.Interaction_v_Additive <- results(ddsObj.LRT, alpha = 0.05)

table(results.Interaction_v_Additive$padj < 0.05)
# 449 

# explain here that the 449 genes are in line with the PCA results
# we see that there is an interaction and these 449 genes are liklely very 
# exciting biologically 


################################################################################
# getting some degs from the interaction model 
resultsNames(ddsObj.interaction)

#“What is the difference in gene expression between Infected and Uninfected at 11 days post infection?”
results.interaction.11 <- results(ddsObj.interaction, 
                                  name="Status_Infected_vs_Uninfected",
                                  alpha=0.05)


#“What is the difference in gene expression between Infected and Uninfected at 33 days post infection?”
results.interaction.33 <- results(ddsObj.interaction, 
                                  contrast = list(c("Status_Infected_vs_Uninfected", 
                                                    "TimePointd33.StatusInfected")),
                                  alpha=0.05)

sum(results.interaction.11$padj < 0.05, na.rm = TRUE)
sum(results.interaction.33$padj < 0.05, na.rm = TRUE)



################################################################################
# Exercise 5

#Extract the results for d33 v d11 for Uninfected mice.
#How many genes have an adjusted p-value less than 0.05?
 #   Is this remarkable?

results.d33_v_d11_uninfected <- results(ddsObj.interaction,
                                        name = "TimePoint_d33_vs_d11",
                                        alpha = 0.05)
table(results.d33_v_d11_uninfected$padj < 0.05)

#Extract the results for d33 v d11 for Infected mice.
#How many genes have an adjusted p-value less than 0.05?
results.d33_v_d11_infected <- results(ddsObj.interaction,
                                      contrast = list(c("TimePoint_d33_vs_d11",
                                                        "TimePointd33.StatusInfected")),
                                      alpha = 0.05)
table(results.d33_v_d11_infected$padj < 0.05)


#Do these results suggest another approach to analysing this data set?
# the uninfected mice just live another few days their gene expression does not change


################################################################################
# Saving results
write_tsv(sampleinfo, "results/samplesheet_corrected.tsv")
saveRDS(ddsObj.interaction, "results/DESeqDataSet.interaction.rds")
saveRDS(results.interaction.11, "results/DESeqResults.interaction_d11.rds")
saveRDS(results.interaction.33, "results/DESeqResults.interaction_d33.rds")
