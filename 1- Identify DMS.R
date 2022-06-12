library(splitTools)
library(varSelRF)
library(tidyverse)
library(lamisc)
library(ROCR)
library(Rtsne)
library(M3C)
library(gplots)
library(GOfuncR)

ann = read.table("450k_Annotation.csv", sep = ",", header = TRUE) # BEDfile containing curated annotation of 450k

setwd("E:/BRCA/DiNome/1- Lymph Node Invasion/N0vsN1")
data = read.table("BRCA_450k.csv", sep = ",", header = TRUE, row.names = 1)
colnames(data) = substr(colnames(data), 1, 12)
data = data[which(rownames(data) %in% ann$probe),]

meta = read.table("TCGA_Clinical_Data.txt", sep = "\t", header = TRUE, row.names = 1) # File containing curated Clinical Data from TCGA
meta = meta[which(meta$CPE > 0.65),]

barc = read.table("Barcodes to exclude.txt", sep = "\t") # Barcodes of patients with aberrant DNAm pattern previously obtained
meta = meta[which(!rownames(meta) %in% barc$V1),] 
meta = meta[which(rownames(meta) %in% colnames(data)),]

meta$LN = "LN" # All non-N0 are considered LN
meta$LN[which(substr(meta$AJCC_N,1,2) == "N0")] = "N0" # All types of N0 are considered N0
data = data[which(colnames(data) %in% rownames(meta))]
data = data[,order(match(colnames(data), rownames(meta)))]

data = na.omit(data) 
data_t = as.data.frame(as.matrix(t(data)))

#### Identify DEG ####
meta1 = meta[which(meta$LN == "LN"),]
meta0 = meta[which(meta$LN == "N0"),]

data1 = data[which(colnames(data) %in% rownames(meta1))]
data0 = data[which(colnames(data) %in% rownames(meta0))]

mean1 = apply(data1, 1, mean)
mean0 = apply(data0, 1, mean)

fold = (mean1-mean0)

# Compute statistical significance #
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data1)) {
  x = data1[i,]
  y = data0[i,]
  
  t = wilcox.test(as.numeric(x), as.numeric(y)) 
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}

combined  = cbind(pvalue, fold, abs(fold), data)
combined_t = as.data.frame(as.matrix(t(combined)))

# Filter per fold and p-value #
fold_cutoff = 0.1
pvalue_cutoff = 0.01

# fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
#dim(data_t[filter_by_fold, ]) # Check how many probes have DM > 0.1

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
#dim(data_t[filter_by_pvalue, ]) # Check how many probes have p-value < 0.01

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue 

filtered = data[filter_combined,]
filtered_t = as.data.frame(as.matrix(t(filtered)))
dim(filtered)
filtered_combined = combined[rownames(combined) %in% rownames(filtered),]
Probes_with_changes = filtered_combined[,1:3]
UpProbes = Probes_with_changes$Row.names[which(Probes_with_changes$fold>0)]
DownProbes = Probes_with_changes$Row.names[which(Probes_with_changes$fold<0)]

# Create dataframe with LN data #
meta$LN = meta$LN
plottingColors = c("Red", "cadetblue")
names(plottingColors) = unique(meta$LN)
