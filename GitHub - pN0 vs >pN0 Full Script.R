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

#### Identify DMS ####
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

#### Apply Machine Learning ####
part = partition(meta$LN, p = c(train = 2/3, valid = 1/3))
BarcTrain = rownames(meta)[part$train]
BarcVal = rownames(meta)[which(!rownames(meta) %in% BarcTrain)]

meta$Cohort = NA
meta$Cohort[which(rownames(meta) %in% BarcTrain)] = "Training"
meta$Cohort[which(rownames(meta) %in% BarcVal)] = "Validation"
table(meta$Cohort, meta$LN) # Check that the partition maintains the ratios in pN0 and >pN0

metaTrain = meta[which(meta$Cohort == "Training"),]
meta1 = metaTrain[which(metaTrain$LN == "LN"),]
meta0 = metaTrain[which(metaTrain$LN == "N0"),]
dataTrain = data[,which(colnames(data) %in% BarcTrain)]

metaVal = meta[which(meta$Cohort == "Validation"),]
dataVal = data[,which(colnames(data) %in% BarcVal)]

filteredTrain = dataTrain[which(rownames(dataTrain) %in% rownames(Probes_with_changes)),]
filteredTraint = as.data.frame(as.matrix(t(filteredTrain)))

filteredVal = dataVal[which(rownames(dataVal) %in% rownames(Probes_with_changes)),]
filteredValt = as.data.frame(as.matrix(t(filteredVal)))

# Use RandomForest #
RF_Genes = lapply(1:500,
                  function(x){
                    RF_Cl = varSelRF(filteredTraint, as.factor(metaTrain$LN), c.sd = 2, mtryFactor = 1, ntree = 1000,
                                     ntreeIterat = 500, vars.drop.num = NULL, vars.drop.frac = 0.2,
                                     whole.range = TRUE, recompute.var.imp = TRUE, verbose = FALSE,
                                     returnFirstForest = FALSE, fitted.rf = NULL, keep.forest = FALSE)
                    
                    sel_hist = RF_Cl$selec.history
                    sel_hist=sel_hist[order(sel_hist$Number.Variables),]
                    sel_hist=sel_hist[sel_hist$Number.Variables>=10,] # Select at least 10 probes to avoid over-fitting
                    sel_hist=sel_hist[order(sel_hist$OOB),]
                    sel_hist2=data.frame(x = sel_hist$Vars.in.Forest)
                    All_variables = sel_hist2 %>% separate(x, as.character(1:ncol(filteredTraint)))
                    sel_genes = na.omit(as.factor(All_variables[1,]))
                    sel_genes = droplevels.factor(sel_genes)
                    RF_genes = as.vector(sel_genes)
                    print(x)
                    
                    RF_genes # List containing the selected probes in each iteration of the process
                  })


# Select the most repeated combination #
RF_Genes_Sorted = lapply(1:500,
                         function(k){
                           sort(RF_Genes[[k]])
                         })  

n.obs = sapply(RF_Genes_Sorted, length)
seq.max = seq_len(max(n.obs))
RFComb = as.data.frame(t(sapply(RF_Genes_Sorted, "[", i = seq.max)))

repeated = count_duplicates(RFComb)
repeated = repeated[order(-as.numeric(repeated$dupe_count)),]
RFgenes = (repeated[1,-(ncol(repeated))])


#filtered_RF = filteredTraint[,which(colnames(filteredTraint) %in% RFgenes)]
#filtered_RF$Label = metaTrain$LN

#### Select combination with highest accuracy ####
# Re-train models #
AUCTraining = NA
AUCValidation = NA
AUCWhole = NA
lengthGenes = NA
repetitions = NA

for(i in 1:nrow(repeated)){ # Compute AUC of all combinations 
  RFgenes = (repeated[i,-(ncol(repeated))])
  lengthGenes[i] = length(RFgenes[which(!RFgenes == "<NA>")])
  repetitions[i] = repeated$dupe_count[i]
  filtered_RF = filteredTraint[,which(colnames(filteredTraint) %in% RFgenes)]
  filtered_RF = filtered_RF[order(match(rownames(filtered_RF), metaTrain$LN)),]
  filtered_RF$Label = metaTrain$LN
  
  fit.def <- randomForest(factor(Label) ~ ., data = filtered_RF, importance = TRUE, type = "classification", ntree = 10000)
  
  ROCR.simple_Train = list(predictions=fit.def$votes[,2], labels=filtered_RF$Label)
  pred_Train = prediction(ROCR.simple_Train$predictions, ROCR.simple_Train$labels)
  perf_Train = performance(pred_Train,"tpr","fpr")
  perfAUC_Train = performance(pred_Train,"tpr","fpr", measure="auc")
  AUCTraining[i] = perfAUC_Train@y.values
  
  # Validate #
  filteredVal = dataVal[which(rownames(dataVal) %in% RFgenes),]
  filteredValt = as.data.frame(as.matrix(t(filteredVal)))
  filteredValt = filteredValt[,order(match(colnames(filteredValt), colnames(filtered_RF)))]
  colnames(filteredValt) == colnames(filtered_def) # Check that probes are in the same order
  rownames(filteredValt) == rownames(metaVal) # Check that patients are in the same order

  GSEresult = predict(fit.def, filteredValt, type = "prob")
  ROCR.simple_Val = list(predictions=GSEresult[,2], labels=metaVal$LN)
  pred_Val = prediction(ROCR.simple_Val$predictions, ROCR.simple_Val$labels)
  perf_Val = performance(pred_Val,"tpr","fpr")
  perfAUC_Val = performance(pred_Val,"tpr","fpr", measure="auc")
  
  AUCValidation[i] = perfAUC_Val@y.values
  
  # Whole cohort
  ROCR.simple_Whole = list(predictions=c(fit.def$votes[,2], GSEresult[,2]), labels=c(metaTrain$LN, metaVal$LN))
  pred_Whole = prediction(ROCR.simple_Whole$predictions, ROCR.simple_Whole$labels)
  perf_Whole = performance(pred_Whole,"tpr","fpr")
  perfAUC_Whole = performance(pred_Whole,"tpr","fpr", measure="auc")
  AUCWhole[i] = perfAUC_Whole@y.values
  
  print(i)
}

summ = as.data.frame(cbind("N" = as.numeric(lengthGenes), # Summary of all combinations
                           "Repetitions" = as.numeric(repetitions),
                           "AUC Training" = as.numeric(AUCTraining), 
                           "AUC Validation" = as.numeric(AUCValidation),
                           "AUC Whole" = as.numeric(AUCWhole)))

summ = summ[order(-summ$Repetitions),]
summ = summ[order(-summ$`AUC Whole`),]

colorsAUC = c("Red", "cadetblue", "green", "gold", "darkorchid1") # Select colors to plot ROC Curves
position = c("bottomright", "bottomleft", "topright", "right", "topleft") # Add legends

pdf("5 ROC Curves Whole Cohort.pdf")
for(i in 1:5){
  l = as.numeric(rownames(summ)[i]) # "Barcode" directing to the combinations of genes
  RFgenes = repeated[which(rownames(repeated) == l),]
  RFgenes = (repeated[l,-(ncol(repeated))])
  lengthGenes[i] = length(RFgenes[which(!RFgenes == "<NA>")])
  filtered_RF = filteredTraint[,which(colnames(filteredTraint) %in% RFgenes)]
  filtered_RF = filtered_RF[order(match(rownames(filtered_RF), metaTrain$LN)),]
  filtered_RF$Label = metaTrain$LN
  
  fit.def <- randomForest(factor(Label) ~ ., data = filtered_RF, importance = TRUE, type = "classification", ntree = 10000)
  
  ROCR.simple_Train = list(predictions=fit.def$votes[,2], labels=filtered_RF$Label)
  pred_Train = prediction(ROCR.simple_Train$predictions, ROCR.simple_Train$labels)
  perf_Train = performance(pred_Train,"tpr","fpr")
  perfAUC_Train = performance(pred_Train,"tpr","fpr", measure="auc")
  
  # Validate #
  filteredVal = dataVal[which(rownames(dataVal) %in% RFgenes),]
  filteredValt = as.data.frame(as.matrix(t(filteredVal)))
  filteredValt = filteredValt[,order(match(colnames(filteredValt), colnames(filtered_RF)))]
  colnames(filteredValt) == colnames(filtered_def) # Check that probes are in the same order
  rownames(filteredValt) == rownames(metaVal) # Check that patients are in the same order
  
  GSEresult = predict(fit.def, filteredValt, type = "prob")
  ROCR.simple_Val = list(predictions=GSEresult[,2], labels=metaVal$LN)
  pred_Val = prediction(ROCR.simple_Val$predictions, ROCR.simple_Val$labels)
  perf_Val = performance(pred_Val,"tpr","fpr")
  perfAUC_Val = performance(pred_Val,"tpr","fpr", measure="auc")
  
  # Whole cohort
  ROCR.simple_Whole = list(predictions=c(fit.def$votes[,2], GSEresult[,2]), labels=c(metaTrain$LN, metaVal$LN))
  pred_Whole = prediction(ROCR.simple_Whole$predictions, ROCR.simple_Whole$labels)
  perf_Whole = performance(pred_Whole,"tpr","fpr")
  perfAUC_Whole = performance(pred_Whole,"tpr","fpr", measure="auc")
  plot(perf_Whole, col = colorsAUC[i], pce = 2, lwd = 4)
  legend(position[i], legend = paste("P", summ$N[i], " AUC=", round(as.numeric(perfAUC_Whole@y.values), 2), sep = ""), 
         fill=colorsAUC[i], border = T)
  
  par(new=TRUE)
}
dev.off()
par(new=FALSE)

par(mfrow=c(1,2))

#### AUC plots using the selected combination ####
RFgenes = repeated[which(rownames(repeated) == as.numeric(rownames(summ)[i])),]
par(mfrow=c(2,2))

filtered_RF = filteredTraint[,which(colnames(filteredTraint) %in% RFgenes)]
filtered_RF = filtered_RF[order(match(rownames(filtered_RF), metaTrain$LN)),]
filtered_RF$Label = metaTrain$LN

fit.def <- randomForest(factor(Label) ~ ., data = filtered_RF, importance = TRUE, type = "classification", ntree = 10000)
ROCR.simple_Train = list(predictions=fit.def$votes[,2], labels=filtered_RF$Label)
pred_Train = prediction(ROCR.simple_Train$predictions, ROCR.simple_Train$labels)
perf_Train = performance(pred_Train,"tpr","fpr")
perfAUC_Train = performance(pred_Train,"tpr","fpr", measure="auc")
perfAUC_Train@y.values
plot(perf_Train, main = paste("Training cohort - AUC = ", perfAUC_Train@y.values))


filteredVal = dataVal[which(rownames(dataVal) %in% RFgenes),]
filteredValt = as.data.frame(as.matrix(t(filteredVal)))
filteredValt = filteredValt[,order(match(colnames(filteredValt), colnames(filtered_RF)))]
colnames(filteredValt) == colnames(filtered_RF) # Check that probes are in the same order
rownames(filteredValt) == rownames(metaVal) # Check that patients are in the same order

GSEresult = predict(fit.def, filteredValt, type = "prob")
ROCR.simple_Val = list(predictions=GSEresult[,2], labels=metaVal$LN)
pred_Val = prediction(ROCR.simple_Val$predictions, ROCR.simple_Val$labels)
perf_Val = performance(pred_Val,"tpr","fpr")
perfAUC_Val = performance(pred_Val,"tpr","fpr", measure="auc")
plot(perf_Val, main = paste("Validation cohort - AUC = ", perfAUC_Val@y.values))

# Whole cohort #
ROCR.simple_Whole = list(predictions=c(fit.def$votes[,2], GSEresult[,2]), labels=c(metaTrain$LN, metaVal$LN))
pred_Whole = prediction(ROCR.simple_Whole$predictions, ROCR.simple_Whole$labels)
perf_Whole = performance(pred_Whole,"tpr","fpr")
perfAUC_Whole = performance(pred_Whole,"tpr","fpr", measure="auc")
perfAUC_Whole@y.values
plot(perf_Whole, main = paste("Whole cohort - AUC = ", perfAUC_Whole@y.values))



#### Repeat ROC curves using all DMS ####
filtered_RF = filteredTraint
filtered_RF$Label = metaTrain$LN

fit.def <- randomForest(factor(Label) ~ ., data = filtered_RF, importance = TRUE, type = "classification", ntree = 10000)
ROCR.simple_TrainDMS = list(predictions=fit.def$votes[,2], labels=filtered_RF$Label)
pred_TrainDMS = prediction(ROCR.simple_TrainDMS$predictions, ROCR.simple_TrainDMS$labels)
perf_TrainDMS = performance(pred_TrainDMS,"tpr","fpr")
perfAUC_TrainDMS = performance(pred_TrainDMS,"tpr","fpr", measure="auc")
perfAUC_TrainDMS@y.values
plot(perf_TrainDMS, main = paste("Training cohort - AUC = ", perfAUC_TrainDMS@y.values))

filteredValt = filtered_t[rownames(filtered_t) %in% rownames(metaVal),]
GSEresult = predict(fit.def, filteredValt, type = "prob")
ROCR.simple_ValDMS = list(predictions=GSEresult[,2], labels=metaVal$LN)
pred_ValDMS = prediction(ROCR.simple_ValDMS$predictions, ROCR.simple_ValDMS$labels)
perf_ValDMS = performance(pred_ValDMS,"tpr","fpr")
perfAUC_ValDMS = performance(pred_ValDMS,"tpr","fpr", measure="auc")
perfAUC_ValDMS@y.values
plot(perf_ValDMS, main = paste("Validation cohort - AUC = ", perfAUC_ValDMS@y.values))

ROCR.simple_WholeDMS = list(predictions=c(fit.def$votes[,2], GSEresult[,2]), labels=c(metaTrain$LN, metaVal$LN))
pred_WholeDMS = prediction(ROCR.simple_WholeDMS$predictions, ROCR.simple_WholeDMS$labels)
perf_WholeDMS = performance(pred_WholeDMS,"tpr","fpr")
perfAUC_WholeDMS = performance(pred_WholeDMS,"tpr","fpr", measure="auc")
perfAUC_WholeDMS@y.values
#pdf("ROC Curve Whole Cohort All DMS.pdf")
plot(perf_WholeDMS, main = paste("Whole cohort - AUC = ", perfAUC_WholeDMS@y.values))
dev.off()

#### Plots ####
# All DMS #
Tsne = Rtsne(as.matrix(filtered_t),perplexity = 5, theta = 0, max_iter = 5000)

#pdf('tSNE All patients All DMS.pdf')
plot(Tsne$Y, main="All DMS tSNE", xlab="t-SNE Component 1", ylab="t-SNE Component 2", col=plottingColors[meta$LN], pch=19)
legend('bottomleft', legend = unique(meta$LN), fill=plottingColors[unique(meta$LN)], border=T, title='Subtype')
dev.off()

pdf("UMAP All Patients All DMS.pdf")
umap(as.matrix(filtered), labels = meta$LN)
dev.off()

# Selected combination #
filtered_RF = filtered_t[,which(colnames(filtered_t) %in% RFgenes)]

Tsne = Rtsne(as.matrix(filtered_RF),perplexity = 15, theta = 0, max_iter = 5000)

#pdf('tSNE VarSelRF All patients PX.pdf')
plot(Tsne$Y, main="PX tSNE", xlab="t-SNE Component 1", ylab="t-SNE Component 2", col=plottingColors[meta$LN], pch=19)
legend('bottomleft', legend = unique(meta$LN), fill=plottingColors[unique(meta$LN)], border=T, title='Subtype')
dev.off()

pdf("UMAP VarSelRF All patients P37.pdf")
umap(as.matrix(t(filtered_RF)), labels = meta$LN)
dev.off()

# Volcano Plot #
pdf("Volcano_Plot.pdf")

plot(fold, -log10(pvalue), main = "pN0 vs >pN0 - Volcano", 
     xlim=c(-max(abs(fold)), max(abs(fold))), ylim=c(0, max(-log10(pvalue))), 
     xlab="Differential methylation", ylab = "-log10(pvalue)")

points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "blue")

abline(v = fold_cutoff, col = "red", lwd = 3)
abline(v = -fold_cutoff, col = "blue", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "grey", lwd = 3)
dev.off()

# Heatmap All DMS #
colv = as.dendrogram(hclust(as.dist(1-cor(t(as.matrix(filtered_t))))))
rowv = as.dendrogram(hclust(as.dist(1-cor(as.matrix(filtered_t)))))

pdf("Heatmap_All_DMS_Classical_GeneExp.pdf")
heatmap.2(as.matrix(t(filtered_t)), cexCol=0.7,
          col = rev(redblue(256)), Rowv=rowv, Colv=colv, scale = "row", 
          trace=NULL, tracecol=NULL, ColSideColors = plottingColors[meta$LN])
legend('topleft', legend = unique(meta$LN), fill = unique(plottingColors[meta$LN]), border=T, title='Subtype')
dev.off()


#### BEDfile using DMS ####
setwd("D:/R/BRCA/DiNome/1- Lymph Node Invasion/N0vsAll/Updated")
DMSbed = ann[ann$probe %in% rownames(filtered),]
write.table(DMSbed, file = "BEDfile DMS pN0 vs All.bed", sep = "\t", row.names = F)

Hypobed = ann[ann$probe %in% DownProbes,]
write.table(Hypobed, file = "BEDfile Hypo pN0 vs All.bed", sep = "\t", row.names = F)

Hyperbed = ann[ann$probe %in% UpProbes,]
write.table(Hyperbed, file = "BEDfile Hyper pN0 vs All.bed", sep = "\t", row.names = F)

# We charge these BEDfiles in MEME/T-Gene to identify genes affected by these CpG sites #

#### Gene Enrichment ####
hyper_genes = read.table("Hyper T-Gene/links.tsv", sep = "\t", header = T)
hyper_genes = hyper_genes[hyper_genes$CnD_P_Value < 0.01,]
hypo_genes = read.table("Hypo T-Gene/links.tsv", sep = "\t", header = T)
hypo_genes = hypo_genes[hypo_genes$CnD_P_Value < 0.01,]
all_genes = read.table("All DMS T-Gene/links.tsv", sep = "\t", header = T)
all_genes = all_genes[all_genes$CnD_P_Value < 0.01,]

genes1 = as.data.frame(cbind(hyper_genes$Gene_Name))
genes1$V2 = 1
enrHyper = go_enrich(genes1, test = 'hyper', n_randsets = 1000, organismDb = 'Homo.sapiens', gene_len = FALSE,
                     regions = FALSE, circ_chrom = FALSE, silent = T, domains = NULL, orgDb = NULL,
                     txDb = NULL, annotations = NULL, gene_coords = NULL,  godir = NULL)

genes2 = as.data.frame(cbind(hypo_genes$Gene_Name))
genes2$V2 = 1
enrHypo = go_enrich(genes2, test = 'hyper', n_randsets = 1000, organismDb = 'Homo.sapiens', gene_len = FALSE,
                    regions = FALSE, circ_chrom = FALSE, silent = T, domains = NULL, orgDb = NULL,
                    txDb = NULL, annotations = NULL, gene_coords = NULL,  godir = NULL)

genes3 = as.data.frame(cbind(all_genes$Gene_Name))
genes3$V2 = 1
enrAll = go_enrich(genes3, test = 'hyper', n_randsets = 1000, organismDb = 'Homo.sapiens', gene_len = FALSE,
                   regions = FALSE, circ_chrom = FALSE, silent = T, domains = NULL, orgDb = NULL,
                   txDb = NULL, annotations = NULL, gene_coords = NULL,  godir = NULL)

enrHypersign = enrHyper$results[which(enrHyper$results$raw_p_overrep < 0.05),]
enrHyposign = enrHypo$results[which(enrHypo$results$raw_p_overrep < 0.05),]
enrAllsign = enrAll$results[which(enrAll$results$raw_p_overrep < 0.05),]

refineHyper = refine(enrHyper, fwer=0.1)
refineHypo = refine(enrHypo, fwer=0.1)
refineAll = refine(enrAll, fwer=0.1)

refineHypersign = refineHyper[refineHyper$signif == T,]
refineHyposign = refineHypo[refineHypo$signif == T,]
refineAllsign = refineAll[refineAll$signif == T,]
