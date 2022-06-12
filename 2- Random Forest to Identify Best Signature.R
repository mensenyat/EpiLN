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
