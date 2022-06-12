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
