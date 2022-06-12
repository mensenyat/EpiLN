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
