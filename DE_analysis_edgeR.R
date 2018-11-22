# Setup environment and load data
library(edgeR)
library(VennDiagram)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(goseq)
library(GO.db)

setwd("~/Dropbox/coauthored/Jenna_Jones/Analysis_sep2017/")

data <- read.table("Trinity_genes.counts.matrix", header = T, row.names = 1) 
data <- round(data)
data <- data[which(rowSums(cpm(data)>1)>=4),]


# Conduct differential expression analysis
groups <- data.frame(treatment=factor(c("Conc.00","Conc.00","Conc.00","Conc.00","Conc.05","Conc.05","Conc.05","Conc.05","Conc.10","Conc.10","Conc.10","Conc.10","Conc.30","Conc.30","Conc.30","Conc.30")))
group <- relevel(groups$treatment, ref = "Conc.00")
design <- model.matrix(~group)

y <- DGEList(counts = data, group = group)
y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)

fit <- glmFit(y, design)
colnames(fit)
Cv5.lrt <- glmLRT(fit, coef=2);Cv5_de <- row.names(Cv5.lrt$table)[which(p.adjust(Cv5.lrt$table$PValue, method = "BH")<0.05)]
Cv10.lrt <- glmLRT(fit, coef=3);Cv10_de <- row.names(Cv10.lrt$table)[which(p.adjust(Cv10.lrt$table$PValue, method = "BH")<0.05)]
Cv30.lrt <- glmLRT(fit, coef=4);Cv30_de <- row.names(Cv30.lrt$table)[which(p.adjust(Cv30.lrt$table$PValue, method = "BH")<0.05)]

DE_GENES <- unique(c(Cv5_de, Cv10_de, Cv30_de))
DE_TABLE <- matrix(ncol = 3, nrow = length(DE_GENES), dimnames = list(DE_GENES, c("LogFC_5", "LogFC_10", "LogFC_30")))
DE_TABLE[Cv5_de, 1] <- as.matrix(Cv5.lrt$table[Cv5_de,c(1)])
DE_TABLE[Cv10_de, 2] <- as.matrix(Cv10.lrt$table[Cv10_de,c(1)])
DE_TABLE[Cv30_de, 3] <- as.matrix(Cv30.lrt$table[Cv30_de,c(1)])

write.table(DE_TABLE, "DE_gene_table.txt", quote = F, sep = "\t")

# Produce venn diagram of shared DE genes among catgories
de_all <- intersect(intersect(Cv5_de, Cv10_de), Cv30_de)
de_5_10 <- intersect(Cv5_de, Cv10_de)
de_5_30 <- intersect(Cv5_de, Cv30_de)
de_10_30 <- intersect(Cv10_de, Cv30_de)

pdf("DE_venn_diagram.pdf", height = 3, width = 3)
draw.triple.venn(area1 = length(Cv5_de), area2 = length(Cv10_de), area3 = length(Cv30_de), n12 = length(de_5_10), n23 = length(de_10_30), n13 = length(de_5_30), 
                 n123 = length(de_all), category = c("5 ug/L", "10 ug/L", "30 ug/L"), lty = 1, reverse = T)
dev.off()


# # Calculate parametric means and confidence intervals for plotting
# cpm_counts <- cpm(data, prior.count = 2, log = T)
# means <- t(sapply(rownames(cpm_counts), function(i){tapply(cpm_counts[i,], groups$treatment, mean)}))
# sd <- t(sapply(rownames(cpm_counts), function(i){tapply(cpm_counts[i,], groups$treatment, sd)}))
# ci <- 1.96 * (sd/sqrt(4))
# uci <- means + ci
# lci <- means - ci


# Plot heatmaps of differential expression and candidate genes

GOI_all <- intersect(intersect(Cv5_de, Cv10_de), Cv30_de)
GOI_10_and_30 <- setdiff(intersect(Cv10_de, Cv30_de), GOI_all)

my.dist <- function(x) {
    as.dist(1-cor(t(x)))
}

cpm_counts <- cpm(data, prior.count = 2, log = T)
cu_anno <- HeatmapAnnotation(df = data.frame(Cu = rep(c(0,5,10,30), each = 4)), col = list(Cu = colorRamp2(c(0, 30), c("grey80", "grey10"))))
col <- colorRamp2(c(-10, 0, 10), c("#5e3c99", "#f7f7f7", "#e66101"))
hm1 <- Heatmap(cpm_counts[c(GOI_10_and_30, GOI_all),], 
               split = rep(c("A", "B"), c(length(GOI_10_and_30), length(GOI_all))),
               gap = unit(5, "mm"),
               col = col, 
               clustering_distance_rows = my.dist, 
               cluster_columns = F, 
               show_row_names = F, 
               top_annotation = cu_anno)

pdf("heatmap_de_genes.pdf", width = 4, height = 6)
hm1
dev.off()

# Heatmap for each concentration

hm5 <- Heatmap(cpm_counts[setdiff(Cv5_de, GOI_all),], 
               col = col, 
               clustering_distance_rows = my.dist, 
               cluster_columns = F, 
               show_row_names = F, 
               top_annotation = cu_anno)

pdf("heatmap_de_genes_cu5.pdf", width = 4, height = 6)
hm5
dev.off()

hm10 <- Heatmap(cpm_counts[setdiff(Cv10_de, GOI_10_and_30),], 
               col = col, 
               clustering_distance_rows = my.dist, 
               cluster_columns = F, 
               show_row_names = F, 
               top_annotation = cu_anno)

pdf("heatmap_de_genes_cu10.pdf", width = 4, height = 6)
hm10
dev.off()

hm30 <- Heatmap(cpm_counts[setdiff(Cv30_de, GOI_10_and_30),], 
               col = col, 
               clustering_distance_rows = my.dist, 
               cluster_columns = F, 
               show_row_names = F, 
               top_annotation = cu_anno)

pdf("heatmap_de_genes_cu30.pdf", width = 4, height = 6)
hm30
dev.off()


# Heatmap of candidate genes 

candidate_ids <- read.table("lamprey_candidate_gene_ids.txt", header = T)
c_vec <- which(row.names(data) %in% as.character(candidate_ids$trans_id))

cand_anno <- HeatmapAnnotation(df = data.frame(c05 = as.numeric(Cv5.lrt$table$PValue[c_vec] < 0.05),
                                               c10 = as.numeric(Cv10.lrt$table$PValue[c_vec] < 0.05),
                                               c30 = as.numeric(Cv30.lrt$table$PValue[c_vec] < 0.05)), 
                               col = list(c05 = c("1" = "black", "0" = "grey"),
                                          c10 = c("1" = "black", "0" = "grey"),
                                          c30 = c("1" = "black", "0" = "grey")),
                               which = "row")

hm2 <- Heatmap(cpm_counts[c_vec,], 
               col = col, 
               #clustering_distance_rows = my.dist, 
               cluster_columns = F, 
               show_row_names = T, 
               top_annotation = cu_anno
)

pdf("heatmap_candidate_genes.pdf", width = 6, height = 6)
hm2 + cand_anno
dev.off()


#### goseq annotation enrichment tests ####

if (!file.exists("goanno")) {
    anno <- read.table("lamprey_go_anno.txt", header = F, stringsAsFactors = F)
    colnames(anno) <- c("trans_id", "GO_term")
    PM.anno <- split(as.character(anno$GO_term), anno$trans_id, drop = T)
    
    # List GO DAG to facilitate extracting ancestors
    BP=as.list(GOBPANCESTOR)
    
    # Eactract ancestors for all anotated GO terms for all annotated genes
    PM.go.all = list()[names(PM.anno)]
    names(PM.go.all) = names(PM.anno)
    keys <- keys(GO.db)
    for(i in names(PM.anno)) {
        go_terms = PM.anno[[i]]
        for(j in go_terms) {
            if(j %in% keys && suppressMessages(select(GO.db, keys = j, columns = "ONTOLOGY")[,2] == "BP")) {
                PM.go.all[[i]] = append(PM.go.all[[i]], c(BP[[j]]))
            }
        }
        PM.go.all[[i]] = grep("GO:", PM.go.all[[i]], value = T)
        PM.go.all[[i]] = unique(PM.go.all[[i]])
    }
    saveRDS(PM.go.all, "goanno")
} else {
    PM.go.all <- readRDS("goanno")
}

lengths <- read.table("lamprey_gene_lengths.txt", header = F, row.names = 1)

Cv5_de_bool <- (p.adjust(Cv5.lrt$table$PValue, method = "BH")<0.05)
names(Cv5_de_bool) <- rownames(Cv5.lrt$table)
Cv5_pwf <- nullp(Cv5_de_bool, bias.data = lengths[,1])
goCv5 <- goseq(Cv5_pwf, gene2cat = PM.go.all)
goCv5 <- goCv5[which(goCv5$ontology == "BP" & goCv5$numInCat > 5),]
goCv5$over_represented_FDR <- p.adjust(goCv5$over_represented_pvalue, method = "BH")
write.table(goCv5[which(goCv5$over_represented_pvalue < 0.05),], "GO_over_represented_Cu05.txt", quote = F, sep = "\t", row.names = F)

Cv10_de_bool <- (p.adjust(Cv10.lrt$table$PValue, method = "BH")<0.05)
names(Cv10_de_bool) <- rownames(Cv10.lrt$table)
Cv10_pwf <- nullp(Cv10_de_bool, bias.data = lengths[,1])
goCv10 <- goseq(Cv10_pwf, gene2cat = PM.go.all)
goCv10 <- goCv10[which(goCv10$ontology == "BP" & goCv10$numInCat > 5),]
goCv10$over_represented_FDR <- p.adjust(goCv10$over_represented_pvalue, method = "BH")
write.table(goCv10[which(goCv10$over_represented_pvalue < 0.05),], "GO_over_represented_Cu10.txt", quote = F, sep = "\t", row.names = F)

Cv30_de_bool <- (p.adjust(Cv30.lrt$table$PValue, method = "BH")<0.05)
names(Cv30_de_bool) <- rownames(Cv30.lrt$table)
Cv30_pwf <- nullp(Cv30_de_bool, bias.data = lengths[,1])
goCv30 <- goseq(Cv30_pwf, gene2cat = PM.go.all)
goCv30 <- goCv30[which(goCv30$ontology == "BP" & goCv30$numInCat > 5),]
goCv30$over_represented_FDR <- p.adjust(goCv30$over_represented_pvalue, method = "BH")
write.table(goCv30[which(goCv30$over_represented_pvalue < 0.05),], "GO_over_represented_Cu30.txt", quote = F, sep = "\t", row.names = F)


Cv10_30_de_bool <- (p.adjust(Cv10.lrt$table$PValue, method = "BH")<0.05) * (p.adjust(Cv30.lrt$table$PValue, method = "BH")<0.05)
names(Cv10_30_de_bool) <- rownames(Cv30.lrt$table)
Cv10_30_pwf <- nullp(Cv10_30_de_bool, bias.data = lengths[,1])
goCv10_30 <- goseq(Cv10_30_pwf, gene2cat = PM.go.all)
goCv10_30 <- goCv10_30[which(goCv10_30$ontology == "BP" & goCv10_30$numInCat > 5),]
goCv10_30$over_represented_FDR <- p.adjust(goCv10_30$over_represented_pvalue, method = "BH")
write.table(goCv10_30[which(goCv10_30$over_represented_pvalue < 0.05),], "GO_over_represented_Cu10_30.txt", quote = F, sep = "\t", row.names = F)

CvAll_de_bool <- (p.adjust(Cv5.lrt$table$PValue, method = "BH")<0.05) * (p.adjust(Cv10.lrt$table$PValue, method = "BH")<0.05) * (p.adjust(Cv30.lrt$table$PValue, method = "BH")<0.05)
names(CvAll_de_bool) <- rownames(Cv30.lrt$table)
CvAll_pwf <- nullp(CvAll_de_bool, bias.data = lengths[,1])
goCvAll <- goseq(CvAll_pwf, gene2cat = PM.go.all)
goCvAll <- goCvAll[which(goCvAll$ontology == "BP" & goCvAll$numInCat > 5),]
goCvAll$over_represented_FDR <- p.adjust(goCvAll$over_represented_pvalue, method = "BH")
write.table(goCvAll[which(goCvAll$over_represented_pvalue < 0.05),], "GO_over_represented_all.txt", quote = F, sep = "\t", row.names = F)

## Compare GO over-representation overlap betwen doses
all_go_terms <- unique(as.character(unlist(PM.go.all)))
    go_overlap <- data.frame(Cv5 = all_go_terms %in% as.character(goCv5[goCv5$over_represented_FDR <= 0.05, "category"]),
                             Cv10 = all_go_terms %in% as.character(goCv10[goCv10$over_represented_FDR <= 0.05, "category"]),
                             Cv30 = all_go_terms %in% as.character(goCv30[goCv30$over_represented_FDR <= 0.05, "category"])
    )

vc <- vennCounts(go_overlap)

pdf("GO_venn_diagram.pdf", height = 3, width = 3)
vennDiagram(vc)
dev.off()



# # top logFC
# 
# Cv5_topFC_bool <- (Cv5.lrt$table$logFC < quantile(Cv5.lrt$table$logFC, 0.025)) | (Cv5.lrt$table$logFC > quantile(Cv5.lrt$table$logFC, 0.975))
# names(Cv5_topFC_bool) <- rownames(Cv5.lrt$table)
# Cv5_pwf <- nullp(Cv5_topFC_bool, bias.data = lengths[,1])
# goCv5 <- goseq(Cv5_pwf, gene2cat = PM.go.all)
# goCv5 <- goCv5[which(goCv5$ontology == "BP" & goCv5$numInCat > 5),]
# goCv5$over_represented_FDR <- p.adjust(goCv5$over_represented_pvalue, method = "BH")
# write.table(goCv5[which(goCv5$over_represented_pvalue < 0.05),], "GO_over_represented_Cu05_topFC.txt", quote = F, sep = "\t")
# 
# Cv10_topFC_bool <- (Cv10.lrt$table$logFC < quantile(Cv10.lrt$table$logFC, 0.025)) | (Cv10.lrt$table$logFC > quantile(Cv10.lrt$table$logFC, 0.975))
# names(Cv10_topFC_bool) <- rownames(Cv10.lrt$table)
# Cv10_pwf <- nullp(Cv10_topFC_bool, bias.data = lengths[,1])
# goCv10 <- goseq(Cv10_pwf, gene2cat = PM.go.all)
# goCv10 <- goCv10[which(goCv10$ontology == "BP" & goCv10$numInCat > 5),]
# goCv10$over_represented_FDR <- p.adjust(goCv10$over_represented_pvalue, method = "BH")
# write.table(goCv10[which(goCv10$over_represented_pvalue < 0.05),], "GO_over_represented_Cu10_topFC.txt", quote = F, sep = "\t")
# 
# Cv30_topFC_bool <- (Cv30.lrt$table$logFC < quantile(Cv30.lrt$table$logFC, 0.025)) | (Cv30.lrt$table$logFC > quantile(Cv30.lrt$table$logFC, 0.975))
# names(Cv30_topFC_bool) <- rownames(Cv30.lrt$table)
# Cv30_pwf <- nullp(Cv30_topFC_bool, bias.data = lengths[,1])
# goCv30 <- goseq(Cv30_pwf, gene2cat = PM.go.all)
# goCv30 <- goCv30[which(goCv30$ontology == "BP" & goCv30$numInCat > 5),]
# goCv30$over_represented_FDR <- p.adjust(goCv30$over_represented_pvalue, method = "BH")
# write.table(goCv30[which(goCv30$over_represented_pvalue < 0.05),], "GO_over_represented_Cu30_topFC.txt", quote = F, sep = "\t")

