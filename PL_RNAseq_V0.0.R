#load input data
if (length(commandArgs(trailingOnly=TRUE))>0) {
  args <- commandArgs(trailingOnly=TRUE)
}
# library(phangorn)
print(length(args))
if (length(args) < 3) {
  stop(" Usage: bloop", call.=FALSE)
}

cat("You ran the program with ", args,"\n")

raw_dir <- (as.character(args[2]))
wor_dir <- (as.character(args[1]))
strain <- (as.character(args[5]))
setwd(raw_dir)
dir.create("output_plots")
dir.create("output_tables")
dir.create("pathview")
setwd(wor_dir)
input_table_DESEQ2 <- read.table(as.character(args[3]), header=T, sep="\t")

padj.cutoff <- 0.05
lfc.cutoff <- 0.58

library(DESeq2)
library(apeglm)
library(pheatmap)
library(ggplot2)
library(pathview)
library(clusterProfiler)
library(EnhancedVolcano)
library(GOfuncR)
library(tibble)
library(dplyr)
library(BiocGenerics)
library(data.table)
library(Biobase)
library(MatrixGenerics)
library(S4Vectors)
library(ggrepel)

rownames(input_table_DESEQ2) <- input_table_DESEQ2$SampleName

head(input_table_DESEQ2)
nrow(input_table_DESEQ2)
ncol(input_table_DESEQ2)
nmcolz <- colnames(input_table_DESEQ2)


if (length(unique(input_table_DESEQ2$Gene)) > 1 && (length(unique(input_table_DESEQ2$Condition)) > 1)){
  input_star <- DESeqDataSetFromHTSeqCount(sampleTable = input_table_DESEQ2, directory = (args[4]), design =  formula(paste("~" ,nmcolz[4:4],"+",nmcolz[3:3],"+",nmcolz[4:4],":",nmcolz[3:3])))
}else {
  input_star <- DESeqDataSetFromHTSeqCount(sampleTable = input_table_DESEQ2, directory = (args[4]), design = formula(paste("~",nmcolz[4:4])))
}


nrow(input_star)

input_star <- input_star[rowSums(counts(input_star)) > 10, ]

nrow (input_star)

input_star2 <- DESeq(input_star, betaPrior = FALSE)

print("results names")
resultsNames(input_star2)

DESEQ2_norm_counts <- log2(counts(input_star2, normalized = TRUE)+1)

head(DESEQ2_norm_counts)

setwd(raw_dir)

write.table(DESEQ2_norm_counts, "output_tables/DESEQ2_norm_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")

setwd(wor_dir)

resultsNames(input_star2)

extractrial <- unique(input_table_DESEQ2[,4])
extractrial2 <- unique(input_table_DESEQ2[,3])

paste(nmcolz[4:4],"_", extractrial[1:1],"_vs_",extractrial[2:2], sep = "")
paste(nmcolz[4:4],"_", extractrial2[1:1],"_vs_",extractrial2[2:2], sep = "")

if (length(unique(input_table_DESEQ2$Gene)) > 1 && (length(unique(input_table_DESEQ2$Condition)) > 1)){
  DESEQ2_DEG <- results(object = input_star2, contrast= paste(c(nmcolz[4:4], extractrial[1:1], extractrial[2:2])))
  DESEQ2_DEG_alt <- results(object = input_star2, contrast= paste(c(nmcolz[3:3], extractrial2[1:1], extractrial2[2:2])))
  head(DESEQ2_DEG)
  head(DESEQ2_DEG_alt)
  setwd(raw_dir)
  write.table(DESEQ2_DEG, "output_tables/DESEQ2_DEG.txt", quote=F, col.names=T, row.names=T, sep="\t")
  write.table(DESEQ2_DEG_alt, "output_tables/DESEQ2_DEG_alt.txt", quote=F, col.names=T, row.names=T, sep="\t")
}else {
  DESEQ2_DEG <- results(object = input_star2, contrast= paste(c(nmcolz[4:4], extractrial[1:1], extractrial[2:2])))
  head(DESEQ2_DEG)
  setwd(raw_dir)
  write.table(DESEQ2_DEG, "output_tables/DESEQ2_DEG.txt", quote=F, col.names=T, row.names=T, sep="\t")
}

setwd(wor_dir)

DESEQ2_var_stabl <- vst(input_star2)

DESEQ2_DistMatrix <- as.matrix(dist(t(assay(DESEQ2_var_stabl))))

setwd(raw_dir)

pdf("output_plots/sample_distance_heatmap_DESEQ.pdf")
pheatmap(DESEQ2_DistMatrix)

dev.off()

pdf("output_plots/gene_condition1.pdf")
plotCounts(input_star2, gene=which.min(DESEQ2_DEG$padj), intgroup=c(nmcolz[4:4], nmcolz[3:3]))

dev.off()


select <- order(rowMeans(counts(input_star2,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(input_star2)[,c(nmcolz[4:4], nmcolz[3:3])])

pdf("output_plots/gene_heatmap.pdf", height = 20, width = 20)
pheatmap(assay(DESEQ2_var_stabl)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

pdf("output_plots/gene_heatmap_clustered.pdf", height = 20, width = 20)
pheatmap::pheatmap(assay(DESEQ2_var_stabl)[select,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)
dev.off()


if (length(unique(input_table_DESEQ2$Gene)) > 1 && (length(unique(input_table_DESEQ2$Condition)) > 1)){
  p <- plotPCA(object = DESEQ2_var_stabl, intgroup = c(nmcolz[3:3], nmcolz[4:4]))
  p <- p + geom_label_repel(aes(label =DESEQ2_var_stabl@colData@rownames), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + theme_bw()
  ggsave(filename = "output_plots/DESEQ_conditions_all_PCA.pdf", height = 20, width = 20)
}else {
  p <-plotPCA(object = DESEQ2_var_stabl, intgroup = nmcolz[4:4])
  p <- p + geom_label_repel(aes(label =DESEQ2_var_stabl@colData@rownames), box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50') + theme_bw()
  ggsave(filename = "output_plots/DESEQ_conditions_PCA.pdf", height = 20, width = 20)
}


if (length(unique(input_table_DESEQ2$Gene)) > 1 && (length(unique(input_table_DESEQ2$Condition)) > 1)){
  EnhancedVolcano(DESEQ2_DEG, lab = rownames(DESEQ2_DEG), x='log2FoldChange', y='padj', title = 'volcano plot', pCutoff = 0.01, FCcutoff = 1.5)
  ggsave(filename = "output_plots/volcano_plot.pdf")
  EnhancedVolcano(DESEQ2_DEG_alt, lab = rownames(DESEQ2_DEG_alt), x='log2FoldChange', y='padj', title = 'volcano plot', pCutoff = 0.01, FCcutoff = 1.5)
  ggsave(filename = "output_plots/volcano_plot_alt_compare.pdf")
}else {
  EnhancedVolcano(DESEQ2_DEG, lab = rownames(DESEQ2_DEG), x='log2FoldChange', y='padj', title = 'volcano plot', pCutoff = 0.01, FCcutoff = 1.5)
  ggsave(filename = "output_plots/volcano_plot.pdf")
}

setwd(wor_dir)

DESEQ2_table_dge_tb <- DESEQ2_DEG %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  as_tibble()


DESEQ2_DEGX <- DESEQ2_table_dge_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

setwd(raw_dir)

if (strain == "PA14"){
  PA_info_genes <- data.table::fread("~/Desktop/PA_info_genes.txt")
  merge_attempt <- merge(DESEQ2_table_dge_tb, PA_info_genes, by.x='gene', by.y='PA14_num')
  head(merge_attempt)
}else {
  PA_info_genes <- data.table::fread("~/Desktop/PA_info_genes.txt")
  merge_attempt <- merge(DESEQ2_table_dge_tb, PA_info_genes, by.x='gene', by.y='PA_num')
  head(merge_attempt)
}

pdf("output_plots/volcano_plot_named.pdf")
EnhancedVolcano(merge_attempt, lab = merge_attempt$col_use, x='log2FoldChange', y='padj', title = 'volcano plot', pCutoff = 0.01, FCcutoff = 1.5)

dev.off()

setwd(wor_dir)

sig_genes_pth <- c(DESEQ2_DEGX$gene) #take names of significantly differentially expressed genes
sig_scores_pth <- c(DESEQ2_DEGX$padj) #take their associated adjusted pval
sig_fc_pth <- c(DESEQ2_DEGX$log2FoldChange) #tack on the log2 fold changes

gene_ratio_input <- data.frame(gene_ids = c(sig_genes_pth), gene_scores = c(sig_scores_pth), log2FC = as.numeric(c(sig_fc_pth))) #merge into format required for downstream

geneList <- gene_ratio_input[,2] #list out genes
names(geneList) <- as.character(gene_ratio_input[,1]) #file formatting for downstream steps

gene_ratio_inputUP <- subset(gene_ratio_input, log2FC > 0)
geneListUP <- gene_ratio_inputUP[,3]

gene_ratio_inputDOWN <- subset(gene_ratio_input, log2FC < 0)
geneListDOWN <- gene_ratio_inputDOWN[,3]

names(geneListUP) <- as.character(gene_ratio_inputUP[,1])
names(geneListDOWN) <- as.character(gene_ratio_inputDOWN[,1])

geneList <- sort(geneList, decreasing=TRUE) #sorting geneList relative to pval

geneListUP <- sort(geneListUP, decreasing = TRUE)
geneListDOWN <- sort(geneListDOWN, decreasing = TRUE)

gene_up <- names(geneListUP)
gene_down <- names(geneListDOWN)
gene <- names(geneList)[abs(geneList)<0.05] #extracting names of genes with pval less than 0.05


#import GO annotations from pseudomonas.com gene ontology - need to find a way to get this to be uniform to all gene ontologies, but currently works well

if (strain == "PA14"){
  PAO1_GO_all<-data.table::fread("~/Desktop/PA14_gene_ontology_csv.csv") #don't hardcode table - let user import
}else {
  PAO1_GO_all<-data.table::fread("~/Desktop/gene_ontology_csv.csv") #don't hardcode table - let user import
  PAO1_custom_list <- data.table::fread("~/Desktop/R_Annotations_Daniel_20190409.csv")
  custom_term2gene <- subset(PAO1_custom_list, select=c(2,1))
  custom_term2name <- subset(PAO1_custom_list, select=c(2,3,1))
  cluster_profiler_custom <- enricher(gene, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = custom_term2gene, TERM2NAME = custom_term2name)
  cluster_profiler_customUP <- enricher(gene_up, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = custom_term2gene, TERM2NAME = custom_term2name)
  cluster_profiler_customDOWN <- enricher(gene_down, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = custom_term2gene, TERM2NAME = custom_term2name)

  if (is.na(cluster_profiler_custom[1,5]) == "TRUE"){
    print("no significant enrichment terms")
  }else {
    setwd(raw_dir)
    cluster_profiler_custom@result$cp_GeneRatio <- cluster_profiler_custom@result$GeneRatio
    cluster_profiler_custom@result$cp_Description <- cluster_profiler_custom@result$Description
    cluster_profiler_custom@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_custom@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_custom@result$BgRatio, "/"), `[[`, 1), sep = "")
    cluster_profiler_custom@result$Description <- paste(sapply(strsplit(cluster_profiler_custom@result$ID, "\t"), `[[`, 1), sep = "")
    dotplot(cluster_profiler_custom, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("Custom gene list enrichment")
    ggsave(filename = "output_plots/GO_gene_ratio_custom.pdf")
    setwd(wor_dir)
  }
  if (is.na(cluster_profiler_customUP[1,5]) == "TRUE"){
    print("no significant enrichment terms")
  }else {
    setwd(raw_dir)
    cluster_profiler_customUP@result$cp_GeneRatio <- cluster_profiler_customUP@result$GeneRatio
    cluster_profiler_customUP@result$cp_Description <- cluster_profiler_customUP@result$Description
    cluster_profiler_customUP@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_customUP@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_customUP@result$BgRatio, "/"), `[[`, 1), sep = "")
    cluster_profiler_customUP@result$Description <- paste(sapply(strsplit(cluster_profiler_customUP@result$ID, "\t"), `[[`, 1), sep = "")
    dotplot(cluster_profiler_customUP, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("Custom gene list enrichment (up-regulated)")
    ggsave(filename = "output_plots/GO_gene_ratio_custom_upregulated.pdf")
    setwd(wor_dir)
  }
  if (is.na(cluster_profiler_customDOWN[1,5]) == "TRUE"){
    print("no significant enrichment terms")
  }else {
    setwd(raw_dir)
    cluster_profiler_customDOWN@result$cp_GeneRatio <- cluster_profiler_customDOWN@result$GeneRatio
    cluster_profiler_customDOWN@result$cp_Description <- cluster_profiler_customDOWN@result$Description
    cluster_profiler_customDOWN@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_customDOWN@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_customDOWN@result$BgRatio, "/"), `[[`, 1), sep = "")
    cluster_profiler_customDOWN@result$Description <- paste(sapply(strsplit(cluster_profiler_customDOWN@result$ID, "\t"), `[[`, 1), sep = "")
    dotplot(cluster_profiler_customDOWN, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("Custom gene list enrichment (down-regulated)")
    ggsave(filename = "output_plots/GO_gene_ratio_custom_downregulated.pdf")
    setwd(wor_dir)
  }
}


term2name <- subset(PAO1_GO_all, select=c(5,6,7)) #subsetting the relevant columns
term2gene <- subset(PAO1_GO_all, select=c(5,1)) #extracting relevant geneA = GO:X for next step)
term2name <- subset(PAO1_GO_all, select=c(5,6,7)) #extracting relevant GO# = function = category eg: GO:0005524 _ ATP binding _ molecular_function

#enrichment analysis on all GO

cluster_profiler_enriched <- enricher(gene, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene, TERM2NAME = term2name) #gene ontology enrichment analysis, bonferoni correction, all GOs)
cluster_profiler_enrichedUP <- enricher(gene_up, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene, TERM2NAME = term2name) #gene ontology enrichment analysis, bonferoni correction, all GOs)
cluster_profiler_enrichedDOWN <- enricher(gene_down, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene, TERM2NAME = term2name) #gene ontology enrichment analysis, bonferoni correction, all GOs)



if (is.na(cluster_profiler_enriched[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  }else {
  setwd(raw_dir)
  cluster_profiler_enriched@result$cp_GeneRatio <- cluster_profiler_enriched@result$GeneRatio
  cluster_profiler_enriched@result$cp_Description <- cluster_profiler_enriched@result$Description
  cluster_profiler_enriched@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enriched@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enriched, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (all ontologies)")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched.pdf")
  setwd(wor_dir)
  }
if (is.na(cluster_profiler_enrichedUP[1,5]) == "TRUE"){
  print("no significant enrichment terms")
}else {
  setwd(raw_dir)
  cluster_profiler_enrichedUP@result$cp_GeneRatio <- cluster_profiler_enrichedUP@result$GeneRatio
  cluster_profiler_enrichedUP@result$cp_Description <- cluster_profiler_enrichedUP@result$Description
  cluster_profiler_enrichedUP@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enrichedUP@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enrichedUP@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enrichedUP@result$Description <- paste(sapply(strsplit(cluster_profiler_enrichedUP@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enrichedUP@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enrichedUP, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (all ontologies) up-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_upregulated.pdf")
  setwd(wor_dir)
}
if (is.na(cluster_profiler_enrichedDOWN[1,5]) == "TRUE"){
  print("no significant enrichment terms")
}else {
  setwd(raw_dir)
  cluster_profiler_enrichedDOWN@result$cp_GeneRatio <- cluster_profiler_enrichedDOWN@result$GeneRatio
  cluster_profiler_enrichedDOWN@result$cp_Description <- cluster_profiler_enrichedDOWN@result$Description
  cluster_profiler_enrichedDOWN@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enrichedDOWN@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enrichedDOWN@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enrichedDOWN@result$Description <- paste(sapply(strsplit(cluster_profiler_enrichedDOWN@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enrichedDOWN@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enrichedDOWN, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (all ontologies) down-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_downregulated.pdf")
  setwd(wor_dir)
}
#However, doing an enrichment analysis on all GO is silly, you can do them on specific categories this example is molecular_function (fortunately all of this info is in the pseudomonas.com ontology sheet)

PAO1_GO_MF <- PAO1_GO_all[ which(PAO1_GO_all$Namespace=='molecular_function')] #subsetting all ontologies belonging to molecular function
term2name_MF <- subset(PAO1_GO_MF, select=c(5,6,7)) #same as above but on subsetted data
term2gene_MF <- subset(PAO1_GO_MF, select=c(5,1)) # "   "   "   "

#term2gene_MF_aschar <- as.character(term2gene_MF$Accession) #To make sure your GOs are ranking properly you can run this part here, it will show you the GO relationships - however, cluster_profiler actually uses the highest parent in the plots
#tmp <- get_parent_nodes(term2gene_MF_aschar)
#head(tmp)

cluster_profiler_enriched_MF <- enricher(gene, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene_MF, TERM2NAME = term2name_MF) #new GO analysis on subsetted data
cluster_profiler_enriched_MFUP <- enricher(gene_up, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene_MF, TERM2NAME = term2name_MF) #new GO analysis on subsetted data
cluster_profiler_enriched_MFDOWN <- enricher(gene_down, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene_MF, TERM2NAME = term2name_MF) #new GO analysis on subsetted data

if (is.na(cluster_profiler_enriched_MF[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  }else {
    setwd(raw_dir)
    cluster_profiler_enriched_MF@result$cp_GeneRatio <- cluster_profiler_enriched_MF@result$GeneRatio
    cluster_profiler_enriched_MF@result$cp_Description <- cluster_profiler_enriched_MF@result$Description
    cluster_profiler_enriched_MF@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_MF@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_MF@result$BgRatio, "/"), `[[`, 1), sep = "")
    cluster_profiler_enriched_MF@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_MF@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_MF@result$ID, "\t"), `[[`, 1), ")", sep = "")
    dotplot(cluster_profiler_enriched_MF, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Molecular Function)")
    ggsave(filename = "output_plots/GO_gene_ratio_enriched_MF.pdf")
    setwd(wor_dir)
  }
if (is.na(cluster_profiler_enriched_MFUP[1,5]) == "TRUE"){
  print("no significant enrichment terms")
}else {
  setwd(raw_dir)
  cluster_profiler_enriched_MFUP@result$cp_GeneRatio <- cluster_profiler_enriched_MFUP@result$GeneRatio
  cluster_profiler_enriched_MFUP@result$cp_Description <- cluster_profiler_enriched_MFUP@result$Description
  cluster_profiler_enriched_MFUP@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_MFUP@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_MFUP@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enriched_MFUP@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_MFUP@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_MFUP@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enriched_MFUP, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Molecular Function) up-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_MF_upregulated.pdf")
  setwd(wor_dir)
}
if (is.na(cluster_profiler_enriched_MFDOWN[1,5]) == "TRUE"){
  print("no significant enrichment terms")
}else {
  setwd(raw_dir)
  cluster_profiler_enriched_MFDOWN@result$cp_GeneRatio <- cluster_profiler_enriched_MFDOWN@result$GeneRatio
  cluster_profiler_enriched_MFDOWN@result$cp_Description <- cluster_profiler_enriched_MFDOWN@result$Description
  cluster_profiler_enriched_MFDOWN@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_MFDOWN@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_MFDOWN@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enriched_MFDOWN@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_MFDOWN@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_MFDOWN@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enriched_MFDOWN, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Molecular Function) down-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_MF_downregulated.pdf")
  setwd(wor_dir)
}

#BP only
PAO1_GO_BP <- PAO1_GO_all[ which(PAO1_GO_all$Namespace=='biological_process')] #same as above but with biological process - this will likely be the one you want to see.
term2name_BP <- subset(PAO1_GO_BP, select=c(5,6,7))
term2gene_BP <- subset(PAO1_GO_BP, select=c(5,1))
#term2gene_BP_aschar <- as.character(term2gene_BP$Accession)
#par_nodes_BP_PAO1 <- get_parent_nodes(term2gene_BP_aschar)
cluster_profiler_enriched_BP <- enricher(gene, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene_BP, TERM2NAME = term2name_BP)
cluster_profiler_enriched_BPUP <- enricher(gene_up, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene_BP, TERM2NAME = term2name_BP)
cluster_profiler_enriched_BPDOWN <- enricher(gene_down, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene_BP, TERM2NAME = term2name_BP)

if (is.na(cluster_profiler_enriched_BP[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  } else {
    setwd(raw_dir)
    cluster_profiler_enriched_BP@result$cp_GeneRatio <- cluster_profiler_enriched_BP@result$GeneRatio
    cluster_profiler_enriched_BP@result$cp_Description <- cluster_profiler_enriched_BP@result$Description
    cluster_profiler_enriched_BP@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_BP@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_BP@result$BgRatio, "/"), `[[`, 1), sep = "")
    cluster_profiler_enriched_BP@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_BP@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_BP@result$ID, "\t"), `[[`, 1), ")", sep = "")
    dotplot(cluster_profiler_enriched_BP, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Biological Process)")
    ggsave(filename = "output_plots/GO_gene_ratio_enriched_BP.pdf")
    setwd(wor_dir)
  }
if (is.na(cluster_profiler_enriched_BPUP[1,5]) == "TRUE"){
  print("no significant enrichment terms")
} else {
  setwd(raw_dir)
  cluster_profiler_enriched_BPUP@result$cp_GeneRatio <- cluster_profiler_enriched_BPUP@result$GeneRatio
  cluster_profiler_enriched_BPUP@result$cp_Description <- cluster_profiler_enriched_BPUP@result$Description
  cluster_profiler_enriched_BPUP@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_BPUP@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_BPUP@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enriched_BPUP@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_BPUP@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_BPUP@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enriched_BPUP, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Biological Process) up-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_BP_upregulated.pdf")
  setwd(wor_dir)
}
if (is.na(cluster_profiler_enriched_BPDOWN[1,5]) == "TRUE"){
  print("no significant enrichment terms")
} else {
  setwd(raw_dir)
  cluster_profiler_enriched_BPDOWN@result$cp_GeneRatio <- cluster_profiler_enriched_BPDOWN@result$GeneRatio
  cluster_profiler_enriched_BPDOWN@result$cp_Description <- cluster_profiler_enriched_BPDOWN@result$Description
  cluster_profiler_enriched_BPDOWN@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_BPDOWN@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_BPDOWN@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enriched_BPDOWN@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_BPDOWN@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_BPDOWN@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enriched_BPDOWN, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Biological Process) down-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_BP_downregulated.pdf")
  setwd(wor_dir)
}
#CC only

PAO1_GO_CC <- PAO1_GO_all[ which(PAO1_GO_all$Namespace=='cellular_component')] # and again with cellular components
term2name_CC <- subset(PAO1_GO_CC, select=c(5,6,7))
term2gene_CC <- subset(PAO1_GO_CC, select=c(5,1))
#term2gene_CC_aschar <- as.character(term2gene_CC$Accession)
#par_nodes_CC_PAO1 <- get_parent_nodes(term2gene_CC_aschar)

cluster_profiler_enriched_CC <- enricher(gene, qvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = term2gene_CC, TERM2NAME = term2name_CC)
cluster_profiler_enriched_CCUP <- enricher(gene_up, qvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = term2gene_CC, TERM2NAME = term2name_CC)
cluster_profiler_enriched_CCDOWN <- enricher(gene_down, qvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = term2gene_CC, TERM2NAME = term2name_CC)

if (is.na(cluster_profiler_enriched_CC[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  } else {
    setwd(raw_dir)
    cluster_profiler_enriched_CC@result$cp_GeneRatio <- cluster_profiler_enriched_CC@result$GeneRatio
    cluster_profiler_enriched_CC@result$cp_Description <- cluster_profiler_enriched_CC@result$Description
    cluster_profiler_enriched_CC@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_CC@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_CC@result$BgRatio, "/"), `[[`, 1), sep = "")
    cluster_profiler_enriched_CC@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_CC@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_CC@result$ID, "\t"), `[[`, 1), ")", sep = "")
    dotplot(cluster_profiler_enriched_CC, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Celluar Component)")
    ggsave(filename = "output_plots/GO_gene_ratio_enriched_CC.pdf")
    setwd(wor_dir)
  }
if (is.na(cluster_profiler_enriched_CCUP[1,5]) == "TRUE"){
  print("no significant enrichment terms")
} else {
  setwd(raw_dir)
  cluster_profiler_enriched_CCUP@result$cp_GeneRatio <- cluster_profiler_enriched_CCUP@result$GeneRatio
  cluster_profiler_enriched_CCUP@result$cp_Description <- cluster_profiler_enriched_CCUP@result$Description
  cluster_profiler_enriched_CCUP@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_CCUP@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_CCUP@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enriched_CCUP@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_CCUP@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_CCUP@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enriched_CCUP, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Celluar Component) up-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_CC_upregulated.pdf")
  setwd(wor_dir)
}
if (is.na(cluster_profiler_enriched_CCDOWN[1,5]) == "TRUE"){
  print("no significant enrichment terms")
} else {
  setwd(raw_dir)
  cluster_profiler_enriched_CCDOWN@result$cp_GeneRatio <- cluster_profiler_enriched_CCDOWN@result$GeneRatio
  cluster_profiler_enriched_CCDOWN@result$cp_Description <- cluster_profiler_enriched_CCDOWN@result$Description
  cluster_profiler_enriched_CCDOWN@result$GeneRatio <- paste(sapply(strsplit(cluster_profiler_enriched_CCDOWN@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(cluster_profiler_enriched_CCDOWN@result$BgRatio, "/"), `[[`, 1), sep = "")
  cluster_profiler_enriched_CCDOWN@result$Description <- paste(sapply(strsplit(cluster_profiler_enriched_CCDOWN@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(cluster_profiler_enriched_CCDOWN@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(cluster_profiler_enriched_CCDOWN, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("GO gene ratio (Celluar Component) down-regulated")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched_CC_downregulated.pdf")
  setwd(wor_dir)
}

#cluster_profiler can also easily use kegg, which is a massive load off

if (strain == "PA14"){
  clustprof_kegg <- enrichKEGG(gene = gene, organism = "pau", pvalueCutoff = 0.05)
  clustprof_keggUP <- enrichKEGG(gene = gene_up, organism = "pau", pvalueCutoff = 0.05)
  clustprof_keggDOWN <- enrichKEGG(gene = gene_down, organism = "pau", pvalueCutoff = 0.05)
}else {
  clustprof_kegg <- enrichKEGG(gene = gene, organism = "pae", pvalueCutoff = 0.05)
  clustprof_keggUP <- enrichKEGG(gene = gene_up, organism = "pae", pvalueCutoff = 0.05)
  clustprof_keggDOWN <- enrichKEGG(gene = gene_down, organism = "pae", pvalueCutoff = 0.05)
}

if (is.na(clustprof_kegg[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  } else {
  setwd(raw_dir)
  clustprof_kegg@result$cp_GeneRatio <- clustprof_kegg@result$GeneRatio
  clustprof_kegg@result$cp_Description <- clustprof_kegg@result$Description
  clustprof_kegg@result$GeneRatio <- paste(sapply(strsplit(clustprof_kegg@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(clustprof_kegg@result$BgRatio, "/"), `[[`, 1), sep = "")
  clustprof_kegg@result$Description <- paste(sapply(strsplit(clustprof_kegg@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(clustprof_kegg@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(clustprof_kegg, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("Gene ratio (KEGG pathways)")
  ggsave(filename = "output_plots/KEGG_gene_ratio_enriched.pdf")
  setwd(wor_dir)
  }
if (is.na(clustprof_keggUP[1,5]) == "TRUE"){
  print("no significant enrichment terms")
} else {
  setwd(raw_dir)
  clustprof_keggUP@result$cp_GeneRatio <- clustprof_keggUP@result$GeneRatio
  clustprof_keggUP@result$cp_Description <- clustprof_keggUP@result$Description
  clustprof_keggUP@result$GeneRatio <- paste(sapply(strsplit(clustprof_keggUP@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(clustprof_keggUP@result$BgRatio, "/"), `[[`, 1), sep = "")
  clustprof_keggUP@result$Description <- paste(sapply(strsplit(clustprof_keggUP@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(clustprof_keggUP@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(clustprof_keggUP, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("Gene ratio (KEGG pathways) up-regulated")
  ggsave(filename = "output_plots/KEGG_gene_ratio_enriched_upregulated.pdf")
  setwd(wor_dir)
}
if (is.na(clustprof_keggDOWN[1,5]) == "TRUE"){
  print("no significant enrichment terms")
} else {
  setwd(raw_dir)
  clustprof_keggDOWN@result$cp_GeneRatio <- clustprof_keggDOWN@result$GeneRatio
  clustprof_keggDOWN@result$cp_Description <- clustprof_keggDOWN@result$Description
  clustprof_keggDOWN@result$GeneRatio <- paste(sapply(strsplit(clustprof_keggDOWN@result$GeneRatio, "/"), `[[`, 1), "/", sapply(strsplit(clustprof_keggDOWN@result$BgRatio, "/"), `[[`, 1), sep = "")
  clustprof_keggDOWN@result$Description <- paste(sapply(strsplit(clustprof_keggDOWN@result$Description, "/t"), `[[`, 1), " ","\n", "(", sapply(strsplit(clustprof_keggDOWN@result$ID, "\t"), `[[`, 1), ")", sep = "")
  dotplot(clustprof_keggDOWN, x = "GeneRatio", orderBy = "x", showCategory=20) + ggtitle("Gene ratio (KEGG pathways) down-regulated")
  ggsave(filename = "output_plots/KEGG_gene_ratio_enriched_downregulated.pdf")
  setwd(wor_dir)
}


data_fold_changes <- DESEQ2_DEG$log2FoldChange
names(data_fold_changes) <- rownames(DESEQ2_DEG)

setwd(raw_dir)
setwd("pathview")

if (strain == "PA14"){
  tmp <- pathview(gene.data = data_fold_changes, pathway.id = clustprof_kegg[,1], species = "pau", gene.idtype = "kegg", limit = list(gene= 2)) #this will output all of the differentially regulated pathways from pathview showing fancy pathway maps
}else {
  tmp <- pathview(gene.data = data_fold_changes, pathway.id = clustprof_kegg[,1], species = "pae", gene.idtype = "kegg", limit = list(gene= 2)) #this will output all of the differentially regulated pathways from pathview showing fancy pathway maps
}

setwd(wor_dir)

#quit(save="no")
