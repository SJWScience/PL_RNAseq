#load input data
args = commandArgs(trailingOnly=TRUE)
args

raw_dir <- (as.character(args[2]))
wor_dir <- (as.character(args[1]))
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

rownames(input_table_DESEQ2) <- input_table_DESEQ2$SampleName

head(input_table_DESEQ2)
nrow(input_table_DESEQ2)
ncol(input_table_DESEQ2)
nmcolz <- colnames(input_table_DESEQ2)


if (length(unique(input_table_DESEQ2$Gene)) > 1 && (length(unique(input_table_DESEQ2$Condition)) > 1)){
  input_star <- DESeqDataSetFromHTSeqCount(sampleTable = input_table_DESEQ2, directory = (args[4]), design =  formula(paste("~" ,nmcolz[4:4],"+",nmcolz[3:3])))
}else {
  input_star <- DESeqDataSetFromHTSeqCount(sampleTable = input_table_DESEQ2, directory = (args[4]), design = formula(paste("~",nmcolz[4:4])))
}


nrow(input_star)

input_star <- input_star[rowSums(counts(input_star)) > 10, ]

nrow (input_star)

input_star2 <- DESeq(input_star, betaPrior = TRUE)

DESEQ2_norm_counts <- log2(counts(input_star2, normalized = TRUE)+1)

head(DESEQ2_norm_counts)

setwd(raw_dir)

write.table(DESEQ2_norm_counts, "output_tables/DESEQ2_norm_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")

setwd(wor_dir)

resultsNames(input_star2)

extractrial <- unique(input_table_DESEQ2[,4])
paste(nmcolz[4:4],"_", extractrial[1:1],"_vs_",extractrial[2:2], sep = "")

DESEQ2_DEG <- results(object = input_star2, contrast= paste(c(nmcolz[4:4], extractrial[1:1], extractrial[2:2])))
##DESEQ2_DEG_alt <- results(object = input_star2, name="Gene_spoT_vs_relA")
##this part needs to not be hardcoded##
#DESEQ2_DEG <- results(object = input_star2, contrast = c("Condition","In","Un"))
#paste(nmcolz[4:4], extractrial[1:1],extractrial[2:2], sep = ",")
head(DESEQ2_DEG)
##head(DESEQ2_DEG_alt)
#DESEQ2_DEG_shrink <- lfcShrink(dds = input_star2, coef= paste(nmcolz[4:4],"_", extractrial[2:2],"_vs_",extractrial[1:1], sep = ""), type="apeglm") #needs to not be hardcoded
#head(DESEQ2_DEG_shrink)

##DESEQ2_alt_DEG_shrink <- lfcShrink(dds = input_star2, coef="Gene_spoT_vs_relA", type="apeglm")
##head(DESEQ2_alt_DEG_shrink)

setwd(raw_dir)

#write.table(DESEQ2_DEG_shrink, "output_tables/DESEQ2_DEG_shrink.txt", quote=F, col.names=T, row.names=T, sep="\t")
##write.table(DESEQ2_alt_DEG_shrink, "DESEQ2_alt_DEG_shrink.txt", quote=F, col.names=T, row.names=T, sep="\t")

write.table(DESEQ2_DEG, "output_tables/DESEQ2_DEG.txt", quote=F, col.names=T, row.names=T, sep="\t")
#write.table(DESEQ2_DEG_alt, "DESEQ2_DEG_alt.txt", quote=F, col.names=T, row.names=T, sep="\t")

setwd(wor_dir)

DESEQ2_var_stabl <- vst(input_star2)

DESEQ2_DistMatrix <- as.matrix(dist(t(assay(DESEQ2_var_stabl))))

setwd(raw_dir)

pdf("output_plots/sample_distance_heatmap_DESEQ.pdf")
pheatmap(DESEQ2_DistMatrix)

dev.off()

pdf("output_plots/DESEQ_gene_PCA.pdf")
plotPCA(object = DESEQ2_var_stabl, intgroup = nmcolz[3:3])

dev.off()

pdf("output_plots/DESEQ_PCA.pdf")
plotPCA(object = DESEQ2_var_stabl, intgroup = nmcolz[4:4])

dev.off()

pdf("output_plots/volcano_plot.pdf")
EnhancedVolcano(DESEQ2_DEG, lab = rownames(DESEQ2_DEG), x='log2FoldChange', y='padj', title = 'volcano plot', pCutoff = 0.01, FCcutoff = 1.5)

dev.off()

setwd(wor_dir)

DESEQ2_table_dge_tb <- DESEQ2_DEG %>%
  data.frame() %>%
  tibble::rownames_to_column(var="gene") %>%
  as_tibble()


DESEQ2_DEGX <- DESEQ2_table_dge_tb %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

setwd(raw_dir)

PA_info_genes <- data.table::fread("~/Desktop/PA_info_genes.txt")
merge_attempt <- merge(DESEQ2_table_dge_tb, PA_info_genes, by.x='gene', by.y='PA_num')

head(merge_attempt)

pdf("output_plots/volcano_plot_named.pdf")
EnhancedVolcano(merge_attempt, lab = merge_attempt$col_use, x='log2FoldChange', y='padj', title = 'volcano plot', pCutoff = 0.01, FCcutoff = 1.5)

dev.off()

setwd(wor_dir)

sig_genes_pth <- c(DESEQ2_DEGX$gene) #take names of significantly differentially expressed genes
sig_scores_pth <- c(DESEQ2_DEGX$padj) #take their associated adjusted pval
gene_ratio_input <- data.frame(gene_ids = c(sig_genes_pth), gene_scores = c(sig_scores_pth)) #merge into format required for downstream

#import GO annotations from pseudomonas.com gene ontology - need to find a way to get this to be uniform to all gene ontologies, but currently works well

PAO1_GO_all<-data.table::fread("~/Desktop/gene_ontology_csv.csv") #don't hardcode table - let user import

term2name <- subset(PAO1_GO_all, select=c(5,6,7)) #subsetting the relevant columns
geneList <- gene_ratio_input[,2] #list out genes
names(geneList) <- as.character(gene_ratio_input[,1]) #file formatting for downstream steps
geneList <- sort(geneList, decreasing=TRUE) #sorting geneList relative to pval
gene <- names(geneList)[abs(geneList)<0.05] #extracting names of genes with pval less than 0.05 (error in here i think - as i only used sig genes above for the relA_gene_ration_appmpt df)
term2gene <- subset(PAO1_GO_all, select=c(5,1)) #extracting relevant geneA = GO:X for next step)
term2name <- subset(PAO1_GO_all, select=c(5,6,7)) #extracting relevant GO# = function = category eg: GO:0005524 _ ATP binding _ molecular_function

#enrichment analysis on all GO

cluster_profiler_enriched <- enricher(gene, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene, TERM2NAME = term2name) #gene ontology enrichment analysis, bonferoni correction, all GOs)

if (is.na(cluster_profiler_enriched[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  }else {
  setwd(raw_dir)
  dotplot(cluster_profiler_enriched, showCategory=15) + ggtitle("GO gene ratio (all ontologies)")
  ggsave(filename = "output_plots/GO_gene_ratio_enriched.pdf")
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

if (is.na(cluster_profiler_enriched_MF[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  }else {
  setwd(raw_dir)
  dotplot(cluster_profiler_enriched_MF, showCategory=15) + ggtitle("GO gene ratio (Molecular function)")
  ggsave(filename = "output_plots/MF_gene_ratio_enriched.pdf")
  setwd(wor_dir)
  }


#BP only
PAO1_GO_BP <- PAO1_GO_all[ which(PAO1_GO_all$Namespace=='biological_process')] #same as above but with biological process - this will likely be the one you want to see.
term2name_BP <- subset(PAO1_GO_BP, select=c(5,6,7))
term2gene_BP <- subset(PAO1_GO_BP, select=c(5,1))
#term2gene_BP_aschar <- as.character(term2gene_BP$Accession)
#par_nodes_BP_PAO1 <- get_parent_nodes(term2gene_BP_aschar)
cluster_profiler_enriched_BP <- enricher(gene, qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE = term2gene_BP, TERM2NAME = term2name_BP)

if (is.na(cluster_profiler_enriched_BP[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  } else {
  setwd(raw_dir)
  dotplot(cluster_profiler_enriched_BP, showCategory=15) + ggtitle("GO gene ratio (Biological process)")
  ggsave(filename = "output_plots/BP_gene_ratio_enriched.pdf")
  setwd(wor_dir)
  }

#CC only

PAO1_GO_CC <- PAO1_GO_all[ which(PAO1_GO_all$Namespace=='cellular_component')] # and again with cellular components
term2name_CC <- subset(PAO1_GO_CC, select=c(5,6,7))
term2gene_CC <- subset(PAO1_GO_CC, select=c(5,1))
#term2gene_CC_aschar <- as.character(term2gene_CC$Accession)
#par_nodes_CC_PAO1 <- get_parent_nodes(term2gene_CC_aschar)
#pdf("output_plots/CC_gene_ratio_enriched.pdf")
#dotplot(cluster_profiler_enriched_CC, showCategory=15) + ggtitle("GO gene ratio (cellular component)") #will create generic cool looking dotplot
#dev.off()

cluster_profiler_enriched_CC <- enricher(gene, qvalueCutoff = 0.1, pAdjustMethod = "BH", TERM2GENE = term2gene_CC, TERM2NAME = term2name_CC)

if (is.na(cluster_profiler_enriched_CC[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  } else {
  setwd(raw_dir)
  dotplot(cluster_profiler_enriched_CC, showCategory=15) + ggtitle("GO gene ratio (all ontologies)")
  ggsave(filename = "output_plots/CC_gene_ratio_enriched.pdf")
  setwd(wor_dir)
  }

#cluster_profiler can also easily use kegg, which is a massive load off

clustprof_kegg <- enrichKEGG(gene = gene, organism = "pae", pvalueCutoff = 0.05)
#organism currently hardcoded which needs to change, this will give you pathways which are enriched from your significantly differentially expressed genes
#head(clustprof_kegg) #will show you what they are
#head(clustprof_kegg[,1]) #will show you specifically what pathways you are looking at significantly enriched

if (is.na(clustprof_kegg[1,5]) == "TRUE"){
  print("no significant enrichment terms")
  } else {
  setwd(raw_dir)
  dotplot(clustprof_kegg, showCategory=15) + ggtitle("GO gene ratio (KEGG pathways)")
  ggsave(filename = "output_plots/KEGG_gene_ratio_enriched.pdf")
  setwd(wor_dir)
  }

data_fold_changes <- DESEQ2_DEG$log2FoldChange
names(data_fold_changes) <- rownames(DESEQ2_DEG)

setwd(raw_dir)
setwd("pathview")

tmp <- pathview(gene.data = data_fold_changes, pathway.id = clustprof_kegg[,1], species = "pae", gene.idtype = "kegg", limit = list(gene= 2)) #this will output all of the differentially regulated pathways from pathview showing fancy pathway maps

setwd(wor_dir)

#pdf("output_plots/KEGG_gene_ratio_enriched.pdf")
#dotplot(clustprof_kegg, showCategory=15) + ggtitle("KEGG gene ratio (cellular component)") #will create generic cool looking dotplot
#dev.off()



quit(save="no")
