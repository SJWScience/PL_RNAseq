#load input data
args = commandArgs(trailingOnly=TRUE)
args

raw_dir <- (as.character(args[2]))
wor_dir <- (as.character(args[1]))
setwd(wor_dir)
input_table_DESEQ2 <- read.table(as.character(args[3]), header=T, sep="\t")

rownames(input_table_DESEQ2) <- input_table_DESEQ2$SampleName

head(input_table_DESEQ2)
nrow(input_table_DESEQ2)
ncol(input_table_DESEQ2)

library(DESeq2)
library(apeglm)
library(pheatmap)

#input_star <- DESeqDataSetFromHTSeqCount(sampleTable = input_table_DESEQ2, directory = (args[4]), design =  ~ Condition + Gene)
input_star <- DESeqDataSetFromHTSeqCount(sampleTable = input_table_DESEQ2, directory = (args[4]), design =  ~ Condition)

nrow(input_star)

input_star <- input_star[rowSums(counts(input_star)) > 10, ]

nrow (input_star)

input_star2 <- DESeq(input_star)

DESEQ2_norm_counts <- log2(counts(input_star2, normalized = TRUE)+1)

head(DESEQ2_norm_counts)

write.table(DESEQ2_norm_counts, "DESEQ2_norm_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")

resultsNames(input_star2)

DESEQ2_DEG <- results(object = input_star2, name="Condition_Un_vs_In")
##DESEQ2_DEG_alt <- results(object = input_star2, name="Gene_spoT_vs_relA")
##this part needs to not be hardcoded##

head(DESEQ2_DEG)
##head(DESEQ2_DEG_alt)

DESEQ2_DEG_shrink <- lfcShrink(dds = input_star2, coef="Condition_Un_vs_In", type="apeglm")
head(DESEQ2_DEG_shrink)

##DESEQ2_alt_DEG_shrink <- lfcShrink(dds = input_star2, coef="Gene_spoT_vs_relA", type="apeglm")
##head(DESEQ2_alt_DEG_shrink)

write.table(DESEQ2_DEG_shrink, "../DESEQ2_DEG_shrink.txt", quote=F, col.names=T, row.names=T, sep="\t")
##write.table(DESEQ2_alt_DEG_shrink, "../DESEQ2_alt_DEG_shrink.txt", quote=F, col.names=T, row.names=T, sep="\t")

write.table(DESEQ2_DEG, "../DESEQ2_DEG.txt", quote=F, col.names=T, row.names=T, sep="\t")
#write.table(DESEQ2_DEG_alt, "../DESEQ2_DEG_alt.txt", quote=F, col.names=T, row.names=T, sep="\t")

DESEQ2_var_stabl <- vst(input_star2)

DESEQ2_DistMatrix <- as.matrix(dist(t(assay(DESEQ2_var_stabl))))

png("sample_distance_heatmap_DESEQ.png")
pheatmap(DESEQ2_DistMatrix)

png("DESEQ_PCA.png")
plotPCA(object = DESEQ2_var_stabl, intgroup = "Condition")

png("DESEQ_gene_PCA.png")
plotPCA(object = DESEQ2_var_stabl, intgroup = "Gene")

dev.off()
