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

input_star <- DESeqDataSetFromHTSeqCount(sampleTable = input_table_DESEQ2, directory = (args[4]), design = Gene ~ Condition)

nrow(input_star)

input_star <- input_star[rowSums(counts(input_star)) > 10, ]

nrow (input_star)

input_star2 <- DESeq(input_star)

DESEQ2_norm_counts <- log2(counts(input_star2, normalized = TRUE)+1)

head(DESEQ2_norm_counts)

write.table(DESEQ2_norm_counts, "DESEQ2_norm_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")

resultsNames(input_star2)

DESEQ2_DEG <- results(object = input_star2, name="Condition_Un_vs_In")

head(DESEQ2_DEG)
