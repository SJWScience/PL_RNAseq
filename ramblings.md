for QC checks

fastqc outputs a zip file that contains a 'fastqc_data.txt' and 'summary.txt'

fastqc_data.txt contains a line with total sequences and %GC - these can be extracted and averaged across all samples in the experiment and used to identify if there is an outlier (not same species or too few reads)

summary.txt contains pass, fail, warning for each of the standard fastqc metrics

add in an option to be verbose about the QC data and output a QC log file for the user to check. - from the summary.txt of each fastqc output zip files. Add in pre- and post- trimming QC??? not sure yet.


# 7th April 2021

STAR is quite finicky and seems overkill for bacterial RNAseq. However, things to be mindful of in the future.

  - Making reference genomes in STAR only uses the 3rd column for what they call "exon" identification. This is a problem with Bacterial GTF files as typically they will only contain gene and CDS fields. So you need to use additional fields to tell it not to use exon, rather to use CDS. This also extends to "transcript ID/transcript parent" that mammalian genomes have. Those need to be changes relative to the GTF file you have. See code below for how the reference genome was generated

```bash

star --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /usr/local/bin/STARgenomes/Pseudomonas_aeruginosa_PAO1_070421 \
--genomeFastaFiles /usr/local/bin/STARgenomes/Pseudomonas_aeruginosa_PAO1_070421/GCF_000006765.1_ASM676v1_genomic.fna \
--sjdbGTFfeatureExon gene \
--sjdbGTFfile /usr/local/bin/STARgenomes/Pseudomonas_aeruginosa_PAO1_070421/GCF_000006765.1_ASM676v1_genomic.gff \
--genomeSAindexNbases 10 \
--sjdbGTFtagExonParentTranscript Name \
--sjdbGTFtagExonParentGene locus_tag

```

  - As you can see above the sjdbGTFfeatureExon had to be changes to "gene" and the sjdbGTFtagExonParentTranscript and sjdbGTFtagExonParentGene also needed to be changed. This give an output reference file that has (in Pseudomonas_aeruginosa_PAO1) the PA numbers as gene annotations. I should/could refine to make it show common gene names as well as PA numbers but i am not sure the use of this right now in the early stages.

  - star itself is quite easy to install however. I am going to test its results against something a little quicker and simpler to use like Kallisto or Salmon. Just because the downstream DGE steps are a little easier using those programs as the outputs can be parsed into DEseq2 without much editing.

  - Must be noted. Cannot use the GTF files from pseudomonas.com they aren't in a format that star can handle, the GTF files from them has additional fields before the info, and beyond recoding star or heavily modifying the pseudomonas.com GTF file its worth just collecting the files from NCBI and going from there. Its not perfect but it is a standardised format that works with these programs without much headache. Consider emailing Geoff Windsor and seeing why they are in that format, he might not know that there is an issue.

Running star is a little bit strange too. Super functional and adjustable but. . .well.

```bash

star --runThreadN 6 \
--genomeDir /usr/local/bin/STARgenomes/Pseudomonas_aeruginosa_PAO1_070421 \
--readFilesCommand gunzip -c \
--readFilesIn ~/Desktop/local_storage/PL_RNAseq_testing/trimmed/trimmed_spoT_Un_1.fastq.gz \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes Standard

```


  - As above shows, little things to keep in mind. If the fq files are compressed (when are they not!?) you need to add the gunzip command in there. Additionally for downstream processes like DEseq2 you need to output GeneCounts those files require a little bit of extra work after though (see example below).

  - One other thing i had to do for this to work (it was working without the --outSAMtype and --outSAMattributes) was to run the command 'ulimit -n 10000'. before i did that the command was throwing up an error of "BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file ./_STARtmp//BAMsort/4/49 SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files."


```R
    N_unmapped	977798	977798	977798
    N_multimapping	1130303	1130303	1130303
    N_noFeature	102649	3135934	362250
    N_ambiguous	273035	9444	206518
    PA0001	2502	24	2478
    PA0002	1846	45	1801
    PA0003	1227	10	1217
    PA0004	3803	89	3714
    PA0005	371	7	365
    PA0006	166	8	158
    PA0007	190	29	163
    PA0008	1346	7	1346
    PA0009	890	6	886
    PA0010	53	3	51
    PA0011	330	17	313
    PA0012	78	2	76
    PA0013	308	9	299
    PA0014	40	20	20
    PA0015	109	11	98
    PA0016	1026	21	1010
    PA0017	531	4	528
    PA0018	377	2	377
    PA0019	1166	15	1156
```


  This is an example output from one sample (from code above) according to documentation thise columns represent
- column 1: gene ID
- column 2: counts for unstranded RNA-seq.
- column 3: counts for the 1st read strand aligned with RNA
- column 4: counts for the 2nd read strand aligned with RNA (the most common protocol nowadays). Because this was single end RNAseq its got low numbers in the '1st read strand'. For SE reads i will take Col4.

  - notes little finicky things like if you are putting in more than one sample as part of an experiment you typically sperate them by a comma with no space (sample1,sample2). The problem with this is you need absolute file paths for this, if you use shortcuts like $HOME or ~ it will error out because it doesnt parse it properly. so make sure you use absolute paths in files num_alignments


# 8th April 2021

Adding in additional files to the analysis

```bash
star --runThreadN 6\
--genomeDir /usr/local/bin/STARgenomes/Pseudomonas_aeruginosa_PAO1_070421 \
--readFilesCommand gunzip -c \
--readFilesIn ~/Desktop/local_storage/PL_RNAseq_testing/trimmed/trimmed_spoT_Un_1.fastq.gz,/Users/sam_2021/Desktop/local_storage/PL_RNAseq_testing/trimmed/trimmed_spoT_Un_2.fastq.gz,/Users/sam_2021/Desktop/local_storage/PL_RNAseq_testing/trimmed/trimmed_spoT_Un_3.fastq.gz \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes Standard \
--limitBAMsortRAM 1238776520

```

  - limitBAMsortRAM had to be added because without it the command would error out due to not enough RAM (frustrating).

  - The above code allows you to essentially concat files together (like you would do for technical replicates essentially). Not much use for independent replicates (i think). So run them seperate and just loop it to get it tearing through them.

  ```bash

for i
in *fastq.gz;
do
STAR --runMode alignReads --readFilesCommand gunzip -c --outSAMtype BAM Unsorted --genomeDir /usr/local/bin/STARgenomes/Pseudomonas_aeruginosa_PAO1_070421 --readFilesIn $i --runThreadN 6 --outFileNamePrefix /Users/sam_2021/Desktop/local_storage/PL_RNAseq_testing/STAR_alignments/first_iteration_defaults_basic/${i%fastq.gz} --limitBAMsortRAM 1238776520 --outSAMattributes Standard --quantMode GeneCounts
done

```

  - This loop works. Next step formatting outputs for input into DEseq2


# 9th April 2021

  - first steps are to extract relevant info from the STAR alignment gene GeneCounts, can be done with a simple loop.

```bash

for i in *ReadsPerGene.out.tab                   
do echo $i
cut -f1,4 $i | grep -v "_" > ../../deseq2_inputs/`basename $i ReadsPerGene.out.tab`_counts.txt
done

```

  - Next is to make the "experimental design" document that will load into R for DESeq2, you can extract all this info from your samples if you wanted to and have named them accordingly (this is what i did). Or you can just make your own table and save it and import it into R.

```bash

cat <(echo -e "SampleName\tFileName\tGene\tInduced") <(paste <(ls ../deseq2_inputs | cut -d"_" -f2-4) <(ls ../deseq2_inputs) <(ls ../deseq2_inputs | cut -d"_" -f2 | awk '{print $0}') <(ls ../deseq2_inputs | cut -d"_" -f3 | awk '{print $0}')) > ../sample_sheet_spoT.txt

```

  - that results in a file looking like this

```R
  SampleName	FileName	Gene	Induced
  spoT_In_1.	trimmed_spoT_In_1._counts.txt	spoT	In
  spoT_In_2.	trimmed_spoT_In_2._counts.txt	spoT	In
  spoT_In_3.	trimmed_spoT_In_3._counts.txt	spoT	In
  spoT_Un_1.	trimmed_spoT_Un_1._counts.txt	spoT	Un
  spoT_Un_2.	trimmed_spoT_Un_2._counts.txt	spoT	Un
  spoT_Un_3.	trimmed_spoT_Un_3._counts.txt	spoT	Un

```

  - not perfect, but it works for this purpose (a little messy with the period in the sample names etc. can refine)

  - Next is parsing it into R. I prefer to use Rstudio for this, however to make life easier you can parse it into R straight away through terminal.


```R


library(DESeq2)

spoT_sampletable <- read.table("../sample_sheet_spoT.txt", header=T, sep="\t")

rownames(spoT_sampletable) <- spoT_sampletable$SampleName

head(spoT_sampletable)
nrow(spoT_sampletable)
ncol(spoT_sampletable)

spoT_star <- DESeqDataSetFromHTSeqCount(sampleTable = spoT_sampletable, directory = "../deseq2_inputs/", design = ~ Induced)

nrow(spoT_star)

spoT_star <- spoT_star[rowSums(counts(spoT_star)) > 10, ]

nrow (spoT_star)

spoT_star2 <- DESeq(spoT_star)

spoT_norm_counts <- log2(counts(spoT_star2, normalized = TRUE)+1)

head(spoT_norm_counts)

write.table(spoT_norm_counts, "../spoT_normalized_counts.txt", quote=F, col.names=T, row.names=F, sep="\t")

resultsNames(spoT_star2)

spoT_DEG <- results(object = spoT_star2, name="Induced_Un_vs_In")

head(spoT_DEG)

library(apeglm)

spoT_DEG_shrink <- lfcShrink(dds = spoT_star2, coef="Induced_Un_vs_In", type="apeglm")
head(spoT_DGE)

head(spoT_DEG_shrink)

write.table(spoT_DEG_shrink, "../spoT_deseq2_results_shrink.txt", quote=F, col.names=T, row.names=T, sep="\t")

write.table(spoT_DEG, "../spoT_deseq2_results.txt", quote=F, col.names=T, row.names=T, sep="\t")

library(pheatmap)

spoT_var_stabl <- vst(spoT_star2)

spoT_DistMatrix <- as.matrix(dist(t(assay(spoT_var_stabl))))

png("sample_distance_heatmap_star.png")
pheatmap(spoT_DistMatrix)

png("../spoT_PCA.png")
plotPCA(object = spoT_var_stabl, intgroup = "Induced")

dev.off()

```

![PCA Plot](https://github.com/SJWScience/PL_RNAseq/blob/main/plots/spoT_PCA.png)
![Sample Heatmap](https://github.com/SJWScience/PL_RNAseq/blob/main/plots/sample_distance_heatmap_star.png)


  - This will give you a PCA plot, sample heat map, and differentially expressed genes. Its not pretty though and needs further refining. BUT it currently works. (see memo on In sample 2 being poor amount of reads - i could see this in the original fq.gz files, however i persisted to see where it would go, but it clearly is an outlier)
