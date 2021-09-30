#! /bin/bash
#args: input_dir($1) file_suffix($2) output_dir($3) reference($4) threads($5) lib_type($6)
#echo $1
#echo ""
#echo $2
#echo ""
#echo $3
#echo ""
#echo $4

usage="

HELP

PL_RNAseq version 0.0 12-apr-21

Programs required for operation -
: FastQC
: Trimmomatic
: STAR
: R
: DESeq2 	(R package, installed via bioconductor not CRAN)
: apeglm 	(R package, installed via bioconductor not CRAN)
: pheatmap 	(R package)


Usage: `basename $0` -i <arg> -s <arg> -o <arg> -r <arg> -t <arg> -l <arg>

	-i --INPUT_DIR	Input directory of raw fastq files for this project
			Eg: ~/project1/raw_data


	-s --SUFFIX	Suffix of gene files to be used in this analysis
			Eg: .fq.gz
			Make sure to include the file extention including period.

	-o --OUTPUT_DIR	Output directory that outputs from this analysis will be put
			Eg: ~/project1/output_RNA_pipeline

	-r --REFERENCE	STAR reference genome directory to be used in this analysis
			Eg: /usr/local/bin/STAR_genomes/Pseudomonas_aeruginosa_PAO1_070421
			Make sure not to include the actual files in the directory, just the path to the directory

	-t --THREADS number of threads (virtual cores) you want this pipeline to be run on. Must be set.

	-l --LIB_TYPE Illumina library type (Single end (SE) or Paired end (PE))
			Eg: -l SE
			Make sure to use capital letters for PE or SE


full example: `basename $0` -i ~/project1/raw_data  -s .fq.gz -o ~/project1/output_RNA_pipeline -r /usr/local/bin/STAR_genomes/Pseudomonas_aeruginosa_PAO1_070421 -t 6 -l SE

Aditional options

Usage

`basename $0` -refcreate

Will call a script to create a reference genome for STAR to align reads to

`basename $0` -spliced

Will call a script to perform RNAseq mapping on eukaryote samples

`basename $0` -mixed

Will allow for RNA samples which contain material from multiple species of prokaryotes

`basename $0` -whywouldyoudothistome

Will eventually allow for RNA sameples which contain material from multiple species of prokaryotes or eukaryotes "

if [ "$1" == "-h" ]; then
  echo "$usage"
  exit 0
fi
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -i|--input_dir)
    INPUT_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--suffix)
    SUFFIX="$2"
    shift # past argument
    shift # past value
    ;;
    -o|--output_dir)
    OUTPUT_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    -r|--reference)
    REFERENCE="$2"
    shift # past argument
    shift # past value
    ;;
		-t|--threads)
		THREADS="$2"
		shift # past argument
		shift # past value
		;;
		-l|--lib_type)
		LIB_TYPE="$2"
		shift # past argument
		shift # past value
		;;
		-refcreate|--reference_create)
		REFERENCE_CREATE="$2"
		shift # past argument
		shift # past value
		;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done


command -v star >/dev/null 2>&1 || { echo >&2 "This script requires STAR but it's not installed, or not in your working environment.  Aborting."; exit 1; }
export PATH=/usr/local/bin/FastQC:$PATH
command -v fastqc >/dev/null 2>&1 || { echo >&2 "This script requires fastqc but it's not installed, or not in your working environment.  Aborting."; exit 1; }
export PATH=/usr/local/bin/Trimmomatic-0.39:$PATH
command -v java -jar trimmomatic-0.39.jar >/dev/null 2>&1 || { echo >&2 "This script requires Trimmomatic but it's not installed, or not in your working environment.  Aborting."; exit 1; }

#if [ -n ${REFERENCE_CREATE}]
#then
#  echo "1"
#  exit 0
#else
#  echo "2"
#fi

## Code for parsing arguments, allows for them to be placed in any order within the run command ##
echo ""
if [ -z ${INPUT_DIR} ]
then
echo "$usage"
echo ""
echo "ERROR - No input directory set, please assign an input directory using the -i argument (example: -i ~/project1/raw_data/"
exit 0
else
echo "INPUT_DIR = ${INPUT_DIR}"
fi
if [ -z ${SUFFIX} ]
then
echo "$usage"
echo ""
echo "ERROR - No file suffix set, please assign a file suffix using the -s argument (example: -s .fq.gz)"
exit 0
else
echo "SUFFIX = ${SUFFIX}"
fi
if [ -z ${OUTPUT_DIR} ]
then
echo "$usage"
echo ""
echo "ERROR - No output directory set, please assign an output directory using the -o argument (example: -o ~/project1/output_RNA_pipeline/)"
exit 0
else
echo "OUTPUT_DIR = ${OUTPUT_DIR}"
fi
if [ -z ${REFERENCE} ]
then
echo "$usage"
echo ""
echo "ERROR - No STAR reference genome set, please assign an STAR reference genome directory using the -r argument (example: -r /usr/local/bin/STAR_genomes/Pseudomonas_aeruginosa_PAO1_070421/)"
exit 0
else
echo "REFERENCE = ${REFERENCE}"
fi
if [ -z ${THREADS} ]
then
echo "$usage"
echo ""
echo "ERROR - No threads designated, please use the -t argument and add the number of cores you want this pipeline run on (example: -t 6)"
exit 0
else
echo "THREADS = ${THREADS}"
fi
if [ -z ${LIB_TYPE} ]
then
echo "$usage"
echo ""
echo "ERROR - No library type assigned, please use the -l argument and add the library type of your experience (paired end or single end) (example: -l SE)"
exit 0
else
echo "LIB_TYPE = ${LIB_TYPE}"
fi
echo ""
echo "input files being tested:"
echo ""

## Making folders for the outputs - using date as a prefix to prevent accidental overwriting of runs done previously ##
ls ${INPUT_DIR}*${SUFFIX}
lgth=`echo -n ${SUFFIX} | wc -c`
lgthx=$((lgth+1))
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)raw_data_fastQC
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed_data_fastQC
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)deseq2_inputs

## start of for loop to perform mapping and read counting on each sample within your folder matching the suffix provided ##
for i
in $(ls ${INPUT_DIR}*${SUFFIX} | xargs -n 1 basename | cut -f 1 -d '.' -) ## Takes all files within your input directory matching your suffix provided. Cuts off the suffix (eg .fq.gz) and then sets the file name as the basename (assigned to variable $i. For example a sample called "sample1.fq.gz" becomes "sample1" and when $i variable is used it recognises that as "sample1" ##
do
echo "$i"
## Added so you know exactly which files are being loaded into the pipeline, enabling you to stop the pipeline if you accidentally included files you didn't want ##

fastqc -t ${THREADS} -o ${OUTPUT_DIR}$(date +%Y%m%d_)raw_data_fastQC/ ${INPUT_DIR}"$i"${SUFFIX}
## Prelim fastQC on raw reads ##

java -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar ${LIB_TYPE} -threads ${THREADS} -phred33 ${INPUT_DIR}/"$i"${SUFFIX} ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed/"$i".trimmed.fastq.gz ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:20
## Trimmomatic on raw_reads, in this case it assumes your primers are found in the TruSeq3 SE primer list, this may need to be changed if you have data from different sources, be aware of it ##

fastqc -t ${THREADS} -o ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed_data_fastQC/ ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed/"$i".trimmed.fastq.gz ## Post trim, fastQC ##

## STAR alignment, this can be customised to suit your needs ##
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)"$i"_STAR_alignment
star --runThreadN ${THREADS} \
--genomeDir ${REFERENCE} \
--readFilesCommand gunzip -c \
--readFilesIn  ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed/"$i".trimmed.fastq.gz \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${OUTPUT_DIR}/$(date +%Y%m%d_)"$i"_STAR_alignment/"$i" \

## Extracting the relevant columns from the STAR output, ready for parsing into R and DEseq2 ##
cut -f1,4 ${OUTPUT_DIR}/$(date +%Y%m%d_)"$i"_STAR_alignment/"$i"ReadsPerGene.out.tab | grep -v "_" > ${OUTPUT_DIR}/$(date +%Y%m%d_)deseq2_inputs/`basename ${OUTPUT_DIR}/$(date +%Y%m%d_)"$i"_STAR_alignment/"$i" ReadsPerGene.out.tab`_counts.txt

## Creating an info sheet to be parsed into the R script ##
cat <(echo -e "SampleName\tFileName\tGene\tCondition") <(paste <(ls ${OUTPUT_DIR}/$(date +%Y%m%d_)deseq2_inputs | cut -d"_" -f1-4) <(ls ${OUTPUT_DIR}/$(date +%Y%m%d_)deseq2_inputs) <(ls ${OUTPUT_DIR}/$(date +%Y%m%d_)deseq2_inputs | cut -d"_" -f1 | awk '{print $0}') <(ls ${OUTPUT_DIR}/$(date +%Y%m%d_)deseq2_inputs | cut -d"_" -f2 | awk '{print $0}')) > ${OUTPUT_DIR}/$(date +%Y%m%d_)sample_sheet_DEseq.txt
INPUT_DESEQ=${OUTPUT_DIR}/$(date +%Y%m%d_)sample_sheet_DEseq.txt
done
R_WD=$(pwd ${INPUT_DIR})
IN_R=${OUTPUT_DIR}/$(date +%Y%m%d_)deseq2_inputs/
Rscript PL_RNAseq_V0.0.R ${R_WD} ${1} ${OUTPUT_DIR} ${2} $INPUT_DESEQ ${3} $IN_R ${4} #initialise R script and pass pass arguments
