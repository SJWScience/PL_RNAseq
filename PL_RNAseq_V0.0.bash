#args: input_dir($1) file_suffix($2) output_dir($3) reference($4) threads($5) lib_type($6)
#echo $1
#echo ""
#echo $2
#echo ""
#echo $3
#echo ""
#echo $4

usage="
PL_RNAseq version 0.0 12-apr-21

Programs required for operation -
: FastQC
: Trimmomatic
: STAR
: R
: DESeq2 		(R package, installed via bioconductor not CRAN)
: apeglm 		(R package, installed via bioconductor not CRAN)
: pheatmap 	(R package)


Usage: ./`basename $0` -i <arg> -s <arg> -o <arg> -r <arg> -t <arg> -l <arg>

	-i --INPUT_DIR	Input directory of raw fastq files for this project
			Eg: ~/project1/raw_data/
			Make sure to include all forward slashes ( / )

	-s --SUFFIX	Suffix of gene files to be used in this analysis
			Eg: .fq.gz
			Make sure to include the file extention including period.

	-o --OUTPUT_DIR	Output directory that outputs from this analysis will be put
			Eg: ~/project1/output_RNA_pipeline/
			Make sure to include all /

	-r --REFERENCE	STAR reference genome directory to be used in this analysis
			Eg: /usr/local/bin/STAR_genomes/Pseudomonas_aeruginosa_PAO1_070421/
			Make sure not to include the actual files in the directory, just the path to the director (don't forget to inclue the /)

	-t --THREADS number of threads (virtual cores) you want this pipeline to be run on. Must be set.

	-l --LIB_TYPE Illumina library type (Single end (SE) or Paired end (PE))
			Eg: -l SE
			Make sure to use capital letters for PE or SE


full example: ./`basename $0` -i ~/project1/raw_data/  -s .fq.gz -o ~/project1/output_RNA_pipeline/ -r /usr/local/bin/STAR_genomes/Pseudomonas_aeruginosa_PAO1_070421/ -t 6 -l SE"

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

ls ${INPUT_DIR}*${SUFFIX}
lgth=`echo -n ${SUFFIX} | wc -c`
lgthx=$((lgth+1))
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)raw_data_fastQC
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed
mkdir ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed_data_fastQC
for i
in $(ls ${INPUT_DIR}*${SUFFIX} | xargs -n 1 basename | cut -f 1 -d '.' -)
do
echo "$i"

fastqc -t ${THREADS} -o ${OUTPUT_DIR}$(date +%Y%m%d_)raw_data_fastQC/ ${INPUT_DIR}"$i"${SUFFIX}

java -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar ${LIB_TYPE} -threads ${THREADS} -phred33 ${INPUT_DIR}"$i"${SUFFIX} ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed/"$i".trimmed.fastq.gz ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:20

fastqc -t ${THREADS} -o ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed_data_fastQC/ ${OUTPUT_DIR}/$(date +%Y%m%d_)trimmed/"$i".trimmed.fastq.gz

mkdir ${OUTPUT_DIR}$(date +%Y%m%d_)"$i"_STAR_alignment
star --runThreadN ${THREADS} \
--genomeDir ${REFERENCE} \
--readFilesCommand gunzip -c \
--readFilesIn  ${OUTPUT_DIR}$(date +%Y%m%d_)trimmed/"$i".trimmed.fastq.gz \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${OUTPUT_DIR}$(date +%Y%m%d_)"$i"_STAR_alignment/"$i" \
--outSAMattributes Standard

done
