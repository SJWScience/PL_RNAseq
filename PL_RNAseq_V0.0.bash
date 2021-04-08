#args: input_dir($1) file_suffix($2) output_dir($3) database($4)
#echo $1
#echo ""
#echo $2
#echo ""
#echo $3
#echo ""
#echo $4

usage="
PL_RNAseq version 0.0 07-apr-21

Usage: ./`basename $0` -i <arg> -s <arg> -o <arg> -d <arg>
	-i --INPUT_DIR	Input directory of raw fastq files for this project
			Eg: ~/project1/raw_data/
			Make sure to include all forward slashes ( / )

	-s --SUFFIX	Suffix of gene files to be used in this analysis
			Eg: .fasta
			Make sure to include the file extention including period.

	-o --OUTPUT_DIR	Output directory that outputs from this analysis will be put
			Eg: /Volumes/userdata/staff_groups/lamontlab/output/
			Make sure to include all /

	-d --DATABASE	BLAST database to be used in this analysis
			Eg /Volumes/userdata/staff_groups/lamontlab/output/BLASTdb
			Make sure not to include the suffix for the BLAST databases otherwise BLAST will error

full example: ./`basename $0` -i /Volumes/userdata/staff_groups/lamontlab/genes/ -s .fasta -o /Volumes/userdata/staff_groups/lamontlab/output/ -d /Volumes/userdata/staff_groups/lamontlab/output/BLASTdb"

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
    -d|--database)
    DATABASE="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done


export PATH=/Volumes/userdata/staff_groups/lamontlab/Documents/BLAST/output_files/mview-1.60.1/bin/:$PATH
command -v mview >/dev/null 2>&1 || { echo >&2 "This script requires mview but it's not installed, or not in your working environment.  Aborting."; exit 1; }
export PATH=/Volumes/userdata/student_users/samtaylorwardell/provean-1.1.5/bin:$PATH
command -v provean >/dev/null 2>&1 || { echo >&2 "This script requires provean but it's not installed, or not in your working environment.  Aborting."; exit 1; }

echo ""
if [ -z ${INPUT_DIR} ]
then
echo "$usage"
echo ""
echo "ERROR - No input directory set, please assign an input directory using the -i argument (example: -i ~/myfiles/inputgenes/"
exit 0
else
echo "INPUT_DIR = ${INPUT_DIR}"
fi
if [ -z ${SUFFIX} ]
then
echo "$usage"
echo ""
echo "ERROR - No file suffix set, please assign a file suffix using the -s argument (example: -s .fasta)"
exit 0
else
echo "SUFFIX = ${SUFFIX}"
fi
if [ -z ${OUTPUT_DIR} ]
then
echo "$usage"
echo ""
echo "ERROR - No output directory set, please assign an output directory using the -o argument (example: -o ~/myfiles/output/)"
exit 0
else
echo "OUTPUT_DIR = ${OUTPUT_DIR}"
fi
if [ -z ${DATABASE} ]
then
echo "$usage"
echo ""
echo "ERROR - No BLAST database set, please assign a BLAST database using the -d argument (example: -d ~/myfiles/blastdb/mydatabase)"
exit 0
else
echo "DATABASE = ${DATABASE}"
fi
echo ""
echo "input files being tested:"
echo ""
ls ${INPUT_DIR}*${SUFFIX}
lgth=`echo -n ${SUFFIX} | wc -c`
lgthx=$((lgth+1))
for i
in $(ls ${INPUT_DIR}*${SUFFIX} | xargs -n 1 basename | cut -f 1 -d '.' -)
do
echo "$i"
tmp_dir=$(mktemp -d -t varEV-XXXXXXXXXX --tmpdir=${INPUT_DIR})
cat ${INPUT_DIR}$i${SUFFIX} | tr -d " " > $tmp_dir/$i${SUFFIX}.tmp
tblastn -query $tmp_dir/$i${SUFFIX}.tmp -db ${DATABASE} -out ${OUTPUT_DIR}"$i".blast -num_threads 100 -max_hsps 1 -num_alignments 5000 -num_descriptions 5000 -seg no
echo "$i" MVIEW
mview -top 5000 -in blast ${OUTPUT_DIR}"$i".blast -out fasta > ${OUTPUT_DIR}"$i".blast.fasta
echo "$i" SNP-SITES
/Volumes/userdata/student_users/samtaylorwardell/bin/snp-sites/src/snp-sites  -v ${OUTPUT_DIR}"$i".blast.fasta -o ${OUTPUT_DIR}"$i".output_variants.vcf
awk '/#CHROM/ ,EOF { print $4,$2,$5 }' ${OUTPUT_DIR}"$i".output_variants.vcf | awk '{ n=split($3,arr,","); for(i=1;i<=n;i++) print $1,$2,arr[i] }' | tr -d "[:blank:]" | sed '1d' > ${OUTPUT_DIR}"$i".var
echo "$i" provean
if [ -f /Volumes/userdata/student_users/samtaylorwardell/provean-1.1.5/OUTPUTS/"$i".sss ]; then
provean.sh --num_threads 100 -q $tmp_dir/$i${SUFFIX}.tmp -v ${OUTPUT_DIR}"$i".var --supporting_set /Volumes/userdata/student_users/samtaylorwardell/provean-1.1.5/OUTPUTS/"$i".sss --verbose > ${OUTPUT_DIR}"$i".provean 2>&1
else
provean.sh --num_threads 100 -q $tmp_dir/$i${SUFFIX}.tmp -v ${OUTPUT_DIR}"$i".var --save_supporting_set /Volumes/userdata/student_users/samtaylorwardell/provean-1.1.5/OUTPUTS/"$i".sss --verbose > ${OUTPUT_DIR}"$i".provean 2>&1
fi
rm -rf $tmp_dir
done
