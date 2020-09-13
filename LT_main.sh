#!/usr/bin/env bash

usage() { 
 echo "
Usage: $0
 Required arguments:
 -i|--input            Input directory containing raw reads
 -o|--output           Output directory
 -s|--sketch           Mash sketch of GenomeTrakr and PubMLST Salmonella sequences
 -r|--reference        Reference Sequence for SNP tree
 
 Optional arguments:
 -n|--neighbours       Number of neighbours to include per query [Default: 50]
 -l|--list             List of query genome paths
 -t|--threads          Number of threads [Default: 1]
 --population          Construct the tree in the context of the entire Salmonella population structure
 -h|--help             Display help message

 "
}

# initialize variables
n_threads=1
population="false"
neighbours=50
src_dir=$(dirname $0)
start_time=$(date +%s)

# parse arguments
if [ $# == 0 ]
then
    usage
    exit 0
fi

opts=`getopt -o hi:o:t:s:l:n:r: -l help,input:,output:,sketch:,population,list:,neighbour:,reference: -- "$@"`
if [ $? != 0 ] ; then echo "WARNING: Invalid arguements used, exiting"; usage; exit 1 ; fi
eval set -- "$opts"

while true; do
  case "$1" in
    -i|--input) IN_DIR=$2; shift 2 ;;
    -l|--list) list=$2; shift 2 ;;
    -s|--sketch) ref_sketch=$2; shift 2 ;;
    -r|--reference) reference=$2; shift 2 ;;
    -n|--neighbours) neighbours=$2; shift 2 ;;
    --population) population="true"; shift ;;
    -o|--output) OUT_DIR=$2; shift 2 ;;
    -t) n_threads=$2; shift 2 ;;
    --) shift; break ;;
    -h|--help) usage; exit ;;
  esac
done

getreadnames() {
    input_array=($@)
    for i in ${!input_array[*]}; do
        if test -f ${input_array[${i}]}/*R1*; then
            seq_1=$(echo ${input_array[${i}]}/*R1*)
        else
            echo "Forward reads cannot be found for sample ${i}, exiting"; exit 1
        fi
        
        if test -f ${input_array[$i]}/*R2*; then
            seq_2=$(echo ${input_array[${i}]}/*R2*)
        else
            echo "Reverse reads cannot be found for sample ${i}, exiting"; exit 1
        fi
        
        input_array[${i}]=$(echo ${seq_1},${seq_2})
    done
}

getcontignames() {
    contig_array=($@)
    for i in ${!contig_array[*]}; do
        contig_array[${i}]=$(echo ${contig_array[${i}]},${contig_array[${i}]}"/"$(basename ${contig_array[${i}]}".fa"))
    done
}

# Genome assembly
assembly() {  
  
  # $1 Forward Reads
  # $2 Reverse Reads

  # Get filename  
  filename=$(basename $(dirname $1))
  
  # Read trimming using fastp
  time=$(date +"%T")
  echo "[$time] Read trimming on "$(basename $1)" and "$(basename $2)

  cleaned_fastq_1=$OUT_DIR/LT_Senterica_tmp/fastp_res/$(basename ${1%.*})_clean.fastq
  cleaned_fastq_2=$OUT_DIR/LT_Senterica_tmp/fastp_res/$(basename ${2%.*})_clean.fastq
  fastp -i $1 -I $2 \
        -o $cleaned_fastq_1 \
        -O $cleaned_fastq_2 \
        -w $n_threads
  
  # Genome assembly using shovill
  time=$(date +"%T")
  echo "[$time] Genome assembly on "$(basename $1)" and "$(basename $2)

  mkdir -p $OUT_DIR/LT_Senterica_tmp/shovill_res/$filename
  shovill --outdir $OUT_DIR/LT_Senterica_tmp/shovill_res/$filename \
          --R1 $cleaned_fastq_1 \
          --R2 $cleaned_fastq_2 \
          --gsize 4.5M \
          --cpus $n_threads \
          --force
  
  # rename shovill output
  mv $OUT_DIR/LT_Senterica_tmp/shovill_res/$filename/contigs.fa $OUT_DIR/LT_Senterica_tmp/shovill_res/$filename/$filename.fa
}

# Reference sequence comparisons using mash
mash_exe() {

    mash sketch -p $n_threads -s 10000 -o $OUT_DIR/LT_Senterica_tmp/mash_res/sketch/$1.msh $2
    mash dist $OUT_DIR/LT_Senterica_tmp/mash_res/sketch/$1.msh $ref_sketch -t -p $n_threads > $OUT_DIR/LT_Senterica_tmp/mash_res/results/$1.tab
}

# identify top hits
process_ref_seq_mash() {

    # $1 Query assembly path

    # identify top hits
    $src_dir/src/top_hit.R $1 $2 $3 > $OUT_DIR/LT_Senterica_tmp/mash_res/process/$(basename ${1%.*})_neighbours.csv
}

### MAIN

## check if required arguments are present
if [[ -z $OUT_DIR ]]; then usage; echo "Required argument -o is missing, exiting"; exit 1; fi
if [[ -z $IN_DIR ]]; then usage; echo "Required argument -i is missing, exiting"; exit 1; fi
if [[ -z $ref_sketch ]]; then usage; echo "Required argument -s is missing, exiting"; exit 1; fi
if [[ -z $reference ]]; then usage; echo "Required argument -r is missing, exiting"; exit 1; fi
if ! test -f $reference; then echo "Specified reference sequence cannot be found, exiting"; exit 1; fi

## create LT_Senterica_tmp directory structure
mkdir -p $OUT_DIR/LT_Senterica_tmp/fastp_res
mkdir -p $OUT_DIR/LT_Senterica_tmp/shovill_res
mkdir -p $OUT_DIR/LT_Senterica_tmp/mash_res/sketch
mkdir -p $OUT_DIR/LT_Senterica_tmp/mash_res/results
mkdir -p $OUT_DIR/LT_Senterica_tmp/mash_res/process
mkdir -p $OUT_DIR/LT_Senterica_tmp/genomes
mkdir -p $OUT_DIR/LT_Senterica_tmp/phame_input
mkdir -p $OUT_DIR/LT_Senterica_tmp/phame_ref

## declare input file array
if test -d $IN_DIR; then
    declare -a input_array=($IN_DIR/*)
    samples_n=${#input_array[@]}
    echo "Total number of samples found: ${samples_n}"
    getreadnames ${input_array[@]}
else
    echo "Specified input raw reads directory does not exist, exiting"; exit 1
fi

## genome assembly
for i in ${!input_array[*]}; do
    seq_1=$(echo ${input_array[${i}]} | cut -d, -f1)
    seq_2=$(echo ${input_array[${i}]} | cut -d, -f2)
    assembly $seq_1 $seq_2
done

## declare contig file array
declare -a contig_array=($OUT_DIR/LT_Senterica_tmp/shovill_res/*)
getcontignames ${contig_array[@]}

## reference sequence comparisons
time=$(date +"%T")
echo "[${time}] Query to reference sequence comparisons using mash..."

# shovill-assembled sequences
for i in ${!contig_array[*]}; do
    filename=$(basename $(echo ${contig_array[${i}]} | cut -d, -f1))
    filepath=$(echo ${contig_array[${i}]} | cut -d, -f2)
    mash_exe $filename $filepath
done

# pre-assembled sequences
if [[ -z $list ]]; then
    time=$(date +"%T")
    echo "[${time}] List of assembled sequences was not provided, skipping..."
elif test -f $list; then
    while read lines; do
        if test -f $lines; then
            filename=$(basename $lines)
            mash_exe ${filename%.*} $lines
        else
            echo "The query genome cannot be found at ${lines}, skipping..."
        fi
    done < $list
else
    echo "[${time}] The given list of genomes does not exist, skipping..."
fi


## process reference sequence comparisons mash results
time=$(date +"%T")
echo "[${time}] Identifying top reference sequence hits..."

declare -a ref_seq_mash_array=($OUT_DIR/LT_Senterica_tmp/mash_res/results/*)
for i in ${!ref_seq_mash_array[*]}; do
    process_ref_seq_mash ${ref_seq_mash_array[${i}]} $src_dir/metadata/Senterica_population_metadata_V2.tsv $neighbours
done

## download candidate genomes

# identify all unique candidate genomes for download
if [[ $population == "false" ]]; then
    cat $OUT_DIR/LT_Senterica_tmp/mash_res/process/* | sort | uniq > $OUT_DIR/LT_Senterica_tmp/candidate_genome_list.csv
else
    cat $OUT_DIR/LT_Senterica_tmp/mash_res/process/* $src_dir/metadata/ref_population_ftp.csv | sort | uniq > $OUT_DIR/LT_Senterica_tmp/candidate_genome_list.csv
fi

genomes_n=$(wc -l $OUT_DIR/LT_Senterica_tmp/candidate_genome_list.csv | cut -f1 -d' ')
time=$(date +"%T")
echo "[$time] Downloading a total of $genomes_n genomes from NCBI and BIGSdb..."

# call GNU Parallel
cat $OUT_DIR/LT_Senterica_tmp/candidate_genome_list.csv | parallel -j 10 $src_dir/src/wget.sh {} $OUT_DIR

## Construct SNP tree with Phame
# set up phame input
for i in $OUT_DIR/LT_Senterica_tmp/genomes/*; do # downloaded genomes
    ln -s $(realpath $i) $OUT_DIR/LT_Senterica_tmp/phame_input/$(basename $i)
done

for i in $OUT_DIR/LT_Senterica_tmp/shovill_res/*; do # shovill assembled genomes
    ln -s $(realpath $i)/$(basename $i).fa $OUT_DIR/LT_Senterica_tmp/phame_input/$(basename $i).fa
done

if [[ ! -z $list ]]; then # pre-assembled genomes
    while read lines; do
        if test -f $lines; then
            ln -s $lines $OUT_DIR/LT_Senterica_tmp/phame_input/$(basename $lines)
        fi
    done < $list
fi

ln -s $(realpath $reference) $OUT_DIR/LT_Senterica_tmp/phame_ref/Reference.fa # reference genome

# call phame
time=$(date +"%T")
echo "[${time}] Constructing cgSNP tree"

$src_dir/src/generatePhameCtl.sh $OUT_DIR/LT_Senterica_tmp/phame_input $OUT_DIR/LT_Senterica_tmp/phame_ref $OUT_DIR/LT_Senterica_tmp/phame_ref/Reference.fa $n_threads $OUT_DIR/LT_Senterica_tmp
phame $OUT_DIR/LT_Senterica_tmp/phame.ctl
cp $OUT_DIR/LT_Senterica_tmp/phame_input/results/trees/*.fasttree $OUT_DIR/tree.nwk

time=$(date +"%T")
echo "[${time}] Phylogenetic tree file written to: $OUT_DIR/tree.nwk"

## write microreact metadata file for tree annotations
for i in $OUT_DIR/LT_Senterica_tmp/shovill_res/*; do echo $(basename $i) >> $OUT_DIR/LT_Senterica_tmp/query.list; done
if [[ ! -z $list ]]; then while read lines; do echo $(basename $lines) >> $OUT_DIR/LT_Senterica_tmp/query.list; done < $list; fi
for i in $OUT_DIR/LT_Senterica_tmp/genomes/*; do echo $(basename $i .fa) >> $OUT_DIR/LT_Senterica_tmp/neighbours.list; done
echo "Reference" >> $OUT_DIR/LT_Senterica_tmp/reference.list

if [[ $population == "true" ]]; then
    $src_dir/src/microreact.R $OUT_DIR/LT_Senterica_tmp/query.list \
        $OUT_DIR/LT_Senterica_tmp/neighbours.list \
        $OUT_DIR/LT_Senterica_tmp/reference.list \
        $src_dir/metadata/Senterica_population_metadata_V2.tsv \
        $src_dir/metadata/ref_population_ftp.tsv > $OUT_DIR/microreact_metadata.csv
else
    $src_dir/src/microreact.R $OUT_DIR/LT_Senterica_tmp/query.list \
        $OUT_DIR/LT_Senterica_tmp/neighbours.list \
        $OUT_DIR/LT_Senterica_tmp/reference.list \
        $src_dir/metadata/Senterica_population_metadata_V2.tsv \
        false > $OUT_DIR/microreact_metadata.tsv
fi

time=$(date +"%T")
echo "[${time}] Microreact tree annotation metadata file written to: $OUT_DIR/microreact_metadata.tsv"

## Say Goodbye
end_time=$(date +%s)
time=$(date +"%T")
echo "[$time] Pipeline finished!"
echo "[${time}] Run time: $((($end_time - $start_time)/60)) minutes"