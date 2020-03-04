#!/usr/bin/env bash

usage() { 
 echo "
Usage: $0
 Required:
 -i|--input     specify input directory
 -o|--output    specify output directory
 -s|--sketch    specify reference sequence mash sketch path
 -tmp           specify directory for storing temporary files
 
 Options:
 -t             Number of threads [Default: 1]
 -h|--help      Display help message

 "
}

# initialize variables
n_threads=1
tmp_dir=""

# parse arguments
if [ $# == 0 ]
then
    usage
    exit
fi

opts=`getopt -o hi:o:t:s: -l tmp,help,input,output,sketch -- "$@"`
if [ $? != 0 ] ; then echo "WARNING: Invalid arguements used, exiting"; usage; exit 1 ; fi
eval set -- "$opts"

while true; do
  case "$1" in
    -i|--input) IN_DIR=$2; shift 2 ;;
    -s|--sketch) ref_seq=$2; shift 2 ;;
    -o|--output) OUT_DIR=$2; shift 2 ;;
    -t) n_threads=$2; shift 2 ;;
    --tmp) tmp_dir=$2; shift 2 ;;
    --) shift; break ;;
    -h|--help) usage; exit ;;
  esac
done


getfilenames() {
    for i in "${!$1[@]}"; do
        seq_1=$(echo ${$1[${i}]}/*R1*)
        seq_2=$(echo ${$1[${i}]}/*R2*)
        $1[${i}]=$(echo ${seq_1},${seq_2})
    done
}

getcontignames() {
    for i in "${!$1[@]}"; do
        $1[${i}]=$(echo $1[${i}]"/contigs.fa")
    done
}

# Genome assembly
assembly() {
  time=$(date +"%T")
  echo "[$time] Genome assembly on "$(basename $1)" and "$(basename $2)
  filename=$(basename ${1%_R1*})
  # Read trimming using fastp  
  cleaned_fastq_1=$tmp_dir/fastp_res/$(basename ${$1%.*})_clean.fastq
  cleaned_fastq_2=$tmp_dir/fastp_res/$(basename ${$2%.*})_clean.fastq
  fastp -i $1 -I $2 \
        -o $cleaned_fastq_1 \
        -O $cleaned_fastq_2 \
        -w $n_threads &> ${tmp_dir}/log_files/fastp.log
  
  # Genome assembly using shovill
  mkdir $tmp_dir/shovill_res/$filename
  shovill --outdir $tmp_dir/shovill_res/$filename \
          --R1 $cleaned_fastq_1 \
          --R2 $cleaned_fastq_2 \
          --gsize 4.5M \
          --cpus $n_threads \
          --force &> $tmp_dir/log_files/shovill.log
}

# Reference sequence comparisons using mash
 mash_exe() {
    time=$(date +"%T")
    echo "[$time] Reference sequence comparisons using mash"

    mash sketch -p $n_threads -s 10000 -o $tmp_dir/$filename
    mash dist $tmp_dir/$filename $ref_seq -p $n_threads > $tmp_dir/mash_res.tab
}

# MAIN

# create tmp directory structure
mkdir -p $tmp_dir/log_files
mkdir -p $tmp_dir/shovill_res
mkdir -p $tmp_dir/fastp_res

# declare input file array
declare -a input_array=($INPUT_DIRECTORY/*)
getfilenames input_array

# genome assembly
for i in ${!input_array[@]}; do
    seq_1=$(echo ${input_array[${i}]} | cut -d, -f1)
    seq_2=$(echo ${input_array[${i}]} | cut -d, -f2)
    assembly $seq_1 $seq_2
done

# declare contig file array
declare -a contig_array=($tmp_dir/shovill_res/*)
getcontignames contig_array

# reference sequence comparisons