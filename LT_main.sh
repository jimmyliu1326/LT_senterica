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
IN_DIR=""
OUT_DIR=""
ref_seq=""

# parse arguments
if [ $# == 0 ]
then
    usage
    exit
fi

opts=`getopt -o hi:o:t:s: -l tmp:,help,input:,output:,sketch: -- "$@"`
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


getreadnames() {
    input_array=($@)
    for i in ${!input_array[*]}; do
        seq_1=$(echo ${input_array[${i}]}/*R1*)
        seq_2=$(echo ${input_array[${i}]}/*R2*)
        input_array[${i}]=$(echo ${seq_1},${seq_2})
    done
}

getcontignames() {
    contig_array=($@)
    for i in ${!contig_array[*]}; do
        contig_array[${i}]=$(echo ${contig_array[${i}]},${contig_array[${i}]}"/contigs.fa")
    done
}

getfilename_single() {
    for i in "${!$1[*]}"; do
        $1[${i}]=$(echo $1[${i}]/*.$2)
    done
}

getfilename_multi() {
    for i in "${!$1[*]}"; do
        $1[${i}]=$(echo $1[${i}]/*.$2)
    done
}

# Genome assembly
assembly() {  
  filename=$(basename $(dirname $1))
  
  # Read trimming using fastp
  time=$(date +"%T")
  echo "[$time] Read trimming on "$(basename $1)" and "$(basename $2)

  source /opt/galaxy/tool_dependencies/_conda/bin/activate /home/$USER/.conda/envs/fastp 

  cleaned_fastq_1=$tmp_dir/fastp_res/$(basename ${1%.*})_clean.fastq
  cleaned_fastq_2=$tmp_dir/fastp_res/$(basename ${2%.*})_clean.fastq
  fastp -i $1 -I $2 \
        -o $cleaned_fastq_1 \
        -O $cleaned_fastq_2 \
        -w $n_threads &>> $tmp_dir/log_files/fastp.log
  
  # Genome assembly using shovill
  time=$(date +"%T")
  echo "[$time] Genome assembly on "$(basename $1)" and "$(basename $2)

  source /opt/galaxy/tool_dependencies/_conda/bin/activate /opt/miniconda2/envs/shovill-1.0.4

  mkdir -p $tmp_dir/shovill_res/$filename
  shovill --outdir $tmp_dir/shovill_res/$filename \
          --R1 $cleaned_fastq_1 \
          --R2 $cleaned_fastq_2 \
          --gsize 4.5M \
          --cpus $n_threads \
          --force &>> $tmp_dir/log_files/shovill.log
}

# Reference sequence comparisons using mash
mash_exe() {
    time=$(date +"%T")
    echo "[$time] Reference sequence comparisons using mash"

    source /opt/galaxy/tool_dependencies/_conda/bin/activate /opt/miniconda2/envs/mash-2.1

    mash sketch -p $n_threads -s 10000 -o $tmp_dir/mash_res/1/sketch/$1.msh $2
    mash dist $tmp_dir/mash_res/1/sketch/$1.msh $ref_seq -p $n_threads > $tmp_dir/mash_res/1/results/$1.tab
}

# process reference sequence comparison mash results
process_ref_seq_mash() {
    time=$(date +"%T")
    echo "[$time] Processing mash results"

    source /opt/galaxy/tool_dependencies/_conda/bin/activate /home/$USER/.conda/envs/r_env

    # identify top hits
    src/top_hit.R $1 > $tmp_dir/mash_res/1/process/$(basename ${1%.*})_top_ref.tsv

    # identify candidates
    src/identify_candidates.R $tmp_dir/mash_res/1/process/$(basename ${1%.*})_top_ref.tsv > $tmp_dir/mash_res/1/process/$(basename ${1%.*})_candidates_ftp_path.csv
}

# identify all unique candidate genomes for download
unique_genomes() {
    cat $tmp_dir/mash_res/1/process/* | sort | uniq -u > $tmp_dir/candidate_genome_list.csv
}

# download genomes
download() {
    filename=$($1 | cut -d, -f1)
    url=$($1 | cut -d, -f2)
    wget -nc $url -O $tmp_dir/genomes/$filename.fa
}

# candidate sequence comparisons


# MAIN

# create tmp directory structure
mkdir -p $tmp_dir/log_files
mkdir -p $tmp_dir/fastp_res
mkdir -p $tmp_dir/shovill_res
mkdir -p $tmp_dir/mash_res/1/sketch
mkdir -p $tmp_dir/mash_res/1/results
mkdir -p $tmp_dir/mash_res/1/process
mkdir -p $tmp_dir/mash_res/2/sketch
mkdir -p $tmp_dir/mash_res/2/results
mkdir -p $tmp_dir/mash_res/2/process
mkdir -p $tmp_dir/genomes

# declare input file array
declare -a input_array=($IN_DIR/*)
getreadnames ${input_array[@]}

# genome assembly
for i in ${!input_array[*]}; do
    seq_1=$(echo ${input_array[${i}]} | cut -d, -f1)
    seq_2=$(echo ${input_array[${i}]} | cut -d, -f2)
    assembly $seq_1 $seq_2
done

# declare contig file array
declare -a contig_array=($tmp_dir/shovill_res/*)
getcontignames ${contig_array[@]}

# reference sequence comparisons
for i in ${!contig_array[*]}; do
    filename=$(basename $(echo ${contig_array[${i}]} | cut -d, -f1))
    filepath=$(echo ${contig_array[${i}]} | cut -d, -f2)
    mash_exe $filename $filepath
done

# declare reference sequence mash results array
declare -a ref_seq_mash_array=($tmp_dir/mash_res/1/results/*)

# process reference sequence comparisons mash results
for i in ${!ref_seq_mash_array[*]}; do
    process_ref_seq_mash ${ref_seq_mash_array[${i}]}
done

# download candidate genomes
unique_genomes

while read lines; do
    download $lines
done < $tmp_dir/candidate_genome_list.csv


