#!/usr/bin/env bash

usage() { 
 echo "
Usage: $0
 Required:
 -1       Input forward fastq
 -2       Input reverse fastq
 -o       specify output directory
 -tmp     specify directory for storing temporary files
 
 Options:
 -t       Number of threads [Default: 1]
 -h       Display help message

 "
}

# initialize variables
n_threads=1

# parse arguments
if [ $# == 0 ]
then
    usage
    exit
fi

opts=`getopt -o h1:2:o:t: -l tmp,help -- "$@"`
if [ $? != 0 ] ; then echo "WARNING: Invalid arguements used, exiting"; usage; exit 1 ; fi
eval set -- "$opts"

while true; do
  case "$1" in
    -1) fastq_1=$2; shift 2 ;;
    -2) fastq_2=$2; shift 2 ;;
    -o) OUT_DIR=$2; shift 2 ;;
    -t) n_threads=$2; shift 2 ;;
    --) shift; break ;;
    -h|--help) usage; exit ;;
  esac
done

fastp_exe() {
  time=$(date +"%T")
  echo "[$time] Read trimming using fastp"

  source /opt/galaxy/tool_dependencies/_conda/bin/activate /home/$USER/.conda/envs/fastp
  cleaned_fastq_1=/scratch/$USER/tmp/$(basename ${fastq_1%.*})_clean.fastq
  cleaned_fastq_2=/scratch/$USER/tmp/$(basename ${fastq_2%.*})_clean.fastq
  fastp -i $fastq_1 -I $fastq_2 \
        -o $cleaned_fastq_1 \
        -O $cleaned_fastq_2 \
        -w $n_threads &> $OUT_DIR/log_files/fastp.log
}