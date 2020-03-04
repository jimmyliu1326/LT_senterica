#!/usr/bin/env bash

usage() { 
 echo "
Usage: $0
 Required:
 -1       Input forward fastq
 -2       Input reverse fastq
 -o       specify output directory
 
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

opts=`getopt -o h1:2:o:t: -l qc,mlst,amr,help \
      -- "$@"`
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