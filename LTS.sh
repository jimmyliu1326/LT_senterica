#!/usr/bin/env bash
#################
# usage message #
#################
usage() {
    echo "
Usage: $0
 Available modules:
 1. molecularlinkage    identify closely related Salmonella genomes in PubMLST and GenomeTrakr using a set of query reads or assemblies
 2. download_db         download reduced Salmonella genome representation database and functional annotation profiles
 3. update_db           update existing Salmonella database
"
}
molecularlinkage_usage() { 
 echo "
Usage: $0 molecularlinkage
 Required arguments:
 -i|--input            Input directory containing raw reads [Must specify either raw reads or assembled contigs using -i or -l]
 -l|--list             List of query genome paths [Must specify either raw reads or assembled contigs using -i or -l]
 -o|--output           Output directory
 
 Optional arguments:
 --notree              Find genomic neighbours in public genome databases for each sample only, skip phlyogenetic tree construction
 --noannotations       Do not use functional variations as clustering dimensions, resulting in lower discriminatory power
 -t|--threads          Number of threads [Default: 32]
 -h|--help             Display help message
"
}
download_usage() {
    echo "
Usage: $0 download_db
    Description:
        Running $0 download_db will retrieve all the required data used by the pipeline to path where the pipeline package is stored
    
    Optional arguments:
    -h|--help           Display help message
"
}
update_usage() {
    echo "
Usage: $0 update_db
    Description:
        Running $0 update_db will update existing database information at where the pipeline package is stored
    
    Optional arguments:
    -h|--help           Display help message
"
}
#############################
# define all variables used #
#############################
n_threads=1
pipeline_dir=$(dirname $(realpath $0))
output_dir=""
input_dir=""
list=""
tree=1
annotations=1
###################
# parse arguments #
###################
# exit if no arguments given
if [ $# == 0 ]; then usage; exit 0; fi
# exit if invalid modules given
if [[ ! $1 =~ ^(molecularlinkage|download_db|update_db)$ ]]; then
    echo "Invalid module name given"
    usage
    exit 1
fi
# parse subcommands
if [[ $1 == "molecularlinkage" ]]; then
    shift 1
    # print help message if no arguments passed and exit
    if [[ $# -eq 0 ]]; then molecularlinkage_usage; exit 0; fi
    # evaluate module arguments passed
    opts=`getopt -o hi:o:t:s:l: -l help,input:,output:,threads:,list:,notree,noannotations -- "$@"`
    if [ $? != 0 ] ; then echo "WARNING: Invalid arguements used, exiting"; molecularlinkage_usage; exit 1 ; fi
    eval set -- "$opts"

    while true; do
        case "$1" in
            -i|--input) input_dir=$(realpath $2); shift 2 ;;
            -l|--list) list=$2; shift 2 ;;
            -o|--output) output_dir=$(realpath $2); mkdir -p $output_dir; shift 2 ;;
            -t|--threads) n_threads=$2; shift 2 ;;
            --notree) tree=0; shift 1 ;;
            --noannotations) annotation=0; shift 1 ;;
            --) shift; break ;;
            -h|--help) molecularlinkage_usage; exit 0;;
        esac
    done
    # evaluate if all required arugments are present
    if [[ $input_dir == "" && $list == "" ]]; then echo "Must specify either raw reads or assemblies using -i or -l, exiting"; exit 1; fi
    if [[ $output_dir == "" ]]; then echo "The output directory was not specified using -o, exiting"; exit 1; fi
    # validate input raw read directory
    if ! test -d $input_dir; then echo "$input_dir does not exist or is not a directory, exiting"; exit 1; fi
    # validate input raw read files
    if [[ $input_dir != "" ]]; then
        for samples_dir in $input_dir; do
            if ! test -d $samples_dir; then echo "${samples_dir} is not a directory that contains fastq files. Remove it and run again"; exit 1; fi
            if [[ $(ls $samples_dir/*R1*.fastq | wc -l) -ne 1 ]]; then echo "${sample_dir} does not contain a single FORWARD read fastq file."; exit 1; fi
            if [[ $(ls $samples_dir/*R2*.fastq | wc -l) -ne 1 ]]; then echo "${sample_dir} does not contain a single REVERSE read fastq file."; exit 1; fi
        done
    fi
    # validate input genome list
    if ! test -f $list; then echo "$list does not exiting, exiting"; exit 1; fi
    # validate genome files in input list
    while read lines; do
        if ! test -f $lines; then echo "$lines cannot be found, verify its path listed in ${list}, exiting"; exit 1; fi
    done < $list
    # validate if genome database is present
    for i in 13 17 21 25 29; do
        if ! test -f $pipeline_dir/database/Senterica_ref.${i}.msh; then
            echo "$pipeline_dir/database/Senterica_ref.${i}.msh cannot be found, download the required database files before running the pipeline, exiting"
            exit 1
        fi
    done
    # validate if genome annotation profiles are present
    for i in "vfdb" "plasmidfinder"; do
        if ! test -f $pipeline_dir/database/Senterica_ref_${i}.RDS; then 
        echo "$pipeline_dir/database/Senterica_ref_${i}.RDS cannot be found, download the required database files before running the pipeline, exiting"
        fi
    done
    # create snakemake samples input file
    echo "Sample,R1,R2,Type" > $output_dir/samples.csv # remove if samples.csv file already exists
    if [[ $input_dir != "" ]]; then # identify samples in input directory if given and add R1/R2 absolute paths to samples.csv
        for samples_dir in $input_dir; do
            sample=$(basename $samples_dir)
            R1=$(realpath $(ls $samples_dir/*R1*))
            R2=$(realpath $(ls $samples_dir/*R2*))
            type="raw"
            echo "$(basename $samples_dir),${R1},${R2},type" >> $output_dir/samples.csv
        done
    fi
    if [[ $list != "" ]]; then # read genome list if given and add absolute paths to samples.csv
        while read lines; do
            sample=$(basename $lines)
            R1=$(realpath $lines)
            R2=""
            type="assembly"
            echo "${sample%.*},${R1},${R2},$type" >> $output_dir/samples.csv
        done < $list
    fi
    # call snakemake
    snakemake --cores $n_threads --snakefile $pipeline_dir/Snakefile \
        --config samples=$output_dir/samples.csv output_dir=$output_dir threads=$n_threads pipeline_dir=$pipeline_dir tree=$tree annotations=$annotations
elif [[ $1 == "download_db" ]]; then
    shift 1
    # evaluate module arguments passed
    opts=`getopt -o h -l help -- "$@"`
    if [ $? != 0 ] ; then echo "WARNING: Invalid arguements used, exiting"; download_usage; exit 1 ; fi
    eval set -- "$opts"

    while true; do
        case "$1" in
            -h|--help) download_usage; exit 0 ;;
            --) shift; break ;;
        esac
    done

    # return
    echo "This function has yet to be implemented"
elif [[ $1 == "update_db" ]]; then
    shift 1
    # evaluate module arguments passed
    opts=`getopt -o h -l help -- "$@"`
    if [ $? != 0 ] ; then echo "WARNING: Invalid arguements used, exiting"; update_usage; exit 1 ; fi
    eval set -- "$opts"

    while true; do
        case "$1" in
            -h|--help) update_usage; exit 0 ;;
            --) shift; break ;;
        esac
    done

    # return
    echo "This function has yet to be implemented"
fi