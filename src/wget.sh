#!/usr/bin/env bash

# download genomes
# $1 lines containing genome ID and url separated by ',' (e.g. GCA_0001,https://:)
# $2 output directory
    
filename=$(echo $1 | cut -d',' -f1)
source=$(echo $filename | cut -d'_' -f1)
url=$(echo $1 | cut -d',' -f2)
OUT_DIR=$2

if [[ $source == "GCA" ]]; then
    var=$(echo $url | sed 's/.*GCA_//')
    wget -q -nc $url"/GCA_"${var}"_genomic.fna.gz" -O $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa.gz
    echo $url"/GCA_"${var}"_genomic.fna.gz"
    # check filesize
    #size=$(stat -c %s $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa.gz)

    # retrieve again if failed
    #while [[ $size -eq 0 ]]; do
    #    rm $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa.gz
    #    wget -q -nc $url"/GCA_"$var"_genomic.fna.gz" -O $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa.gz
    #    size=$(stat -c %s $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa.gz)
    #    sleep 2s
    #done
        
    # unzip file
    gunzip $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa.gz
else
    wget -q -nc -c $url -O $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa

    # check filesize
    size=$(stat -c %s $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa)

    # retrieve again if failed
    while [[ $size -eq 0 ]]; do
        rm $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa
        wget -q -nc -c $url -O $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa
        size=$(stat -c %s $OUT_DIR/LT_Senterica_tmp/genomes/$filename.fa)
        sleep 2s
    done
fi

