#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# load salmonella genome cluster and ftp_path info
cluster <- fread("data/ftp_path.tsv")

# load top hits
top_hit <- fread(args[1], header = F) %>%
    rename(asm_acc = V1)

# clean top hits
top_hit$asm_acc <- gsub(".*/|\\..*","", top_hit$asm_acc)

# identify cluster codes of top hits
target_clusters <- cluster %>%
    filter(asm_acc %in% top_hit$asm_acc) %>%
    pull(ClusterNumber)

# identify all candidates in same clusters as top hits
candidates <- cluster %>%
    filter(ClusterNumber %in% target_clusters)

# write output
write.table(candidates[c("asm_acc", "ftp_path")], quote = F, row.names = F, col.names = F, sep = ",")