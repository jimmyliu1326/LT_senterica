#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# initialize threshold
threshold <- as.numeric(args[2])

# identify top hits
top_ref <- fread(args[1], header = F) %>%
    rename("query" = V1,
           "reference" = V2,
           "identity" = V3,
           "shared-hashes" = V4) %>%
    arrange(desc(identity)) %>%
    filter(identity <= threshold) %>%
    pull(reference)

write.table(top_ref, sep = "\t", quote = F, row.names = F, col.names = F)
