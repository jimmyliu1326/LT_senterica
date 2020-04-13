#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# identify top hits
top_ref <- fread(args[1], header = F) %>%
    rename("query" = V1,
           "reference" = V2,
           "identity" = V3,
           "shared-hashes" = V4) %>%
    arrange(desc(identity)) %>%
    slice(1) %>%
    pull(reference)

write.table(top_ref, sep = "\t", quote = F, row.names = F, col.names = F)
