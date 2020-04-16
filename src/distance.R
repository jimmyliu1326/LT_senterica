#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(phangorn))

# load allele profile
profile <- fread(args[1], header = T, sep = "\t")

# calculate hamming distance
combinations <- expand_grid(x = profile[[1]], y = profile[[1]])

dist <- map2(.x = combinations$x, .y = combinations$y, function(x,y) {
    x_profile <- profile[which(profile$FILE == x), 2:ncol(profile)]
    y_profile <- profile[which(profile$FILE == y), 2:ncol(profile)]
    hamming_dist <- sum(x_profile != y_profile)
    return(hamming_dist)
  }) %>% unlist()

matrix <- cbind(combinations, dist) %>%
  pivot_wider(names_from = y, values_from = dist) %>%
  column_to_rownames("x") %>% 
  as.matrix()

# write output in phylip format
writeDist(matrix)

# write output in tsv format
#write.table(matrix, quote = F, row.names = T, col.names = NA, sep = "\t")