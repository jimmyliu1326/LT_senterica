#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(furrr))

# parse arguments
# args[1] = abricate results directory
# args[2] = output object path
# args[3] = number of threads
args <- commandArgs(trailingOnly = T)

# set up parallel sessions
plan(multisession, workers = as.numeric(args[3]))

# loading profiles
read_profile <- function(path) {
	
	files <- paste0(path, "/", list.files(path))
	profile <- future_map_dfr(files, function(x) {
		df <- fread(x, sep = "\t", header = T) %>% 
				rename("FILE" = `#FILE`,	
							 "Coverage" = `%COVERAGE`,
							 "Identity" = `%IDENTITY`,
							 "Gene" = GENE) %>% 
				select(FILE, Gene, Coverage, Identity) %>% 
				mutate(FILE = as.character(gsub(".fasta|.fna|.fa", "", FILE)),
							 Gene = as.character(Gene))
			
		# check if profile is empty
		if (nrow(df)==0) { 
			df <- data.frame(`FILE` = gsub(".*/|_vfdb_res|_plasmidfinder_res|.tsv", "", x),
											 "Gene" = "None",
											 "Coverage" = NA,
											 "Identity" = NA,
											 stringsAsFactors = F)
		}
			return(df) 
		}) %>% 
		group_by(`FILE`) %>% 
		nest()
		return(profile)
}

# main
main <- function(db_path, output_path) {
	
	# read data
	database_profile <- read_profile(db_path)
	# write output
	saveRDS(database_profile, file = output_path)
}

# call main
main(args[1], args[2])