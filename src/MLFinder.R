#!/usr/bin/env Rscript

# load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(furrr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dbscan))

# parse arguments
args <- commandArgs(trailingOnly = F)
pipeline_dir <- dirname(gsub("--file=", "", args[which(grepl("--file=", args) == T)]))
input_directory <- args[which(grepl("--args", args) == T)+1]
sample_name <- args[which(grepl("--args", args) == T)+2]
output_directory <- args[which(grepl("--args", args) == T)+3]
n_threads <- args[which(grepl("--args", args) == T)+4] %>% as.numeric()
query_amr_res <- args[which(grepl("--args", args) == T)+5]
query_vf_res <- args[which(grepl("--args", args) == T)+6]
query_plasmid_res <- args[which(grepl("--args", args) == T)+7]
annotations <- args[which(grepl("--args", args) == T)+8]

# import annotation distance functions
source(paste0(pipeline_dir, "/annotation_distance.R"))

# set up parallel sessions
options(future.rng.onMisuse="ignore")
plan(multisession, workers = n_threads)

# read and pre-process dataframes
read <- function(filename, kvalue) {
	
	df <- fread(filename, header = F, sep = "\t")
	colnames(df) <- c("reference", "query", "dist", "pval", "jaccard")
	
	df <- df %>% 
		separate(jaccard, into = c("shared_k", "total_k"), sep = "/") %>% 
		mutate(kvalue = kvalue,
					 shared_k = as.numeric(shared_k),
					 total_k = as.numeric(total_k),
					 pmatch = shared_k/total_k) %>% 
		dplyr::select(reference, kvalue, pmatch, dist)
	
	
	return(df)
}

# merge mash results with k=13,15,17,19,21,23
load <- function(dir, sample_name) {
	k <- c(13,17,21,25,29)
	df <- future_map_dfr(k, function(x) read(paste(paste(dir, sample_name, sep = "/"), x, "tsv", sep = "."), x))
	
	# identify irrelevant strains
	ref_ids <- df %>% 
		filter(kvalue == 21, 
					 dist <= 0.005) %>% 
		pull(reference)
	# filter irrelevant strains and create nested dataframes
	df <- df %>% 
		filter(reference %in% ref_ids) %>% 
		dplyr::select(-dist) %>% 
		group_by(reference) %>% 
		nest() %>% 
		mutate(reference = gsub(".*/|.fasta|.fna|.fa", "", reference)) %>%  # remove paths and extensions in reference ids
		ungroup()
	
	return(df)
}

# construct linear fit
linear_fit <- function(df) {
	
	# log transform pmatch
	df <- df %>% 
		mutate(pmatch = log(pmatch))
	# design matrix
	dmatrix <- matrix(df$kvalue, length(df$kvalue), 1)
	# linear fit
	linear.fit <- lm.fit(cbind(1, dmatrix), df$pmatch)
	#linear.fit <- lm(pmatch ~ kvalue, df)
	intercept <- linear.fit$coefficients[[1]]
	slope <- linear.fit$coefficients[[2]]
	
	return(paste(intercept, slope, sep = ","))
}

# compute accessory and core genome distances
compute <- function(df, linear_regression) {
	
	distances <- df %>% 
		cbind(linear_fit = linear_regression) %>% 
		separate(linear_fit, into = c("intercept", "slope"), sep = ",") %>% 
		mutate(intercept = as.numeric(intercept),
					 slope = as.numeric(slope),
					 accessory = 1-exp(1)^intercept,
					 core = 1-exp(1)^slope,
					 core = case_when(core < 0 ~ 0,
					 								 T ~ core),
					 accessory = case_when(accessory < 0 ~ 0,
					 											T ~ accessory))
	
	return(distances)
}

# optimize hdbscan hyperparameters
optimize <- function(df, range) {
	
	optimize_tbl <- tibble(minPts = range,
												 low_proportion = 0)
	optimize_tbl$low_proportion <- future_map_dbl(range, function(x) {
		
		# cluster data
		cluster_res <- hdbscan(df, minPts = x)
		# compute proportion of low confidence assignments (p < 0.05)
		low_count <- length(cluster_res$membership_prob[which(cluster_res$membership_prob < 0.05)])
		total_count <- length(cluster_res$membership_prob)
		low_proportion <- low_count/total_count
		# check if cluster results contain a single cluster to avoid using a given minPts that yield a single cluster
		if (length(unique(cluster_res$cluster)) == 1) {
			low_proportion <- 1
		}
		
		return(low_proportion)
	})
	
	# identify optimal minPts
	opt_minPts <- optimize_tbl %>% 
		filter(low_proportion == min(low_proportion)) %>% 
		arrange(-desc(minPts)) %>% 
		slice(1) %>% 
		pull(minPts)
	
	return(opt_minPts)
}

# plot cluster assignments 
plot_cluster <- function(df, cutoff) {
	
	plot <- ggplot() +
		geom_point(data = df %>%
							 	filter(cluster == 0),
							 mapping = aes(x = core, y = accessory),
							 alpha = 0.1,
							 color = "black") +
		geom_point(data = df %>% 
							 	filter(cluster != "0"),
							 mapping = aes(x = core, y = accessory, color = cluster),
							 alpha = 0.1) +
		#facet_zoom(x = core <= 0.0025, split = T)+
		stat_ellipse(data = df %>%
								 	filter(cluster != "0"),
								 geom = "polygon",
								 alpha = 0.2,
								 level = 0.9,
								 mapping = aes(x = core, y = accessory, color = cluster, fill = cluster)) +
		geom_label(data = df %>%
							 	filter(cluster != "0") %>%
							 	group_by(cluster) %>%
							 	summarize(avg_core = mean(core), # average X and Y to find cluster centroid
							 						avg_accessory = mean(accessory)),
							 mapping = aes(x = avg_core, y= avg_accessory, label = cluster, color = cluster))+
		#geom_abline(slope = cutoff[[1]], intercept = cutoff[[2]], color = "red") +
		#geom_abline(slope = -0.029, intercept = 0.016, color = "red") +
		geom_hline(yintercept = cutoff[[2]], color = "red") +
		geom_vline(xintercept = cutoff[[1]], color = "red") +
		xlab("Core genome distance") +
		ylab("Accessory genome distance") +
		#scale_color_brewer(palette = "Set2") +
		guides(color = F,
					 fill = F) +
		theme_bw() #+
	#theme_minimal() #+
	#scale_x_continuous(limits = c(0,0.005)) +
	#scale_y_continuous(limits = c(0,0.1))
	
	return(plot)
}

# identify threshold cut-offs defined by the cluster closest to origin
core_accessory_threshold <- function(df) {
	
	# determine the two cluster ids closest to origin
	cluster_dist <- df %>% 
		filter(cluster != 0) %>% # remove noise
		group_by(cluster) %>% 
		summarize(avg_core = mean(core), # average X and Y to find cluster centroid
							avg_accessory = mean(accessory)) %>% 
		mutate(distance = sqrt(avg_core^2 + avg_accessory^2))
	
	cluster_id <- cluster_dist %>% 
		arrange(-desc(distance)) %>% 
		slice(1:2) %>% 
		pull(cluster)
	
	# calculate slope, midpoint and intercept between the two closest clusters
	# cluster1 <- as.numeric((cluster_dist %>% filter(cluster == cluster_id[1]) %>% select(avg_core, avg_accessory))[1,])
	# cluster2 <- as.numeric((cluster_dist %>% filter(cluster == cluster_id[2]) %>% select(avg_core, avg_accessory))[1,])
	# 
	# slope <- ((cluster2[2]-cluster1[2])/(cluster2[1]-cluster1[1]))
	# midpoint <- c((cluster2[1]+cluster1[1])/2, (cluster2[2]+cluster1[2])/2)
	# print(midpoint)
	# intercept <- midpoint[2]-slope*midpoint[1]
	
	#return(list(slope, intercept))
	
	# determine max core and accessory boundaries of cluster closest to origin
	max_core <- df %>% 
		filter(cluster == cluster_id[1]) %>% 
		pull(core) %>% 
		max()
	
	max_accessory <- df %>% 
		filter(cluster == cluster_id[1]) %>% 
		pull(accessory) %>%   
		max()
	
	return(list(c(max_core, max_accessory), cluster_id[1]))
}

full_dimensional_threshold <- function(df, features) {
	
	# determine the two cluster ids closest to origin
	cluster_dist <- df %>% 
		filter(cluster != 0) %>% # remove noise
		group_by(cluster) %>% 
		summarize(core = mean(core)^2, # average X and Y to find cluster centroid
							accessory = mean(accessory)^2,
							#amr_dist = mean(amr_dist)^2,
							vf_dist = mean(vf_dist)^2,
							plasmid_dist = mean(plasmid_dist)^2) %>% 
		mutate(distance = rowSums(.[names(.) %in% features]))
	
	cluster_id <- cluster_dist %>% 
		arrange(-desc(distance)) %>% 
		slice(1:2) %>% 
		pull(cluster)
	
	boundaries <- future_map(features, function(x) {
		df %>% filter(cluster == cluster_id[1]) %>% pull(x) %>% max()
	})
	names(boundaries) <- features
	
	return(list(boundaries, cluster_id[1]))
}

# filter entire database given cutoff
database_filter <- function(df, cutoff, features) {
	filter_exp <- paste(features, "<= ", cutoff)
	neighbours <- df %>% 
		filter(!!!parse_exprs(filter_exp))
	return(neighbours)
}

# identify reference lineages
identify_reference <- function(df, origin_cluster, features) {
	# filter sequences of nearest origin cluster and noise sequences
	df <- df %>% 
		dplyr::select(any_of(c(features, "cluster", "reference"))) %>% 
		filter(cluster != origin_cluster,
					 cluster != "0") %>% 
		group_by(cluster) %>% 
		nest()
	# calculate cluster centroids
	centroids <- future_map(features, function(x) future_map_dbl(df$data, ~mean(.[[x]]))) %>% set_names(features)
	# calculate distance of each sequence to the cluster centroid
	# and select the sequence with the least distance to each cluster centroid
	df[["data"]] <- future_map2(df$data, seq_along(1:nrow(df)), function(x,y) {
		future_map(features, function(z) {
			(x[[z]]-centroids[[z]][[y]])^2
		}) %>% 
			set_names(paste0(features,"_2")) %>% 
			bind_cols(x, .) %>% 
			mutate(distance = sqrt(rowSums(.[features]))) %>% 
			filter(distance == min(distance)) %>% 
			slice(1)
	}) %>% unnest()
	
	return(df)
}

# write out list of neighbours/references
list_output <- function(df, output_path) {
	# load database ftp paths
	ftp_paths <- fread(paste0(pipeline_dir, "../database/Senterica_ftp_paths.tsv"), sep = "\t", header = T)
	# filter ftp paths
	output <- ftp_paths %>% filter(asm_acc %in% df$reference)
	# write output
	write.table(output, file = output_path, quote = F, col.names = F, row.names = F)
}

# main
main <- function(dir, sample_name, out_dir, query_amr, query_vfdb, query_plasmid, annotations) {
	# load data
	mash_res <- load(dir, sample_name)
	# load genome annotation database
	#amr_db <- readRDS(paste0(pipeline_dir, "/../database/", "Senterica_ref_amrfinder.RDS"))
	vf_db <- readRDS(paste0(pipeline_dir, "/../database/", "Senterica_ref_vfdb.RDS"))
	plasmid_db <- readRDS(paste0(pipeline_dir, "/../database/", "Senterica_ref_plasmidfinder.RDS"))
	# compute annotation distances
	mash_res <- mash_res %>% mutate(#vf_dist = future_map_dbl(reference, ~annotations_distance(read_query_vfdb(query_vfdb), vf_db, .)),
									plasmid_dist = future_map_dbl(reference, ~annotations_distance(read_query_plasmidfinder(query_plasmid), plasmid_db, .)))
									#amr_dist = future_map_dbl(reference, ~annotations_distance(read_query_amrfinder(query_amr), amr_db, .)))
	# fit linear model
	linear.fit <- future_map_chr(mash_res$data, ~linear_fit(.))
	# compute accessory and core genome distances
	mash_res <- compute(mash_res, linear.fit)
	message(as.character(nrow(mash_res)))
	# identify optimal hyperparameters for hdbscan
	opt_minpts <- optimize(mash_res[c("core", "accessory")], seq(50, 100, 10))
	opt_minpts <- optimize(mash_res[c("core", "accessory")], seq(opt_minpts-9, opt_minpts+9, 1))
	# cluster data using optimal parameters
	cluster_res <- hdbscan(mash_res[c("core", "accessory")], minPts = opt_minpts)
	mash_res <- cbind(mash_res, cluster = as.character(cluster_res$cluster))
	# identify core-accessory threshold cut-off and nearest origin cluster code
	cutoff <- core_accessory_threshold(mash_res)[[1]]
	names(cutoff) <- c("core", "accessory")
	origin_cluster <- as.character(core_accessory_threshold(mash_res)[[2]])
	# plot cluster assignments
	p <- plot_cluster(mash_res, cutoff)
	ggsave(paste0(out_dir, "/plots/", sample_name, "_cluster_res.png"), p, width = 10, height = 10)
	# filter database
	mash_res_filter <- database_filter(mash_res, cutoff, c("core", "accessory"))
	# output neighbours if excluding annotation variations
	if (annotations == "0") {
		# identify reference lineages
		reference_lineages <- identify_reference(mash_res, origin_cluster, c("core", "accessory"))
		# write out neighbour and reference ids
		list_output(mash_res_filter, paste0(out_dir, "/neighbours/", sample_name, ".list"))
		list_output(reference_lineages, paste0(out_dir, "/references/", sample_name, ".list"))
		# exit script prematurely
		stop()
	}
	# define features for secondary clustering
	features <- c("core", "plasmid_dist")
	# further partition the nearest origin cluster by amr, vf and plasmid distances
	if (nrow(mash_res_filter) < 100) {
		opt_range <- seq(5,nrow(mash_res_filter),10)
	} else {
		opt_range <- seq(5,100,10)
	}
	opt_minpts <- optimize(mash_res_filter[features], opt_range)
	if (opt_minpts == 5) {
		opt_range <- seq(2,14,1)
	} else {
		opt_range <- seq(opt_minpts-9,opt_minpts+9,1)
	}
	opt_minpts <- optimize(mash_res_filter[features], opt_range)
	cluster_res <- hdbscan(mash_res_filter[features], minPts = opt_minpts)
	mash_res_filter <- cbind(mash_res_filter %>% dplyr::select(-cluster), cluster = as.character(cluster_res$cluster))
	# identify cutoff and nearest origin cluster by defined features
	full_threshold_res <- full_dimensional_threshold(mash_res_filter, features)
	full_threshold <- full_threshold_res[[1]]
	full_origin_cluster <- full_threshold_res[[2]] %>% as.character()
	# identify reference lineages
	reference_lineages <- identify_reference(mash_res, origin_cluster, c("core", "accessory"))
	# write out neighbour and reference ids
	list_output(neighbours, paste0(out_dir, "/neighbours/", sample_name, ".list"))
	list_output(reference_lineages, paste0(out_dir, "/references/", sample_name, ".list"))
}

# call main
main(input_directory, sample_name, output_directory, query_amr_res, query_vf_res, query_plasmid_res, annotations)