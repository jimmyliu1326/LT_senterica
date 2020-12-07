# load amrfinder query results
read_query_amrfinder <- function(path) {
		profile <- fread(path, sep = "\t", header = T) %>% 
			filter(`Element type` == "AMR") %>% 
			rename("FILE" = "Name",
						 "Gene" = `Gene symbol`,
						 "Coverage" = `% Coverage of reference sequence`,
						 "Identity" = `% Identity to reference sequence`) %>% 
			select(Gene, Coverage, Identity)
		# check if profile is empty
		if (length(profile) == 0) { 
			profile <- data.frame("FILE" = "query",
														"Gene" = "None", stringsAsFactors = F,
														"Coverage" = NA,
														"Identity" = NA)
		}
		return(profile)
}

# load plasmidfinder query results
read_query_plasmidfinder <- function(path, type) {
	
		profile <- fread(path, sep = "\t", header = T) %>% 
			rename("FILE" = "#FILE",
						 "Gene" = `GENE`,
						 "Coverage" = `%COVERAGE`,
						 "Identity" = `%IDENTITY`) %>% 
			select(Gene, Coverage, Identity)
		# check if profile is empty
		if (length(profile) == 0) { 
			profile <- data.frame("FILE" = "query",
														"Gene" = "None", stringsAsFactors = F,
														"Coverage" = NA,
														"Identity" = NA)
		}
		return(profile)
}

# load VFDB query results
read_query_vfdb <- function(path) {
	
	profile <- fread(path, sep = "\t", header = T) %>% 
		rename("FILE" = "#FILE",
					 "Gene" = `GENE`,
					 "Coverage" = `%COVERAGE`,
					 "Identity" = `%IDENTITY`) %>% 
		select(Gene, Coverage, Identity)
	# check if profile is empty
	if (length(profile) == 0) { 
		profile <- data.frame("FILE" = "query",
													"Gene" = "None", stringsAsFactors = F,
													"Coverage" = NA,
													"Identity" = NA)
	}
	return(profile)
}

# compare hits shared between query and reference
compare_hits <- function(query_df, reference_df) {
	# find shared genes
	shared_genes <- intersect(query_df$Gene[which(query_df$Gene != "None")], reference_df$Gene[which(reference_df$Gene != "None")])
	# if no shared genes, return zero distance
	if (length(shared_genes) == 0) {
		distance <- 0
	} else {
		distance <- future_map_dbl(shared_genes, function(x, query=query_df, reference=reference_df) {
			# subset by given gene
			query <- query %>% filter(Gene == x) %>% arrange(desc(Coverage), desc(Identity)) %>% slice(1)
			reference <- reference %>% filter(Gene == x) %>% arrange(desc(Coverage), desc(Identity)) %>% slice(1)
			# calculate difference in coverage and identity
			coverage_diff <- abs(query$Coverage-reference$Coverage)/100
			identity_diff <- abs(query$Identity-reference$Identity)/100
			# compute gene distance
			sum_diff <- 0.5*(coverage_diff+identity_diff)
			return(sum_diff)
		}) %>% sum()
	}
	return(distance)
}

# compute resistome distance relative to query
annotations_distance <- function(query_profile, database_profile, reference_target) {
	
	# identify pan-resistome
	total_loci <- length(unique(c(query_profile$Gene[which(query_profile$Gene != "None")], unlist(future_map(database_profile$data, ~unique(.$Gene[which(.$Gene != "None")]))))))
	# filter database_profile by reference_target
	target_profile <- database_profile %>% filter(FILE == reference_target) %>% pull(data)
	# compute amr distance
	unique_database <- length(setdiff(target_profile[which(target_profile$Gene != "None")], query_profile$Gene[which(query_profile$Gene != "None")]))
  unique_query <- length(setdiff(query_profile$Gene[which(query_profile$Gene != "None")], target_profile$Gene[which(target_profile$Gene != "None")]))
	shared_gene <- compare_hits(query_profile, target_profile)
	distance <- (unique_database+unique_query+shared_gene)/total_loci
	# return resistome distance vector
	return(distance)
}