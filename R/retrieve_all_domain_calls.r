retrieve_all_domain_calls = function(save_dir, chrs, bin_sizes){
	## Retrieve all domain calls for each chromosome and bin size
	for(bin_size in bin_sizes){
	bin_size_kb = paste0(bin_size/1E03, "kb")

	domain_bin = list()
	for(chr in chrs){
			path_topDom <- paste0(save_dir, "/intermediate_data/sub_compartments/", bin_size_kb, "/chr", chr,"_", bin_size_kb, "_TopDom.rds")
			td <- readRDS(path_topDom)
			domain_bin[[chr]] = td
		}
		saveRDS(domain_bin, file = paste0(save_dir, "sub_compartments/TopDom_domains_", bin_size_kb, ".rds"))
	}
}