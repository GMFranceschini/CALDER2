## function to obtain the optimal compartment calling from various resolutions (bin_sizes)

retrieve_opt_domain_calls <- function(save_dir, chrs, bin_sizes, with_ref) {
    `%dopar%` <- foreach::`%dopar%`
    `%do%` <- foreach::`%do%`

    get_opt_index <- function(consist_tab) {
        if (ncol(consist_tab) == 2)
        ## if only one bin_size value
            {
                opt_index <- rep(1, nrow(consist_tab))
                return(opt_index)
            }

        mins <- apply(consist_tab[, -1], 1, max) - 0.05 ## allow at most 0.05 smaller than the maximum
        opt_index <- apply(1 * ((consist_tab[, -1] >= mins) == 1), 1, function(v) {
            min(which(v == 1))
        })
        return(opt_index)
    }

    ################################

    identities <- "bulk"

    {
        consist_tab_li <- foreach::foreach(identity = identities) %do% {
            consist_tab <- foreach::foreach(bin_size = bin_sizes, .combine = merge) %do% {
                bin_size_kb <- sprintf("%skb", bin_size / 1E3)
                save_dir_binsize <- file.path(
                    save_dir,
                    "intermediate_data/sub_compartments",
                    bin_size_kb
                )
                consist_tab_tmp <- data.table::data.table(
                    chr = paste0("chr", chrs), val =
                        0
                )
                consist_tab_tmp$val <- foreach::foreach(chr = chrs, .combine = c) %do% {
                    log_file <- paste0(save_dir_binsize, "/chr", chr, "_log.txt")
                    # print(log_file)
                    cor_val <- as.numeric(strsplit(readLines(log_file)[5], "this chr is:")[[1]][2])
                    # print(cor_val)
                    cor_val
                }
                colnames(consist_tab_tmp)[2] <- bin_size_kb
                consist_tab_tmp
            }

            s <- gsub("chr", "", consist_tab[["chr"]])

            x <-
                suppressWarnings(as.numeric(s)) ## https://stackoverflow.com/questions/70080294/sort-column-in-r-strings-first-alphabetically-then-numbers-numerically
            consist_tab <- consist_tab[order(x, "is.na<-"(s, !is.na(x))), ]
            # print((consist_tab))
        }

        names(consist_tab_li) <- identities

        min_consist <- (sapply(consist_tab_li, function(v) {
            min(v[, 2])
        }))

        ################################ Choose the best bin size (as small as possible, and save the chosen comp to the opt dir)

        {
            save_dir_opt <- file.path(save_dir, "sub_compartments")

            dir.create(save_dir_opt,
                recursive = TRUE,
                showWarnings = FALSE
            )

            consist_tab <- consist_tab_li[[identity]]
            opt_index <- get_opt_index(consist_tab)
            names(opt_index) <- gsub(":", "", consist_tab$chr)
            opt_bin_tab <- data.table::data.table(
                chr = names(opt_index),
                opt_binsize = (colnames(consist_tab)[-1])[opt_index]
            )

            list_doms <- list()
            for (i in 1:nrow(opt_bin_tab)) {
				chrn = opt_bin_tab$chr[i]
                path_topDom <- paste0(save_dir, "intermediate_data/sub_compartments/", opt_bin_tab$opt_binsize[i], "/", chrn, "_topDom.rds")
                td <- readRDS(path_topDom)
                list_doms[[chrn]] <- td
            }

            save_opt_tab_file <- file.path(save_dir_opt, sprintf("TopDom_domains.rds"))
            saveRDS(list_doms, file = save_opt_tab_file)
        }
    }
}
