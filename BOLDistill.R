# Description: This script consolidates taxonomy within each BIN/OTU and works in conjunction with BOLDistill.sh
# Authors: Sean Prosser & Spencer Monckton (July 2025)
# Version: 1.4

library(dplyr)
library(data.table)
library(Biostrings)
library(ggplot2)
library(tidyr)
library(patchwork)


# get argument from bash script
args <- commandArgs(trailingOnly = TRUE)
fasta_file <- args[1]
wd <- args[2]
setwd(wd)

# import curated BOLD List
input <- fread("combined_boldlist.tsv", header = T, sep = "\t", fill = TRUE)

# remove unwanted columns
df <- input[,c("bin", "kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "species", "subspecies")]

# remove BINs (these will be curated from original input boldlist)
df <- df[!grepl("BOLD:", df$bin),]

# import original boldlist
input_files <- list.files(pattern = "^boldlist_INPUT_.*\\.tsv$")
if (length(input_files) > 0) {
  boldlist_file <- input_files[1]
} else if (file.exists("boldlist.tsv")) {
  boldlist_file <- "boldlist.tsv"
} else {
  stop("No boldlist input file found.")
}
boldlist <- fread(boldlist_file, header = TRUE, sep = "\t", fill = TRUE)

# reduce boldlist to only records with BINs and remove unwanted columns
boldlist_w_bins <- boldlist[grepl("BOLD:", boldlist$bin),]
boldlist_w_bins <- boldlist_w_bins[,c("bin", "kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "species", "subspecies")]

# merge OTUs and BINs into one table
merged <- rbind(boldlist_w_bins, df)

#####################################################################################################################################

get_bin_consensus <- function(
    df,
    ranks = c("kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "species", "subspecies"),
    threshold = 1.0,
    min_ids = 2,
    enforce_scientific = TRUE,
    groups = "bin_uri",
    discord_format = c("list", "text")) {

  stopifnot("One or more provided `ranks` is/are missing from `df`." = all(ranks %in% names(df)),
            "Provided `groups` column is missing from `df`." = (groups %in% names(df)),
            "`threshold` value(s) must be one or more real numbers (i.e. doubles) between 0 and 1." = is.double(unlist(threshold)) & all(unlist(threshold) >= 0) & all(unlist(threshold) <= 1),
            "`threshold` must be either a single number, a vector of unnamed numbers equal in length to `ranks`, or a named list or vector of numbers with names corresponding to ranks." = ((length(threshold) == 1) | (length(threshold) == length(ranks)) | (!is.null(names(threshold)))),
            "`min_ids` value(s) must be one or more whole numbers greater than zero." = is.numeric(unlist(min_ids)) & all(unlist(min_ids) > 0) & all(unlist(min_ids) %% 1 == 0),
            "`min_ids` must be either a single number, a vector of unnamed numbers equal in length to `ranks`, or a named list or vector of numbers with names corresponding to ranks." = ((length(min_ids) == 1) | (length(min_ids) == length(ranks)) | (!is.null(names(min_ids)))),
            '`discord_format` must be one of "list" or "text".' = all(discord_format %in% c("list", "text")))

  # Define regex for non-scientific names
  re_int <- "\\.\\Z|\\S{2,}\\.\\S|[0-9]|\\s[A-Z]|[A-Z]\\Z|[A-Z]{2}|[a-z][A-Z]|_(?!(hn|sl|ss)\\Z)|%|\\?|!|\\[|\\]|\\{|\\}|\\(|\\)|,|\\s(?:aff|agg|cf|complex|group|grp|gr|gp|cmplx|pr|ms|cfr|nr|nsp|near|nomen|hybrid|voucher|form|from|ss|ssl|see|spp?|sample)\\.?(?:\\s|\\Z)"
  
  # Parse threshold & min_ids parameters and align them with ranks
  parse_param_vector <- function(param) {
    if((length(param) != 1) | !is.null(names(param))) {
      if(is.null(names(param))) {
        param <- unlist(unname(param))
      } else {
        named <- as.list(param[(names(param) %in% ranks) & (!duplicated(param))])
        default <- ifelse("default" %in% names(param), param[["default"]], max(unlist(param)))
        if((!length(named) %in% c(0, length(ranks))) & (!"default" %in% names(param))) {
          warning(paste0("Only some ranks found among `",substitute(param),"` values, with no default given; highest value applied to all unspecified ranks."))
        }
        param <- rep(default, length(ranks))
        for(r in names(named)) param[match(r, ranks)] <- named[[r]]
      }
    } else {
      param <- rep(unlist(param), length(ranks))
    }
    return(param)
  }

  threshold <- parse_param_vector(threshold)
  min_ids <- parse_param_vector(min_ids)

  # Create a copy of the data to avoid mutating by reference
  dt <- as.data.table(copy(df))

  # Replace NA in taxonomy columns with empty values (if ignoring non-scientific names, replace those too)
  if(enforce_scientific) {
    dt[, (c(ranks)) := lapply(.SD, function(x) data.table::fifelse(is.na(x), "", data.table::fifelse(grepl(re_int, x, perl = TRUE), "", as.character(x)))), .SDcols = c(ranks)]
  } else {
    dt[, (c(ranks)) := lapply(.SD, function(x) data.table::fifelse(is.na(x), "", as.character(x))), .SDcols = c(ranks)]
  }

  # Convert data table to matrix for faster row access
  mat <- as.matrix(dt[, c(groups,ranks), with = FALSE])

  # Replace trailing "" with NA so that they are not counted as alternative names
  for (i in seq_len(nrow(mat))) {
    row_vals <- mat[i, ]
    non_blank_idx <- which(row_vals != "")
    if (length(non_blank_idx) > 0) {
      last <- max(non_blank_idx)
      if (last < ncol(mat)) {
        mat[i, (last + 1):ncol(mat)] <- NA_character_
      }
    } else {
      mat[i, ] <- NA_character_  # Entire row is blank
    }
  }

  # Convert back to data.table and restore column names
  dt <- as.data.table(mat)
  setnames(dt, c(groups,ranks))

  # Core consensus logic
  get_consistent_taxon <- function(sub_dt,
                                   ranks = c("kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "species", "subspecies"),
                                   threshold = 1.0,
                                   min_ids = 2) {
  
    id_hier <- sapply(ranks,function(x) NULL)
    concordant = FALSE
    rank_set <- ranks
  
    # Ensure min_ids does not exceed group size
    if (any(min_ids > nrow(sub_dt))) {
      for(i in seq_along(min_ids)) min_ids[[i]] <- nrow(sub_dt)
    }
  
    # Expand threshold and min_ids parameters into full vectors if applicable
    if (length(threshold) == 1) { threshold <- rep(threshold, length(ranks)) }
    if (length(min_ids) == 1) { min_ids <- rep(min_ids, length(ranks)) }
  
    result <- list(
      member_count = nrow(sub_dt),
      concordant_rank = NA_character_,
      concordant_id = NA_character_,
      concordant_id_count = 0L,
      discordant_rank = NA_character_,
      discordant_ids = list(),
      discordant_id_count = 0L
    )
  
    for (rank_col in rev(ranks)) {  # Step backwards through ranks
  
      rank_threshold <- threshold[which(ranks==rank_col)]
      rank_min_ids <- min_ids[which(ranks==rank_col)]
      vals <- sub_dt[[rank_col]]
      filtered <- table(vals[!is.na(vals)])
      name_vals <- proportions(filtered)
      props <- proportions(filtered)[(proportions(filtered) >= rank_threshold) & (filtered >= rank_min_ids)]
      names(name_vals) <- sub("^$","<None>",names(name_vals))
  
      if ((length(props) != 1) & (length(unique(filtered)) > 0)) {
  
        concordant <- FALSE
  
        if(id_hier[rank_col] != "") {
          rank_set <- ranks[0:(which(ranks==rank_col)-1)]
          id_hier <- id_hier[rank_set]
        }
  
        result$discordant_rank <- rank_col
        result$discordant_ids <- list(stats::setNames(as.vector(name_vals), names(name_vals)))
        result$discordant_id_count <- sum(filtered)
  
      } else if ((length(props) == 1) && (names(props)[1] != "")) {
  
        if(!concordant) {
          result$concordant_rank <- rank_col
          result$concordant_id <- names(props)[1]
          result$concordant_id_count <- unname(filtered[names(props)[1]])
        }
  
        concordant <- TRUE
  
        if (is.null(id_hier[[rank_col]]) || is.na(id_hier[[rank_col]]) || is.na(names(props)[1]) || (names(props)[1] != id_hier[[rank_col]])) {
          rank_set <- ranks[0:which(ranks==rank_col)]
          id_hier <- as.list(sub_dt[get(rank_col) == names(props)[1], .SD, .SDcols = rank_set][1])
        }
  
      }
    }
  
    for (r in setdiff(ranks,names(id_hier))) {
      id_hier[r] = NA_character_
    }
  
    result[ranks] <- id_hier
  
    return(result)
  }
                      
  # Generate summary of consensus by BIN
  consensus <- dt[!is.na(get(groups)), do.call(get_consistent_taxon, list(.SD, ranks, threshold, min_ids)), by = groups, .SDcols = ranks]

  # Convert discordant_ids to text if appropriate
  if(discord_format[1] == "text"){
    data.table::set(consensus,
                    j = "discordant_ids",
                    value = sapply(consensus[["discordant_ids"]], function(x) {
                      if (length(x) == 0) return("")
                      sort(x, decreasing = TRUE)
                      pairs <- paste0(names(x), " (", formatC(as.numeric(x), format = "f", digits = 2), ")")
                      paste(pairs, collapse = ", ")
                      })
                    )
  }

  return(consensus)

}

df_summary <- get_bin_consensus(merged, threshold = list(default = 0.75, species = 0.95, subspecies = 0.95), min_ids = 1, groups = "bin")

#####################################################################################################################################
# import new and previous FASTA file 
input.new.fasta <-readDNAStringSet(fasta_file)
input.previous.fasta <- readDNAStringSet("previous_sequences_to_keep.fasta")

# convert previous private record back to PID
previous_private_map <- fread("previous_privatemap_to_keep.tsv", header = T, sep = "\t", fill = TRUE)
lookup <- setNames(previous_private_map$processid, previous_private_map$new_processid)

# revert any anonymized private_record_XXX names in the previous FASTA
previous_names <- names(input.previous.fasta)
prev_ids <- sapply(strsplit(previous_names, "\\|"), "[", 1)
reverted_ids <- ifelse(prev_ids %in% names(lookup),
                       lookup[prev_ids],
                       prev_ids)
bin_parts_prev <- sapply(strsplit(previous_names, "\\|"), "[", 2)
names(input.previous.fasta) <- paste0(reverted_ids, "|", bin_parts_prev)

input.fasta <- c(input.previous.fasta, input.new.fasta)

# make master taxonomy table that includes PIDs
df.final <- data.frame("processid" = sapply(strsplit(names(input.fasta), "\\|"), "[", 1),
                       "bin" = sapply(strsplit(names(input.fasta), "\\|"), "[", 2))
df.final$processid <- ifelse(df.final$processid %in% names(lookup), lookup[df.final$processid], df.final$processid)

# add BOLDPUBLIC to master table
df.final <- merge(df.final, boldlist[,c("processid", "BOLDPUBLIC")], by = "processid", all.x = TRUE)

# add curated taxonomy to master table based on BIN/OTU name
df.final <- merge(df.final, df_summary, by = "bin", all.x = TRUE)
df.final <- df.final[,c("bin",
                        "processid",
                        "BOLDPUBLIC",
                        "concordant_rank",
                        "discordant_rank",
                        "kingdom",        
                        "phylum",
                        "class",
                        "order",
                        "family",
                        "subfamily",
                        "tribe",
                        "genus",
                        "species",
                        "subspecies"
                        )]

# anonymize private PIDs and update them in corresponding FASTA file
# create mapping between old and new PIDs
private_map <- df.final %>%
  filter(BOLDPUBLIC == "PRIVATE") %>%
  mutate(new_processid = paste0("private_record_", seq_len(n()))) %>%
  select(processid, new_processid)

# apply the mapping to update the original data frame
df.final <- df.final %>%
  left_join(private_map, by = "processid") %>%
  mutate(processid = ifelse(BOLDPUBLIC == "PRIVATE", new_processid, processid)) %>%
  select(-new_processid)
df.final <- df.final[,-3] # remove BOLDPUBLIC colummn
name_map <- setNames(private_map$new_processid, private_map$processid)

# Split names
full_names <- names(input.fasta)
sample_ids <- sub("\\|.*", "", full_names)     
bin_parts  <- sub(".*\\|", "", full_names)   

# Replace sample_id if it exists in name_map
new_sample_ids <- ifelse(
  sample_ids %in% names(name_map),
  name_map[sample_ids],
  sample_ids
)

# Recombine new sample ID with original BIN part
new_names <- paste0(new_sample_ids, "|", bin_parts)

# Assign new names
names(input.fasta) <- new_names
############################################################3

# replace NA with blanks in taxonomy table
df.final[is.na(df.final)] <- ""
df.final$discordant_rank[df.final$discordant_rank == ""] <- "None"
df.final$concordant_rank[df.final$concordant_rank == ""] <- "None"

# reduce to one row per BIN
df.final <- unique(df.final[,c(1,5:14,3,4)])

# output FASTA and corresponding taxonomy file and private map
writeXStringSet(input.fasta, fasta_file)
write.table(df.final, gsub("_SEQUENCES.fasta", "_TAXONOMY.tsv", fasta_file), quote = F, row.names = F, sep = "\t")
write.table(private_map, gsub("_SEQUENCES.fasta", "_PRIVATEMAP.tsv", fasta_file), quote = F, row.names = F, sep = "\t")
