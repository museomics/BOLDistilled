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
    ranks=c("kingdom", "phylum", "class", "order", "family", "subfamily", "tribe", "genus", "species", "subspecies"),
    threshold=0.75,
    enforce_scientific=TRUE,
    groups="bin") {
  
  # Regex pattern to recognize non-scientific names in BOLD
  re_int <- paste0(
    "\\.\\Z",                                   # trailing period
    "|\\S{2,}\\.\\S",                           # internal dot between non-whitespace strings
    "|[0-9]",                                   # any digit
    "|\\s[A-Z]",                                # space followed by capital letter
    "|[A-Z]\\Z",                                # capital letter at end
    "|[A-Z]{2}",                                # two consecutive capitals
    "|[a-z][A-Z]",                              # lowercase followed by uppercase
    "|_(?!(hn|sl|ss)\\Z)",                       # underscore, unless used to designate a homonym or sense
    "|%",                                       # percent
    "|\\(",                                     # open parenthesis
    "|\\)",                                     # close parenthesis
    "|,",                                       # comma
    "|\\s(?:aff|agg|cf|complex|group|grp|gr|gp|cmplx|pr|ms|cfr|nr|nsp|near|nomen|hybrid|voucher|form|from|ss|ssl|see|spp?|sample)\\.?(?:\\s|\\Z)" 
  )
  
  # Replace NA in taxonomy columns with empty values (if ignoring non-scientific names, replace those too)
  if(enforce_scientific) {
    df[, (c(ranks)) := lapply(.SD, function(x) fifelse(is.na(x), "", fifelse(grepl(re_int, x, perl = TRUE), "", as.character(x)))), .SDcols = c(ranks)]
  } else {
    df[, (c(ranks)) := lapply(.SD, function(x) fifelse(is.na(x), "", as.character(x))), .SDcols = c(ranks)]
  }
  
  # Convert data table to matrix for faster row access
  mat <- as.matrix(df[, .SD, .SDcols = c(groups, ranks)])
  
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
  
  # Function to compute the consensus taxon for each group
  get_consistent_taxon <- function(sub_dt, ranks, threshold) {
    
    id_hier <- sapply(ranks,function(x) NULL)
    concordant = FALSE
    rank_set <- ranks
    
    result <- list(
      member_count = nrow(sub_dt),
      concordant_rank = NA_character_,
      concordant_id = NA_character_,
      discordant_rank = NA_character_,
      discordant_ids = list()
    )
    
    for (rank_col in rev(ranks)) {  # Step backwards through ranks
      
      vals <- sub_dt[[rank_col]]
      filtered <- vals[!is.na(vals)]
      name_vals <- proportions(table(filtered))
      props <- name_vals[name_vals >= threshold]
      names(name_vals) <- sub("^$","<None>",names(name_vals))
      
      if ((length(props) != 1) & (length(unique(filtered)) > 0)) {
        
        concordant <- FALSE
        
        if(id_hier[rank_col] != "") {
          rank_set <- ranks[0:(which(ranks==rank_col)-1)] 
          id_hier <- id_hier[rank_set]
        }
        
        result$discordant_rank <- rank_col
        result$discordant_ids <- list(setNames(as.vector(name_vals), names(name_vals)))
        
      } else if ((length(props) == 1)) {
        
        if(!concordant) {
          result$concordant_rank <- rank_col
          result$concordant_id <- names(props)[1]
        }
        
        concordant <- TRUE
        
        if ((names(props)[1] != id_hier[rank_col])) {
          rank_set <- ranks[0:which(ranks==rank_col)]
          id_hier <- as.list(sub_dt[get(rank_col) == names(props)[1],..rank_set][1])
        }
        
      }
    }
    
    for (r in setdiff(ranks,names(id_hier))) {
      id_hier[r] = NA_character_
    }
    
    result[ranks] <- id_hier
    
    return(result)
  }
  
  # Generate and return summary of consensus by BIN
  dt[, get_consistent_taxon(.SD, ..ranks, ..threshold), by = eval(groups)]
  
}

df_summary <- get_bin_consensus(merged)


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
