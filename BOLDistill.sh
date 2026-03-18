#!/bin/bash
# Title: BOLDistill
# Version : 1.0
# Description: This script will take a BOLD snapshot file (formatted to include full taxonomic hierarchy for public and private records) and distill it into the smallest size possible while retaining genetic diversity.
# Author: Sean Prosser (July 2025)

# Input files:
#       1) boldlist.tsv
#       2) whitelist.tsv (this contains a list of private records (processIDs) that are allowed to be inluded in BOLDistilled libraries)

# N.B. This script assumes that a previous BOLDistilled library exists and only distills BINs that have changed since the previous version. If this is the first use (i.e., no previous BOLDistilled library exists), run code manually, skipping parts that pertain to previous libraries.
##################################################################################################################################################################################################################

# set working directory and initialize variables
cd /home/guelph/REFS/
wd=$(pwd)
cores=$(nproc)
#current_date=$(date +%b%Y)
current_date="Jan2026"
report_date=$(date +"%d-%b-%Y")
threshold=0.9925
divergence_threshold=$(printf "%.2f%%" "$(echo "scale=4; (1 - $threshold) * 100" | bc)")

# unzip boldlist file
unzip -o boldlist.tsv.zip
rm boldlist.tsv.zip

# Filter out private records, except those in whitelist
# Extract column indices
BOLDPUBLIC_COL=$(head -1 boldlist.tsv | tr '\t' '\n' | grep -n -x "BOLDPUBLIC" | cut -d: -f1)
PROCESSID_COL=$(head -1 boldlist.tsv | tr '\t' '\n' | grep -n -x "processid" | cut -d: -f1)

# Ensure the whitelist is clean of carriage returns
sed -i 's/\r$//' whitelist.txt

# Keep rows where BOLDPUBLIC != PRIVATE
awk -F'\t' -v bp_col="$BOLDPUBLIC_COL" '
    NR==1 {print; next}
    $0 && $bp_col != "PRIVATE"
' boldlist.tsv > public.tsv

# Keep rows where BOLDPUBLIC == PRIVATE AND processid is in whitelist
awk -F'\t' -v pid_col="$PROCESSID_COL" -v bp_col="$BOLDPUBLIC_COL" '
    NR==FNR {if (FNR > 1) w[$1]; next}
    NR==1 {print; next}
    $bp_col == "PRIVATE" && $pid_col in w
' whitelist.txt boldlist.tsv > allowedprivate.tsv

# Combine results (remove duplicate header)
head -n1 public.tsv > boldlist_public.tsv
tail -n +2 public.tsv >> boldlist_public.tsv
tail -n +2 allowedprivate.tsv >> boldlist_public.tsv

file_date=$(date +%F)  # e.g., 2025-06-16
boldlistname=$(printf "boldlist_INPUT_%s.tsv" "$file_date")
mv boldlist.tsv "$boldlistname"
mv boldlist_public.tsv boldlist.tsv

# import previous BOLDistilled library and get list of BINs that have not changed
previous_boldlist=$(echo ~/REFS/BOLDistilled*/INTERNAL_USE/boldlist_INPUT*.tsv)
awk -F'\t' 'NR > 1 && $3 != "" { print $3 "\t" $1 }' $previous_boldlist > table1.tsv
awk -F'\t' 'NR > 1 && $3 != "" { print $3 "\t" $1 }' boldlist.tsv > table2.tsv

# Sort process IDs within each bin in table1
awk -F'\t' '{ bins[$1] = bins[$1] $2 "\n" }
    END {
        for (b in bins) {
            split(bins[b], arr, "\n")
            n = asort(arr)
            key = ""
            for (i = 1; i <= n; i++) {
                if (arr[i] != "") key = key arr[i] "|"
            }
            print b "\t" key
        }
    }' table1.tsv | sort > table1_fingerprints.tsv

# Sort process IDs within each bin in table2
awk -F'\t' '{ bins[$1] = bins[$1] $2 "\n" }
    END {
        for (b in bins) {
            split(bins[b], arr, "\n")
            n = asort(arr)
            key = ""
            for (i = 1; i <= n; i++) {
                if (arr[i] != "") key = key arr[i] "|"
            }
            print b "\t" key
        }
    }' table2.tsv | sort > table2_fingerprints.tsv

# Compare the fingerprint files — bins with identical processid sets
comm -12 table1_fingerprints.tsv table2_fingerprints.tsv | cut -f1 > matching_bins.txt

# extract unchanged BINs from previous BOLDistilled library and save for later (will be added back into new BOLDistilled library)
previous_boldistill_fasta=$(echo ~/REFS/BOLDistilled*/PUBLIC/BOLDistilled_COI_*_SEQUENCES.fasta)
awk 'BEGIN { while ((getline < "matching_bins.txt") > 0) bins[$0] = 1 }
     /^>/ {
         split($0, a, "\\|");
         bin = a[2];
         keep = (bin in bins);
     }
     keep' "$previous_boldistill_fasta" > previous_sequences_to_keep.fasta

previous_boldistill_taxonomy=$(echo ~/REFS/BOLDistilled*/PUBLIC/BOLDistilled_COI_*_TAXONOMY.tsv)
awk 'BEGIN { while ((getline < "matching_bins.txt") > 0) bins[$0] = 1 }
     NR == 1 { next }  # Skip the header line
     {
         bin = $1;
         if (bin in bins) {
             print $0;
         }
     }' "$previous_boldistill_taxonomy" > previous_taxonomy_to_keep.tsv

previous_boldistill_privatemap=$(echo ~/REFS/BOLDistilled*/INTERNAL_USE/BOLDistilled_COI_*_PRIVATEMAP.tsv)
awk '
    BEGIN {
        # Load BINs of interest
        while ((getline < "matching_bins.txt") > 0)
            bins[$0] = 1
    }
    # Step 1: Parse FASTA to build new_processid-to-BIN map
    FNR==NR && /^>/ {
        sub(/^>/, "")               # remove leading ">"
        split($0, parts, /\|/)      # split header into PID and BIN
        pid = parts[1]
        bin = parts[2]
        pid2bin[pid] = bin          # map PID to BIN
        next
    }
    # Step 2: Process privatemap TSV (skip header)
    NR!=FNR {
        if (FNR == 1) {
            print                  # print header
            next
        }
        new_pid = $2               # column 2 = new_processid
        bin = pid2bin[new_pid]     # get BIN from FASTA map
        if (bin in bins)
            print                  # output matching rows
    }
' "$previous_boldistill_fasta" "$previous_boldistill_privatemap" > previous_privatemap_to_keep.tsv

# remove unchanged BINs from boldlist
awk 'BEGIN {
    FS = OFS = "\t";
    while ((getline < "matching_bins.txt") > 0)
        bins[$0] = 1
}
NR == 1 {
    print $0; next
}
{
    bin = $3;
    if (!(bin in bins)) {
        print $0;
    }
}' boldlist.tsv > boldlist_filtered.tsv

# move current reference library folder in archive
mv ./BOLDistilled* ./ARCHIVE

# clean up sequences (remove '-' and replace inosine and degenerate bases with N)
gawk -F'\t' 'BEGIN {OFS = FS} 
NR == 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "nuc") nuc_col = i;
    }
    print; 
    next;
} 
{
    if (nuc_col) {
        $nuc_col = gensub(/-/, "", "g", $nuc_col);  # Remove hyphens
        $nuc_col = gensub(/[IRYSWKMBDHV]/, "N", "g", $nuc_col);
        $nuc_col = gensub(/[iryswkmbdhv]/, "N", "g", $nuc_col);
    }
    print;
}' boldlist_filtered.tsv > boldlist_cleaned.tsv

# extract records that are COI-5P, have a BIN, and a sequence >= 500 bp, and are not flagged
awk -F'\t' 'NR == 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "marker_code") marker_col = i;
        if ($i == "bin") bin_col = i;
        if ($i == "nuc") nuc_col = i;
        if ($i == "coi_length") length_col = i;
        if ($i == "filtered") flag_col = i;
    }
    print;
    next;
} 
($marker_col == "COI-5P" && $bin_col != "" && $nuc_col != "" && $length_col >= 500 && $flag_col != "t")' boldlist_cleaned.tsv > boldlist_with_bins.tsv

# extract sequences (as fasta file) that belong to Bacteria, Fungi, Protista, Nematoda, are COI-5P, are >= 500 bp and <= 1600 bp, do not have BINs, and are not flagged
awk -F'\t' 'NR == 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "marker_code") marker_col = i;
        if ($i == "bin") bin_col = i;
        if ($i == "nuc") nuc_col = i;
        if ($i == "coi_length") length_col = i;
        if ($i == "phylum") phylum_col = i;
        if ($i == "kingdom") kingdom_col = i;
        if ($i == "processid") processid_col = i;
        if ($i == "filtered") flag_col = i;
    }
    next;
} 
($marker_col == "COI-5P" && $bin_col == "" && $nuc_col != "" && $length_col >= 500 && $length_col <= 1600 && $flag_col != "t" &&
 ($kingdom_col == "Bacteria" || $kingdom_col == "Fungi" || $kingdom_col == "Protista" || $phylum_col == "Nematoda")) {
    if (processid_col != "" && nuc_col != "") {
        print ">" $processid_col "\n" $nuc_col
    }
}' boldlist_cleaned.tsv > bin_ineligible_records.fasta

# cluster BIN-ineligible records into OTUs
vsearch --cluster_fast bin_ineligible_records.fasta \
    --id 0.988 \
    --uc bin_ineligible_records_table.tsv \
    --iddef 3 \
    --threads $((cores - 2))

# reformat OTU output data
awk -F'\t' 'BEGIN { OFS="\t" }
{
    if (NR == 1) next  # Skip the header
    if ($1 != "C") {
        otu_name = "OTU:" $2
        process_id = $9
        otu_map[process_id] = otu_name
    }
}
END {
    for (id in otu_map) {
        print id, otu_map[id]
    }
}' bin_ineligible_records_table.tsv > otu_map.tsv

# extract bin-ineligible records from cleaned BOLDlist and add OTU names in place of BINs
awk -F'\t' 'BEGIN { OFS="\t"
    # Load lookup table into an array
    while ((getline < "otu_map.tsv") > 0) {
        map[$1] = $2
    }
    close("otu_map.tsv")
}
NR == 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "processid") processid_col = i
        if ($i == "bin") bin_col = i
    }
    print
    next
}
($processid_col in map) {
    if ($bin_col == "") {
        $bin_col = map[$processid_col]
    }
    print
}' boldlist_cleaned.tsv > boldlist_with_bin-ineligible_records.tsv

rm otu_map.tsv bin_ineligible_records_table.tsv

# join BINs and BIN-ineligible tables into one table for distillation
awk 'NR == 1 || FNR > 1' boldlist_with_bins.tsv boldlist_with_bin-ineligible_records.tsv > combined_boldlist.tsv

# count number of records before distillation
unchanged_bin_count=$(awk 'NR==FNR {bins[$1]; next} $3 in bins' matching_bins.txt boldlist.tsv | wc -l)
original_records_count=$(tail -n +2 combined_boldlist.tsv | wc -l)
original_records_count=$((unchanged_bin_count + original_records_count))

original_records_with_bins_count=$(tail -n +2 boldlist_with_bins.tsv | wc -l)
original_records_with_bins_count=$((unchanged_bin_count + original_records_with_bins_count))
original_records_with_bin_ineligible_count=$(tail -n +2 boldlist_with_bin-ineligible_records.tsv | wc -l)

# sort by BIN (column 3) and BOLDPUBLIC (column 28), with 'PUBLIC' entries first, then collapse unique sequences (column 25) within each BIN (column 3)
sort -t$'\t' -k3,3 -k28,28r -k25,25 combined_boldlist.tsv | awk -F'\t' 'NR==1{print; next} {key = $3 FS $25; if (!seen[key]++) print}' > combined_boldlist_collapsed.tsv

# extract singleton BINs/OTUs and process the rest through BOLDistill algorithm
input_file="combined_boldlist_collapsed.tsv"
singletons_file="boldlist_singletons.tsv"
multiples_file="boldlist_multiples.tsv"

# extract header
header=$(head -n1 "$input_file")

# find the index of the "bin" column
bin_col=$(awk -F'\t' 'NR==1 {
    for (i=1; i<=NF; i++) 
        if ($i == "bin") print i;
    exit
}' "$input_file")

# count occurrences of each BIN
awk -F'\t' -v bin_col="$bin_col" '
    NR > 1 {count[$bin_col]++} 
    END {
        for (bin in count) {
            if (count[bin] == 1) print bin > "singleton_bins.tmp"
            else print bin > "multiple_bins.tmp"
        }
    }
' "$input_file"

# isolate singletons
{
    echo "$header"
    awk -F'\t' -v bin_col="$bin_col" 'NR==FNR {singletons[$1]=1; next} (singletons[$bin_col])' singleton_bins.tmp "$input_file"
} > "$singletons_file"

# isolate multiples
{
    echo "$header"
    awk -F'\t' -v bin_col="$bin_col" 'NR==FNR {multiples[$1]=1; next} (multiples[$bin_col])' multiple_bins.tmp "$input_file"
} > "$multiples_file"

rm singleton_bins.tmp multiple_bins.tmp

###############################################################################
########################### BOLDistill Algorithm ##############################
###############################################################################
# for BINs with more than one rep, intelligently select reps for each BIN

bin_list="unique_bins.tmp"
distill_file="boldlist_multiples.tsv"

# find column indices
bin_col=$(awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "bin") print i; exit}' "$distill_file")
nuc_col=$(awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "nuc") print i; exit}' "$distill_file")

# extract unique BINs
awk -F'\t' -v bin_col="$bin_col" 'NR > 1 {print $bin_col}' "$distill_file" | sort -u > "$bin_list"

# function to process a single BIN
process_bin() {
  bin="$1"
  master_keep_list="$2"

  # create a FASTA file containing all sequences for a BIN
  fasta_file="$(echo "$bin" | sed 's/:/-/g').fasta"
  public_seqs_tmp="${bin//:/_}_public_seqs.tmp" 
  private_seqs_tmp="${bin//:/_}_private_seqs.tmp"

  awk -F'\t' -v bin="$bin" -v bin_col="$bin_col" -v nuc_col="$nuc_col" '
    NR == 1 { 
        for (i=1; i<=NF; i++) {
            if ($i == "BOLDPUBLIC") public_col = i;
        }
    }
    ($bin_col == bin) {
        if ($public_col == "PUBLIC") 
            print ">"$1"\n"$nuc_col > "'"$public_seqs_tmp"'";
        else 
            print ">"$1"\n"$nuc_col > "'"$private_seqs_tmp"'";
    }
  ' "$distill_file"

  # concatenate public sequences first, followed by private sequences
  cat "$public_seqs_tmp" "$private_seqs_tmp" > "$fasta_file"
  rm "$public_seqs_tmp" "$private_seqs_tmp"

  # use the provided already-collapsed input
  keep_list="keep_${bin//:/_}.txt"

  # iterate through BIN sequences, retaining only divergent seqs
  while [[ -s "$fasta_file" ]]; do
    # extract the first sequence header and sequence
    focal_header=$(awk '/^>/ {if (NR==1) {print; exit}}' "$fasta_file")
    focal_seq=$(awk '/^>/ {p=0} !/^>/ {if (NR==2) {print; exit}}' "$fasta_file")

    # save the focal header to keep.txt
    echo "${focal_header#>}" >> "$keep_list"

    # stop if no sequences remain
    [[ -z "$focal_header" ]] && rm -f "$fasta_file" && break

    # create a new focal FASTA file
    focal_fasta="focal_${bin//:/_}.fasta"
    echo "$focal_header" > "$focal_fasta"
    echo "$focal_seq" >> "$focal_fasta"

    # run VSEARCH to find similar sequences
    vsearch --usearch_global "$focal_fasta" \
      --db "$fasta_file" \
      --id $threshold \
      --blast6out "hits_${bin//:/_}.txt" \
      --self \
      --maxaccepts 0 \
      --maxrejects 0 \
      --query_cov 0.0 \
      --target_cov 0.0 \
      --threads 1

    wait

    # extract sequence headers of similar sequences (to be removed)
    awk '{print $2}' "hits_${bin//:/_}.txt" | sort -u > "to_remove_${bin//:/_}.txt"

    # add the focal header to the remove list
    echo "${focal_header#>}" >> "to_remove_${bin//:/_}.txt"

    # filter out similar sequences
    seqkit grep -v -f "to_remove_${bin//:/_}.txt" "$fasta_file" -w 0 > "temp_${bin//:/_}.fasta" && mv "temp_${bin//:/_}.fasta" "$fasta_file"
    rm "$focal_fasta" "hits_${bin//:/_}.txt" "to_remove_${bin//:/_}.txt"
  done

  # append to the master keep.txt
  cat "$keep_list" >> "$master_keep_list"
  rm -f "$fasta_file" "$keep_list"
}

export -f process_bin

# run the above function in parallel (one core per BIN)
master_keep_list="master_keep.txt"
> "$master_keep_list"
total_bins=$(wc -l < "$bin_list")
i=0  # Initialize counter

for bin in $(cat "$bin_list"); do
    ((i++))

    # calculate percentage completed
    percent_complete=$(( (i * 100) / total_bins ))

    # print progress
    echo ">>>>>>>>>>>>>>>>>>>>>>>> Progress: $percent_complete% ($i of $total_bins bins)"

    # process the BIN
    process_bin "$bin" "$master_keep_list" &

    # limit parallel jobs
    if [[ $(jobs -r -p | wc -l) -ge $(((cores - 2)/1)) ]]; then
        wait -n
    fi
done
wait

# extract only representative records from boldlist_multiples_filtered.tsv
awk 'NR==FNR {keep[$1]; next} FNR==1 || ($1 in keep)' master_keep.txt boldlist_multiples.tsv > boldlist_multiples_filtered.tsv

###############################################################################
###############################################################################
###############################################################################

# recombine singleton BINs/OTUs and distilled multi-rep BINs/OTUs
{ head -n 1 boldlist_singletons.tsv; tail -n +2 boldlist_singletons.tsv; tail -n +2 boldlist_multiples_filtered.tsv; } > temp_bins.tsv

# generate final file name and rename final TSV file
filename=$(printf "BOLDistilled_COI_%s" "$current_date")
mv temp_bins.tsv "$filename".tsv

# extract sequences and headers and convert to FASTA
awk -F'\t' 'BEGIN {OFS = FS} 
NR == 1 {
  for (i=1; i<=NF; i++) colname[$i] = i;  # Create a mapping of column names to their positions
  next
}
{
  processid = $colname["processid"]
  bin = $colname["bin"]
  nuc = $colname["nuc"]

  fasta_header = ">" processid "|" bin
  print fasta_header
  print nuc
}' "$filename".tsv > "$filename"_SEQUENCES.fasta

fasta_file_name="$filename"_SEQUENCES.fasta

# determine best taxonomic hierarchy for each BIN or other kingdoms without BINs,and output table for later use
Rscript /home/guelph/SCRIPTS/SUBSCRIPTS/BOLDistill.R $fasta_file_name $wd

# generate summary report
Rscript -e "rmarkdown::render(
  '/home/guelph/SCRIPTS/SUBSCRIPTS/BOLDistill.Rmd',
  output_file = sprintf('%s/%s_METADATA.pdf','${wd}','${filename}'),
  params = list(
    wd = '${wd}',
    report_date = '${report_date}',
    filename = '${filename}',
    threshold = '${divergence_threshold}',
    count = '${original_records_count}',
    bin_count = '${original_records_with_bins_count}',
    bin_inel_count = '${original_records_with_bin_ineligible_count}'
  ),
  clean = TRUE)"

# convert sequence file to single-line
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $fasta_file_name > "$fasta_file_name"_2
rm $fasta_file_name
mv "$fasta_file_name"_2 $fasta_file_name

# create VSEARCH reference library
mkdir -m 777 VSEARCH
vsearch --makeudb_usearch $fasta_file_name --output ./VSEARCH/"${fasta_file_name%.fasta}_vsearch"

# create BLAST reference library
mkdir -m 777 BLAST
makeblastdb -dbtype 'nucl' -in $fasta_file_name -input_type 'fasta' -out ./BLAST/$fasta_file_name

# create SINTAX reference library
mkdir -m 777 ./SINTAX
python3 /home/guelph/SCRIPTS/SUBSCRIPTS/BOLDistill_sintax.py $fasta_file_name "$filename"_TAXONOMY.tsv
mv "$filename"_SEQUENCES_sintax.fasta ./SINTAX

# tidy up directory
rm boldlist.tsv
rm boldlist_filtered.tsv
rm unique_bins.tmp
rm boldlist_singletons.tsv
rm boldlist_multiples.tsv
rm boldlist_multiples_filtered.tsv
rm bin_ineligible_records.fasta
rm boldlist_with_bin-ineligible_records.tsv
rm combined_boldlist.tsv
rm combined_boldlist_collapsed.tsv
rm boldlist_cleaned.tsv
rm boldlist_with_bins.tsv
rm master_keep.txt
rm "$filename".tsv
rm matching_bins.txt
rm previous_sequences_to_keep.fasta
rm previous_taxonomy_to_keep.tsv
rm previous_privatemap_to_keep.tsv
rm table1.tsv
rm table1_fingerprints.tsv
rm table2.tsv
rm table2_fingerprints.tsv
rm public.tsv
rm allowedprivate.tsv

mkdir -m 777 ./INTERNAL_USE
mv $boldlistname problemseqs.tsv problemtaxa.tsv *PRIVATEMAP.tsv whitelist.txt ./INTERNAL_USE

mkdir -m 777 ./PUBLIC
mv ./BLAST ./SINTAX ./VSEARCH "$filename"_METADATA.pdf "$filename"_SEQUENCES.fasta "$filename"_TAXONOMY.tsv ./PUBLIC

mkdir -m 777 ./$filename
mv ./INTERNAL_USE ./PUBLIC ./$filename