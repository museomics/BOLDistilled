#!/bin/bash
# Title: BOLDistill
# Version : 2.0
# Description: This script will take a BOLD snapshot file (formatted to include full taxonomic hierarchy for public and private records) and distill it into the smallest size possible while retaining genetic diversity.
# Author: Sean Prosser (July 2025)

# Usage:
#       BOLDistill.sh        # then answer the interactive prompts:
#                            #   1) library type: Public or Internal
#                            #   2) library date (e.g., May2026)
#       Input is read from, and all output written to, the working directory set below (~/REFS).

# Input files (in the working directory):
#       1) boldlist.tsv.zip
#       2) whitelist.txt   (public mode only — processIDs of private records allowed into public libraries)

# N.B. This script distills the entire BOLD snapshot from scratch on every run (no BINs are carried forward
#      from a previous BOLDistilled library). Previous library folders are left in place in the working
#      directory (no archiving) — move them to permanent storage manually when no longer needed.
##################################################################################################################################################################################################################

# locate this script's directory so the R/Python sub-scripts are found regardless of install location
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

# preflight: confirm every external dependency is present before doing any work
missing=()
for _cmd in vsearch seqkit makeblastdb Rscript python3 gawk unzip zip bc; do
    command -v "$_cmd" >/dev/null 2>&1 || missing+=("$_cmd")
done
if [[ ${#missing[@]} -gt 0 ]]; then
    echo "FATAL: missing required command(s): ${missing[*]}" >&2
    exit 1
fi

# this script depends on GNU awk semantics (IGNORECASE, gensub); mawk silently ignores them
if ! awk --version 2>/dev/null | head -1 | grep -q 'GNU Awk'; then
    echo "FATAL: 'awk' is not GNU awk. Put gawk earlier on PATH (mawk silently mis-filters records)." >&2
    exit 1
fi

for _sub in BOLDistill.R BOLDistill.Rmd BOLDistill_sintax.py; do
    [[ -f "$script_dir/$_sub" ]] || { echo "FATAL: cannot find $_sub in $script_dir" >&2; exit 1; }
done

# set the working directory (input is read from here and all output is written here)
wd="$HOME/REFS"
cd "$wd" || { echo "Cannot enter working directory: $wd" >&2; exit 1; }

# Keep temp-heavy steps OFF /tmp. On this host /tmp is a small RAM-backed tmpfs, and the
# external-merge sort of the ~25 GB snapshot needs more scratch than it holds — it dies with
# "Disk quota exceeded" and takes the machine's memory with it. Point TMPDIR at disk-backed
# scratch inside the working directory instead (covers sort, vsearch and R).
export TMPDIR="$wd/.tmp"
mkdir -p "$TMPDIR"

# prompt for build mode
echo "Select library type:"
echo "  1) Public"
echo "  2) Internal"
read -rp "Enter choice [1-2]: " mode_choice
case "$mode_choice" in
    1) mode="public" ;;
    2) mode="internal" ;;
    *) echo "Invalid choice; please run again and enter 1 or 2." >&2; exit 1 ;;
esac

# prompt for library date (e.g., May2026)
read -rp "Enter library date (e.g., May2026): " library_date
if [[ -z "$library_date" ]]; then
    echo "Library date is required." >&2
    exit 1
fi

# mode-specific variables
if [[ "$mode" == "internal" ]]; then
    current_date="${library_date}_INTERNAL"
else
    current_date="$library_date"
fi

# initialize shared variables
cores=$(nproc)
report_date=$(date +"%d-%b-%Y")
threshold=0.9925
divergence_threshold=$(printf "%.2f%%" "$(echo "scale=4; (1 - $threshold) * 100" | bc)")
boldlistname=$(printf "boldlist_INPUT_%s.tsv" "$library_date")  # e.g., boldlist_INPUT_May2026.tsv

# obtain the boldlist: prefer an already-extracted boldlist.tsv, else unzip the archive.
# (the script consumes the zip, so re-running after a failed build no longer requires re-zipping)
if [[ -s boldlist.tsv ]]; then
    echo "Using existing boldlist.tsv (skipping unzip)."
elif [[ -s boldlist.tsv.zip ]]; then
    unzip -o boldlist.tsv.zip
    rm boldlist.tsv.zip
else
    echo "FATAL: neither boldlist.tsv nor boldlist.tsv.zip found in $wd" >&2
    exit 1
fi

# PUBLIC MODE ONLY: filter out private records, except those in whitelist
if [[ "$mode" == "public" ]]; then
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

    # preserve the raw input and swap in the public-filtered list as boldlist.tsv
    mv boldlist.tsv "$boldlistname"
    mv boldlist_public.tsv boldlist.tsv
fi

# previous library folders are intentionally left in place (no archiving); move them to
# permanent storage manually once they are no longer needed

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
}' boldlist.tsv > boldlist_cleaned.tsv

# extract records that are COI-5P, have a BIN, and a sequence >= 500 bp, no more than 6 Ns, and are not flagged
# "not flagged" tests all three BOLD flag columns: filtered='t' and stopcodon='t' mark bad records,
# and contaminant holds a numeric code when flagged (blank when clean)
awk -F'\t' 'BEGIN {IGNORECASE=1}
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i=="marker_code") marker_col=i;
        if ($i=="bin") bin_col=i;
        if ($i=="nuc") nuc_col=i;
        if ($i=="coi_length") length_col=i;
        if ($i=="filtered") flag_col=i;
        if ($i=="stopcodon") stop_col=i;
        if ($i=="contaminant") cont_col=i;
    }
    print; next;
}
{
    seq_copy = $nuc_col                  # copy sequence
    n_count = gsub(/N/, "", seq_copy)    # count Ns without changing original

    if ($marker_col=="COI-5P" &&
        $bin_col!="" &&
        $nuc_col!="" &&
        $length_col>=500 &&
        $flag_col!="t" &&
        $stop_col!="t" &&
        $cont_col ~ /^[[:space:]]*$/ &&
        n_count<=6) {
        print
    }
}' boldlist_cleaned.tsv > boldlist_with_bins.tsv

# extract sequences (as fasta file) that belong to Bacteria, Fungi, Protista, Nematoda, are COI-5P, are >= 500 bp and <= 1600 bp, do not have BINs, and are not flagged
# (not flagged = filtered!='t' AND stopcodon!='t' AND contaminant blank)
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
        if ($i == "stopcodon") stop_col = i;
        if ($i == "contaminant") cont_col = i;
    }
    next;
}
($marker_col == "COI-5P" && $bin_col == "" && $nuc_col != "" && $length_col >= 500 && $length_col <= 1600 &&
 $flag_col != "t" && $stop_col != "t" && $cont_col ~ /^[[:space:]]*$/ &&
 ($kingdom_col == "Bacteria" || $kingdom_col == "Fungi" || $kingdom_col == "Protista" || $phylum_col == "Nematoda")) {
    if (processid_col != "" && nuc_col != "") {
        print ">" $processid_col "\n" $nuc_col
    }
}' boldlist_cleaned.tsv > bin_ineligible_records.fasta

# cluster BIN-ineligible records into OTUs
vsearch --cluster_fast bin_ineligible_records.fasta \
    --id 0.977 \
    --uc bin_ineligible_records_table.tsv \
    --iddef 3 \
    --threads $((cores - 10))

# reformat OTU output data
# N.B. vsearch --uc output has NO header line: row 1 is already a real 'S' (seed) record.
# Fields are the fixed .uc spec: $1 = record type (S/H/C), $2 = cluster number, $9 = query label.
awk -F'\t' 'BEGIN { OFS="\t" }
{
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
original_records_count=$(tail -n +2 combined_boldlist.tsv | wc -l)
original_records_with_bins_count=$(tail -n +2 boldlist_with_bins.tsv | wc -l)
original_records_with_bin_ineligible_count=$(tail -n +2 boldlist_with_bin-ineligible_records.tsv | wc -l)

# locate columns by name (never by position: BOLD has reordered/added columns between snapshots)
bin_c=$(head -1 combined_boldlist.tsv | tr '\t' '\n' | grep -nxF 'bin' | cut -d: -f1)
nuc_c=$(head -1 combined_boldlist.tsv | tr '\t' '\n' | grep -nxF 'nuc' | cut -d: -f1)
pub_c=$(head -1 combined_boldlist.tsv | tr '\t' '\n' | grep -nxF 'BOLDPUBLIC' | cut -d: -f1)
for _c in bin_c nuc_c pub_c; do
    [[ -n "${!_c}" ]] || { echo "FATAL: could not locate the column for '${_c%_c}' in combined_boldlist.tsv" >&2; exit 1; }
done

# sort by BIN, then BOLDPUBLIC descending so 'PUBLIC' records are preferred as the kept representative,
# then collapse duplicate *sequences* within each BIN (header held out of the sort so it stays on top)
{
    head -n1 combined_boldlist.tsv
    tail -n +2 combined_boldlist.tsv \
      | LC_ALL=C sort -T "$TMPDIR" -t$'\t' -k"${bin_c},${bin_c}" -k"${pub_c},${pub_c}r" -k"${nuc_c},${nuc_c}" \
      | awk -F'\t' -v b="$bin_c" -v n="$nuc_c" '{key = $b FS $n; if (!seen[key]++) print}'
} > combined_boldlist_collapsed.tsv

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
pid_col=$(awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "processid") print i; exit}' "$distill_file")
pub_col=$(awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i == "BOLDPUBLIC") print i; exit}' "$distill_file")
for _c in bin_col nuc_col pid_col pub_col; do
    [[ -n "${!_c}" ]] || { echo "FATAL: could not locate the '${_c%_col}' column in $distill_file" >&2; exit 1; }
done

# extract unique BINs
awk -F'\t' -v bin_col="$bin_col" 'NR > 1 {print $bin_col}' "$distill_file" | sort -u > "$bin_list"

# Pre-split the distill table into one FASTA per BIN in a SINGLE streaming pass.
# Previously every process_bin re-scanned the whole (multi-GB) $distill_file to pull out its own
# BIN — i.e. one full pass per BIN, which dominated the runtime of the entire script.
# combined_boldlist_collapsed.tsv is sorted by bin, then BOLDPUBLIC descending (PUBLIC before
# PRIVATE), then nuc, and boldlist_multiples.tsv preserves that order. Each BIN's records
# therefore arrive contiguously and already PUBLIC-first, so writing them out sequentially
# reproduces exactly the public-then-private FASTA the old two-bucket code assembled by hand.
bin_fasta_dir="bin_fastas"
rm -rf "$bin_fasta_dir"; mkdir -p "$bin_fasta_dir"
awk -F'\t' -v b="$bin_col" -v n="$nuc_col" -v p="$pid_col" -v dir="$bin_fasta_dir" '
    NR == 1 { next }
    {
        if ($b != cur) {
            if (cur != "") close(f)
            cur = $b; safe = cur; gsub(/:/, "-", safe)
            f = dir "/" safe ".fasta"
        }
        print ">"$p"\n"$n > f
    }
    END { if (cur != "") close(f) }
' "$distill_file"

# function to process a single BIN
process_bin() {
  bin="$1"
  master_keep_list="$2"

  # this BIN's FASTA was built by the one-pass pre-split above (PUBLIC sequences first)
  fasta_file="$bin_fasta_dir/$(echo "$bin" | sed 's/:/-/g').fasta"

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
    if [[ $(jobs -r -p | wc -l) -ge $(((cores - 10)/1)) ]]; then
        wait -n
    fi
done
wait

# extract only representative records from boldlist_multiples_filtered.tsv
awk -F'\t' -v pid_col="$pid_col" 'NR==FNR {keep[$1]; next} FNR==1 || ($pid_col in keep)' master_keep.txt boldlist_multiples.tsv > boldlist_multiples_filtered.tsv

###############################################################################
###############################################################################
###############################################################################

# recombine singleton BINs/OTUs and distilled multi-rep BINs/OTUs
{ head -n 1 boldlist_singletons.tsv; tail -n +2 boldlist_singletons.tsv; tail -n +2 boldlist_multiples_filtered.tsv; } > temp_bins.tsv

# generate final file name and rename final TSV file
filename=$(printf "BOLDistilled_COI_%s" "$current_date")
mv temp_bins.tsv "$filename".tsv

###############################################################################
############### De-duplicate OTU (BIN-analog) reps against BINs ###############
###############################################################################
# The BIN-ineligible path (Bacteria/Fungi/Protista/Nematoda without a BIN) is
# meant to add NOVEL diversity as OTUs. An OTU rep that matches an existing BIN
# rep at >= BIN-level identity is not novel — it's a BIN-less record already
# represented by a real BIN, and disproportionately a contaminant/misID (e.g.
# human COI filed under Nematoda: shipped as an OTU it would sit inside the real
# human BIN's k-mer space but carry a Nematoda label, corrupting SINTAX).
# Screen OTU reps against BIN reps here — on the fully distilled set (smallest
# possible) — and drop any OTU rep already covered by a BIN. Everything
# downstream is rebuilt from "$filename".tsv, so this single filter propagates
# to the FASTA, taxonomy, and all reference libraries.
# Threshold is deliberately loose. The OTU fork exists to capture taxa that do not qualify for
# BINs — bacteria (ineligible) and nematodes/fungi/protists (rarely COI-barcoded) — whereas BINs
# are overwhelmingly arthropods and chordates. Genuine novel OTU diversity should therefore be far
# more than 5% divergent from anything holding a BIN, so an OTU rep landing within 95% of a BIN is
# far better explained by contamination than by novel diversity. A stricter 97.7% (BIN-level)
# cutoff was tried and demonstrably leaked vertebrate contamination: 16 full-length Nematoda-
# labelled reps sat at 95.0-97.5% identity to the human BIN (NUMT-like), surviving the filter.
otu_dedup_id=0.95   # main tuning knob

pid_c=$(head -1 "$filename".tsv | tr '\t' '\n' | grep -nxF processid | cut -d: -f1)
bin_dc=$(head -1 "$filename".tsv | tr '\t' '\n' | grep -nxF bin       | cut -d: -f1)
nuc_dc=$(head -1 "$filename".tsv | tr '\t' '\n' | grep -nxF nuc       | cut -d: -f1)
for _c in pid_c bin_dc nuc_dc; do
    [[ -n "${!_c}" ]] || { echo "FATAL: could not locate the '${_c%_c}' column for OTU de-dup in $filename.tsv" >&2; exit 1; }
done

# split distilled reps into BIN reps (search db) and OTU reps (query)
awk -F'\t' -v p="$pid_c" -v b="$bin_dc" -v n="$nuc_dc" '
    NR>1 { if ($b ~ /^BOLD:/)      print ">"$p"\n"$n > "bin_reps.fasta"
           else if ($b ~ /^OTU:/)  print ">"$p"\n"$n > "otu_reps.fasta" }
' "$filename".tsv
touch bin_reps.fasta otu_reps.fasta

# find OTU reps that already belong to an existing BIN
if [[ -s otu_reps.fasta && -s bin_reps.fasta ]]; then
    vsearch --usearch_global otu_reps.fasta \
        --db bin_reps.fasta \
        --id "$otu_dedup_id" \
        --iddef 3 \
        --query_cov 0.80 \
        --maxaccepts 1 \
        --maxrejects 32 \
        --blast6out otu_vs_bin.b6 \
        --threads $((cores - 10))
else
    : > otu_vs_bin.b6
fi
cut -f1 otu_vs_bin.b6 | sort -u > redundant_otu_pids.txt
redundant_otu_count=$(wc -l < redundant_otu_pids.txt)

# audit report: which BIN each dropped OTU matched, flagging cross-phylum hits
# (built from the pre-filter table while it still holds every rep + taxonomy;
#  cross-phylum matches are the misID signal, e.g. a Nematoda OTU -> Mammalia BIN)
awk -F'\t' -v OFS='\t' '
    NR==FNR {
        if (FNR==1) { for (i=1;i<=NF;i++) h[$i]=i; next }
        pid=$(h["processid"]); binv[pid]=$(h["bin"]); phy[pid]=$(h["phylum"]); sp[pid]=$(h["species"]); next
    }
    { q=$1; t=$2; id=$3
      cross=(phy[q]!=phy[t])?"YES":"no"
      print q, id, cross, binv[q], phy[q], sp[q], t, binv[t], phy[t], sp[t] }
' "$filename".tsv otu_vs_bin.b6 \
  | sort -t$'\t' -k3,3r -k2,2rn \
  | cat <(printf 'otu_processid\tpct_id\tcross_phylum\totu_bin\totu_phylum\totu_species\tmatched_bin_processid\tmatched_bin\tbin_phylum\tbin_species\n') - \
  > "$filename"_REDUNDANT_OTUs.tsv

# drop redundant OTU reps from the distilled table
awk -F'\t' -v p="$pid_c" 'NR==FNR{drop[$1]; next} FNR==1 || !($p in drop)' \
    redundant_otu_pids.txt "$filename".tsv > "$filename".dedup.tsv && mv "$filename".dedup.tsv "$filename".tsv

echo ">>>>>> Removed $redundant_otu_count OTU rep(s) already represented by a BIN (>= $otu_dedup_id); see ${filename}_REDUNDANT_OTUs.tsv"
rm -f bin_reps.fasta otu_reps.fasta otu_vs_bin.b6 redundant_otu_pids.txt

###############################################################################

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
Rscript "$script_dir/BOLDistill.R" "$fasta_file_name" "$wd"

# generate summary report
Rscript -e "rmarkdown::render(
  '${script_dir}/BOLDistill.Rmd',
  output_file = sprintf('%s/%s_METADATA.pdf','${wd}','${filename}'),
  params = list(
    wd = '${wd}',
    report_date = '${report_date}',
    filename = '${filename}',
    threshold = '${divergence_threshold}',
    count = '${original_records_count}',
    bin_count = '${original_records_with_bins_count}',
    bin_inel_count = '${original_records_with_bin_ineligible_count}',
    redundant_otu_count = '${redundant_otu_count}'
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
python3 "$script_dir/BOLDistill_sintax.py" $fasta_file_name "$filename"_TAXONOMY.tsv
mv "$filename"_SEQUENCES_sintax.fasta ./SINTAX

# tidy up working directory
if [[ "$mode" == "public" ]]; then
    rm boldlist.tsv            # filtered public temp; raw input already saved as $boldlistname
    rm -f public.tsv allowedprivate.tsv
else
    mv boldlist.tsv "$boldlistname"   # internal: the raw list becomes the archived input
fi
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
rm -rf "$TMPDIR" "$bin_fasta_dir"

# package outputs
if [[ "$mode" == "public" ]]; then
    # build the website bundle: one zip per tool (each carrying a copy of the metadata PDF)
    # plus a loose top-level metadata.pdf
    mkdir -m 777 "./To Website"
    cp "$filename"_METADATA.pdf "./To Website/metadata.pdf"

    # blast.zip -> blast/{db files + metadata}
    mv BLAST blast
    cp "$filename"_METADATA.pdf blast/
    zip -r "./To Website/blast.zip" blast
    rm -r blast

    # vsearch.zip -> vsearch/{udb + metadata}
    mv VSEARCH vsearch
    cp "$filename"_METADATA.pdf vsearch/
    zip -r "./To Website/vsearch.zip" vsearch
    rm -r vsearch

    # sintax.zip -> sintax/{sintax fasta + metadata}
    mv SINTAX sintax
    cp "$filename"_METADATA.pdf sintax/
    zip -r "./To Website/sintax.zip" sintax
    rm -r sintax

    # source.zip -> source/{sequences + taxonomy + metadata}
    mkdir source
    cp "$filename"_METADATA.pdf "$filename"_SEQUENCES.fasta "$filename"_TAXONOMY.tsv source/
    zip -r "./To Website/source.zip" source
    rm -r source

    # remove loose public artifacts (now only inside the website zips)
    rm "$filename"_METADATA.pdf "$filename"_SEQUENCES.fasta "$filename"_TAXONOMY.tsv

    # stage internal-use files (problemseqs/problemtaxa are optional)
    mkdir -m 777 ./INTERNAL_USE
    mv "$boldlistname" *PRIVATEMAP.tsv whitelist.txt ./INTERNAL_USE
    for f in problemseqs.tsv problemtaxa.tsv "${filename}_REDUNDANT_OTUs.tsv"; do
        [[ -e "$f" ]] && mv "$f" ./INTERNAL_USE
    done

    # assemble final library folder
    mkdir -m 777 ./$filename
    mv ./INTERNAL_USE "./To Website" ./$filename
else
    mkdir -m 777 ./$filename
    mv "$boldlistname" *PRIVATEMAP.tsv ./$filename
    for f in problemseqs.tsv problemtaxa.tsv "${filename}_REDUNDANT_OTUs.tsv"; do
        [[ -e "$f" ]] && mv "$f" ./$filename
    done
    mv ./BLAST ./SINTAX ./VSEARCH "$filename"_METADATA.pdf "$filename"_SEQUENCES.fasta "$filename"_TAXONOMY.tsv ./$filename
fi
