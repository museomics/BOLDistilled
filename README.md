# BOLDistilled
BOLDistilled scripts used at CBG


## Process
```
BOLDistill.sh  (the driver: bash, awk, vsearch, seqkit)
 ├─ 1. Filter and clean the BOLD snapshot, cluster records with no BIN into OTUs, keep one copy of each sequence per BIN
 ├─ 2. Distill: within each BIN/OTU, keep only sequences that differ from those already kept
 ├─ 3. Drop OTU representatives that match a BIN representative at ≥95%
 ├─ 4. Rscript BOLDistill.R <SEQUENCES.fasta> <wd>         → agreed taxonomy per BIN, private IDs hidden
 ├─ 5. rmarkdown::render(BOLDistill.Rmd, params=...)       → METADATA.pdf report
 ├─ 6. vsearch --makeudb_usearch / makeblastdb             → VSEARCH and BLAST databases
 └─ 7. python3 BOLDistill_sintax.py <fasta> <TAXONOMY.tsv> → SINTAX reference FASTA
```

## Set up and Run
1. **Install the dependencies**. The script checks for these first and stops if any are missing:
   - Command-line tools: vsearch, seqkit, makeblastdb (from BLAST+), Rscript, python3, gawk (must be the awk on your PATH), unzip, zip, bc.
   - R packages: dplyr, data.table, Biostrings, ggplot2, tidyr, patchwork, DBI, duckdb and rmarkdown, plus a LaTeX install (pdflatex) for the PDF.
   - Python: Biopython.
  
2. **Put the inputs in ~/REFS**. This location is hard-coded (wd="$HOME/REFS", line 52), but the scripts themselves can live anywhere.
3. **Run it** and answer the prompts:
   - Disk: temporary files go to ~/REFS/.tmp, and the snapshot is around 25 GB, so plan for plenty of space.
   - Cores: the core count is nproc − 10, so on a machine with 10 or fewer cores vsearch gets a thread count of zero or less.
```
bash /path/to/BOLDistill.sh
# 1) Public  2) Internal
# library date, e.g. May2026
```

## Input requirements
> Inputs must be in `~/REFS`

- **boldlist.tsv.zip**, or an already unzipped boldlist.tsv: the BOLD snapshot with full taxonomy for public and private records. Columns are found by name, and these must be present:
  - Record details: processid, bin, nuc, marker_code, coi_length, BOLDPUBLIC
  - Flags: filtered, stopcodon, contaminant
  - Taxonomy: kingdom -> subspecies

- **whitelist.txt (public mode only)**: one processid per line for private records allowed into the public library. The first line is skipped as a header.

## Outputs
> Outputs are sent to `~/REFS/BOLDistilled_COI_<date>[_INTERNAL]/`

Internal mode puts everything loose in that directory:
- _SEQUENCES.fasta: the distilled sequences
- _TAXONOMY.tsv: one row per BIN, with kingdom through subspecies and the ranks where members agree or disagree
- _METADATA.pdf: the summary report
- BLAST/, VSEARCH/, SINTAX/: databases ready to search against
- _PRIVATEMAP.tsv: which private_record_N is which real processid
- boldlist_INPUT_<date>.tsv: the raw input, kept as an archive
- _REDUNDANT_OTUs.tsv: an audit of the dropped OTU representatives, flagging matches to a BIN in a different phylum
- problemseqs.tsv / problemtaxa.tsv, if any were produced

