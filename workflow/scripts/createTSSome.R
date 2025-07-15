#!/usr/bin/env Rscript
### PURPOSE OF THIS SCRIPT
## Create a TSS FASTA and GTF file from the TSScall output.
## 
## Input:
##  1) TSScall output bedgraph file
##  2) Genomic FASTA file
##
## Output:
##  1) TSS FASTA file
##  2) TSS GTF file


# Parse command line args ------------------------------------------------------

library(optparse)

args <- commandArgs(trailingOnly = TRUE)

option_list <- list(
  make_option(
    c("--fasta", type = "character"),
    help = "Path to genome FASTA file"
  ),
  make_option(
    c("--bed", type = "character"),
    help = "Path to TSS bed file"
  ),
  make_option(
    c("--details", type = "character"),
    help = "Path to TSS bed details file"
  ),
  make_option(
    c("--output_fasta", type = "character"),
    help = "Path to output FASTA file"
  ),
  make_option(
    c("--output_gtf", type = "character"),
    help = "Path to output GTF file"
  ),
  make_option(
    c("--keep_uTSS", type = "logical"),
    default = FALSE,
    help = "Keep unobserved TSSs (name: uTSS<ID>)?"
  ),
  make_option(
    c("--merge_clusters", type = "logical"),
    default = TRUE,
    help = "Merge clusters of TSSs into a single TSS?"
  ),
  make_option(
    c("--upstream", type = "numeric"),
    default = 100,
    help = "Number of nucleotides upstream of TSS site call to include
    in FASTA"
  ),
  make_option(
    c("--downstream", type = "numeric"),
    default = 100,
    help = "Number of nucleotides downstream of TSS site call to include
    in FASTA"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser) # Load options from command line.


# Load dependencies ------------------------------------------------------------

library(rtracklayer)
library(Biostrings)
library(dplyr)
library(GenomicFeatures)
library(data.table)
library(GenomicRanges)
library(BSgenome)

# The script -------------------------------------------------------------------

##### Parameters #####

fa_path <- opt$fasta


##### STEP 1: Get TSS sequences #####

### Load input files
TSSbed <- fread(
  opt$bed,
  skip = 1,
  header = FALSE,
  col.names = c("seqnames", "start", "end", "TSSid", "filler", "strand")
)

TSSdeets <- fread(
  opt$details
)


fasta <- readDNAStringSet(
  fa_path
)

### Filter

if(!opt$keep_uTSS){
  
  TSSbed <- TSSbed %>%
    dplyr::filter(
      !grepl("^uTSS", TSSid)
    )
  
  TSSdeets <- TSSdeets %>%
    dplyr::filter(
      !grepl("^uTSS", `TSS ID`)
    )
  
}


### Create TSS GenomicRanges object

# BED files are 0-indexed, so there are some extra +1's here
TSSgr <- GenomicRanges::GRanges(
  seqnames = TSSbed$seqnames,
  ranges = IRanges(
    start = pmax((TSSbed$start + 1) - opt$upstream, 1),
    end = TSSbed$end + opt$downstream + 1
  ),
  strand = Rle(TSSbed$strand),
  TSSid = TSSbed$TSSid
)


# Merge overlapping TSS windows to avoid redundant sequences in the final
# FASTA file
TSSgr_merged <- reduce(
  TSSgr,
  with.revmap = TRUE,
  ignore.strand = FALSE
)

rev <- mcols(TSSgr_merged)$revmap

# Bit inefficient
# Don't need to do this for list elements of length 1 (i.e., 1 TSS went into
# the merged TSS)
mcols(TSSgr_merged)$TSSid <- sapply(
  rev,
  function(idx){
    
    # Don't do anything crazy for length 1 elements
    if(length(idx) == 1){
      
      return(mcols(TSSgr)$TSSid[idx])
      
    }
    
    possible_TSS <- mcols(TSSgr)$TSSid[idx]
    
    # Check to see which overlap with an actual gene
    # NOTE: may be redundant if keep_uTSS == FALSE, as I think
    # observed TSS are those that can be mapped to a gene.
    TSS_table <- TSSdeets %>%
      dplyr::filter(`TSS ID` %in% possible_TSS) %>%
      dplyr::select(`TSS ID`, `Gene ID`) %>%
      na.omit()
    
    if(nrow(TSS_table) == 0){
      
      return(possible_TSS[1])
      
    }else{
      
      return(TSS_table$`TSS ID`[1])
      
    }
    
  }
)


### Make sure there are no illegal ranges (less than 1 and greater than seqlength)

seqs_in_TSS <- seqnames(TSSgr)
seqlens <- seqlengths(fasta)

seqlens_filter <- seqlens[names(seqlens) %in% seqs_in_TSS]

start_names <- names(seqlens_filter)
starts <- rep(1, times = length(start_names))
names(starts) <- start_names

TSSgr_merged <- GenomicRanges::restrict(TSSgr_merged,
                                        start = 1,
                                        end = seqlens_filter)

### Get sequences

TSS_seqs <- BSgenome::getSeq(
  fasta,
  TSSgr_merged
)

##### STEP 2: Create and write final output #####

names(TSS_seqs) <- TSSgr_merged$TSSid

# FASTA file
writeXStringSet(
  TSS_seqs,
  opt$output_fasta,
  format = "fasta"
)


# GTF file
rtracklayer::export(
  TSSgr_merged,
  opt$output_gtf
)