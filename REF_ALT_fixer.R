##### THIS MAKes the bare bones vcf
library("Biostrings")
library(dplyr)
library(tidyverse)
library(tibble)

## Stores the arguments that a user inputs when executing the R script.
args <- commandArgs(trailingOnly=TRUE) 

# # For gnomad
# bed  <- read.csv(file = paste0(args[1], "lifted_filtered_sorted.bed"),
#                  sep = "\t",
#                  header = FALSE,
#                  col.names = c("CHROM", "POS", "END", "ID", "SVTYPE"))

# For real files
## Reading in the bed file
bed  <- read.csv(file = paste0(args[1], "body2.tab"),
                 sep = "\t",
                 header = FALSE,
                 col.names = c("CHROM", "POS", "END", "ID", "SVTYPE"))

## Reading in the FASTA sequences
fastaFile <- readDNAStringSet(paste0(args[1], "temp.tab"))
seq_name = names(fastaFile)
sequence = paste(fastaFile)
df <- data.frame(seq_name, sequence)

## Combining the BED file and the FASTA sequences
out <- cbind(bed, df)

## Making the REF and ALT to be how Truvari wants it and creating the INFO field.
out <- bed %>%
  mutate(REF = case_when(SVTYPE == "INS" ~ substr(sequence, 1, 1),
                         SVTYPE == "DEL" ~ sequence),
         ALT = case_when(SVTYPE == "INS" ~ paste0(substr(sequence, 1, 1), strrep("N", nchar(sequence) - 1)),
                         SVTYPE == "DEL" ~ substr(sequence, 1, 1)),
         QUAL = ".",
         FILTER = "PASS",
         SVLEN = case_when(SVTYPE == "INS" ~ END-POS,
                           TRUE ~ -(END-POS)),
         END = case_when(SVTYPE == "INS" ~ POS,
                         TRUE ~ END),
         INFO = paste0("END=", END, ";SVLEN=", SVLEN, ";SVTYPE=", SVTYPE),
         FORMAT = "GT",
         MYSAMP = "0|1") %>%
  dplyr::select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, MYSAMP)

## Writing the VCF body file
write.table(out, file=paste0(args[1], "new_bod.tab"), sep='\t', row.names=FALSE, col.names = FALSE, quote = FALSE)