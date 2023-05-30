library(dplyr)
library(tidyverse)
library(stringr)

## Stores the arguments that a user inputs when executing the R script.
args <- commandArgs(trailingOnly=TRUE) 

################################# Read in the files ################################# 
tp_call  <- read.csv(file = paste0(args[1],"/tp_call_INFO.tsv"),
                 sep = "\t",
                 header = FALSE,
                 col.names = c("CHROM_call", "POS_call", "SVTYPE_call", "SVLEN_call", "MatchId","TruScore_call"))

tp_base  <- read.csv(file = paste0(args[1],"/tp_base_INFO.tsv"),
                     sep = "\t",
                     header = FALSE,
                     col.names = c("CHROM_base", "POS_base", "SVTYPE_base", "SVLEN_base", "gnomadId", "MatchId", "TruScore_base"))

################################# Splitting the files ################################# 
## Create the baseID and callID from the MatchId
tp_base[c('chunk', 'base', 'call')] <- stringr::str_split_fixed(tp_base$MatchId, '[.]', 3)
tp_base <- tp_base %>%
  mutate(callId = paste0(chunk, ".", call),
         baseId = paste0(chunk, ".", base)) 

## Split the base tab file into the callID's that appear multiple times, and the ones the appear once
n_occur_base <- data.frame(table(tp_base$callId))
n_occur_base_uniq <- n_occur_base[n_occur_base$Freq == 1,]
n_occur_base_multi <- n_occur_base[n_occur_base$Freq > 1,]

## Create the baseID and callID from the MatchId
tp_call[c('chunk', 'base', 'call')] <- stringr::str_split_fixed(tp_call$MatchId, '[.]', 3)
tp_call <- tp_call %>%
  mutate(callId = paste0(chunk, ".", call),
         baseId = paste0(chunk, ".", base)) 

## Split the call tab file into the baseID's that appear multiple times, and the ones the appear once
n_occur_call <- data.frame(table(tp_call$baseId))
n_occur_call_uniq <- n_occur_call[n_occur_call$Freq == 1,]
n_occur_call_multi <- n_occur_call[n_occur_call$Freq > 1,]

############################################## Unique SVs ############################################## 
## Keep inly the SV's that have unique callId's
tp_base_uniq <- tp_base %>%
  dplyr::filter(callId %in% n_occur_base_uniq$Var1)

## Keep inly the SV's that have unique baseId's
tp_call_uniq <- tp_call %>%
  dplyr::filter(baseId %in% n_occur_call_uniq$Var1)

## Merge the unique SV files by MatchId
tp_uniq <- merge(tp_call_uniq, tp_base_uniq, by="MatchId")

# Select only the columns to output
tp_uniq <- tp_uniq %>%
  mutate(TruScore = TruScore_base) %>%
  dplyr::select(c("CHROM_call", "POS_call", "SVTYPE_call", "SVLEN_call", "CHROM_base", "POS_base", "SVTYPE_base", "SVLEN_base", "gnomadId", "TruScore"))

############################################## Multimatched query SVs ############################################## 
## Pull the SV's from th-call that have been multimatched
tp_call_multi1 <- tp_call %>%
  dplyr::filter(baseId %in% n_occur_call_multi$Var1)

## Pull the SV's from tp-base that have been multimatched
tp_call_multi2 <- tp_base %>%
  dplyr::filter(callId %in% tp_call_multi1$callId)

## Merge the files so that all the repitions of the repeated SV's are included
tp_call_multi <- merge(tp_call_multi1, tp_call_multi2, by="baseId") %>%
  mutate(TruScore = TruScore_call) %>%
  dplyr::select(c("CHROM_call", "POS_call", "SVTYPE_call", "SVLEN_call", "CHROM_base", "POS_base", "SVTYPE_base", "SVLEN_base", "gnomadId", "TruScore"))

############################################## Multimatched base SVs ############################################## 
## Pull the SV's from tp-base that have been multimatched
tp_base_multi1 <- tp_base %>%
  dplyr::filter(callId %in% n_occur_base_multi$Var1)

## Pull the SV's from tp-call that have been multimatched
tp_base_multi2 <- tp_call %>%
  dplyr::filter(baseId %in% tp_base_multi1$baseId)

## Merge the files so that all the repitions of the repeated SV's are included
tp_base_multi <- merge(tp_base_multi1, tp_base_multi2, by="callId", all=TRUE) %>%
  mutate(TruScore = TruScore_base) %>%
  dplyr::select(c("CHROM_call", "POS_call", "SVTYPE_call", "SVLEN_call", "CHROM_base", "POS_base", "SVTYPE_base", "SVLEN_base", "gnomadId", "TruScore"))

############################################## Making the output file ############################################## 
## Combine the uniq and multimatcehed SV files into one dataframe
out <- rbind(tp_uniq, tp_call_multi)
out <- rbind(out, tp_base_multi)

## Write the output file
write.table(out, file=paste0(args[1],"/out.bed"), sep='\t', row.names=FALSE, col.names = FALSE, quote = FALSE)
