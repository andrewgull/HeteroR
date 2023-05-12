# run in an in interactive R session
# setwd to where the excel file is located

library(readxl)
library(dplyr)

ecoli_table <- read_excel('500 sepsis E coli overview.xlsx', skip = 12)

new_names <- c("strain", "Illumina", "Nanopore", "Assembly.Status", "DNA.Prep.Status", "PAP.CFX", "PAP.TZP.R.16", "PAP.TZP.R.8", "GM.48h", "GM.24h", 
               "CTX.Stability", "TZP.Stability", "GM.Stability", "CFX.MIC", "CFX.SIR", "TZP.MIC", "TZP.SIR", "GM.MIC", "GM.SIR", "AMP.MIC", "AMP.SIR",
               "TBM.MIC", "TBM.SIR", "TMP.MIC", "TMP.SIR", "NFT.MIC", "NFT.SIR", "MCN.MIC", "MCN.SIR", "CTX.RIOT", "TZP.RIOT", "GM.RIOT", "Antagonism", "Synergy")

ecoli_new_import <- read.csv("ecoli_database_example.csv")

# add new names to xlsx ecoli table, write it to a csv file, import this file to DB4S
ecoli_table <- ecoli_table[,1:34]
names(ecoli_table) <- new_names
ecoli_table <- distinct(ecoli_table)
ecoli_table <- ecoli_table[-c(1,3),]
ecoli_table_cc <- ecoli_table[complete.cases(ecoli_table$strain),]
# this version is fr SQLite data base
write.csv(ecoli_table_cc, "ecoli_database_converted.csv", na = "", row.names = F)

# this one for resistance labeling
ecoli_table_resistance <- select(ecoli_table_cc, -c(Illumina, Nanopore, Assembly.Status, DNA.Prep.Status, 
                                                     PAP.CFX, PAP.TZP.R.16, PAP.TZP.R.8, GM.48h, GM.24h, 
                                                     CTX.Stability, TZP.Stability, GM.Stability, CTX.RIOT, TZP.RIOT, GM.RIOT, Antagonism, Synergy))
write.csv(ecoli_table_resistance, "ecoli_resistance_data.csv", na = "", row.names = F)
