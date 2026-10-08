

library(ggplot2)
library(org.Hs.eg.db)
library(reshape2)
library(dplyr)
library(tibble)
library(biomaRt)
library(tidyr)


df1 = read.table("RNAseq/ENCFF065ZFA.tsv", header = TRUE)
df2 = read.table("RNAseq/ENCFF088XRL.tsv", header = TRUE)
df3 = read.table("RNAseq/ENCFF090JUO.tsv", header = TRUE)
df4 = read.table("RNAseq/ENCFF238QSB.tsv", header = TRUE)
df5 = read.table("RNAseq/ENCFF973CPC.tsv", header = TRUE)

df1 = df1[,c(1,6)]
df2 = df2[,c(1,6)]
df3 = df3[,c(1,6)]
df4 = df4[,c(1,6)]
df5 = df5[,c(1,6)]

merged_df <- merge(df1, df2, by = "gene_id", all = TRUE)
merged_df <- merge(merged_df, df3, by = "gene_id", all = TRUE)
merged_df <- merge(merged_df, df4, by = "gene_id", all = TRUE)
merged_df <- merge(merged_df, df5, by = "gene_id", all = TRUE)

colnames(merged_df) <- c("ENSEMBL", "ENCFF065ZFA", "ENCFF088XRL", "ENCFF090JUO", "ENCFF238QSB", "ENCFF973CPC")
merged_df$ENSEMBL <- gsub("\\.\\d+", "", merged_df$ENSEMBL)
merged_df[900:910, ]
dim(merged_df)

####### FANCF ##############################################
fancf = read.csv("targets/FANCF_expression.csv")
fancf$ENSEMBL <- gsub("\\.\\d+", "", fancf$ENSEMBL)
fancf
fancf_expression <- merge(merged_df, fancf, by = "ENSEMBL")
fancf_expression

####### CCR5 ##############################################
ccr5 = read.csv("targets/CCR5_expression.csv")
ccr5$ENSEMBL <- gsub("\\.\\d+", "", ccr5$ENSEMBL)
ccr5
ccr5_expression <- merge(merged_df, ccr5, by = "ENSEMBL")
ccr5_expression

####### EMX1 ##############################################
emx1 = read.csv("targets/EMX1_expression.csv")
emx1$ENSEMBL <- gsub("\\.\\d+", "", emx1$ENSEMBL)
emx1
emx1_expression <- merge(merged_df, emx1, by = "ENSEMBL")
emx1_expression

####### all ##############################################
ots = rbind(fancf_expression, ccr5_expression, emx1_expression)
ots
dim(ots)
#write.csv(ots, file = "ots_expression.csv", row.names = FALSE)


###########################################
cols_of_interest <- ots[, 2:6]
cols_of_interest
rows_with_two_nonzero <- rowSums(cols_of_interest != 0) >= 2
expressed_at_least_twice <- sum(rows_with_two_nonzero)
expressed_at_least_twice

####### random ##############################################
all_genes <- merged_df[grepl("^ENSG", merged_df$ENSEMBL), ]
all_genes

#set.seed(1)
random_genes <- all_genes[sample(nrow(all_genes), 49), ]
cols_of_interest <- random_genes[, 2:6]
cols_of_interest
rows_with_two_nonzero <- rowSums(cols_of_interest != 0) >= 2
expressed_at_least_twice <- sum(rows_with_two_nonzero)
expressed_at_least_twice











