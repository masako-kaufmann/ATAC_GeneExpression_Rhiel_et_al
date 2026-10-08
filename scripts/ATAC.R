



### all on- and off-target sites #################################################
ccr5 = read.csv("ATACseq/bigWig/ATAC_coordinates_CCR5.csv", header = TRUE)
emx1 = read.csv("ATACseq/bigWig/ATAC_coordinates_EMX1.csv", header = TRUE)
fancf = read.csv("ATACseq/bigWig/ATAC_coordinates_FANCF.csv", header = TRUE)

df = rbind(ccr5,emx1,fancf)
dim(df)

### bigwigs #################################################
library(rtracklayer)
bw1 <- import("ATACseq/bigWig/ENCFF233TXT.bigWig")
bw2 <- import("ATACseq/bigWig/ENCFF413OTJ.bigWig")
bw3 <- import("ATACseq/bigWig/ENCFF552WZY.bigWig")
bw4 <- import("ATACseq/bigWig/ENCFF665KOD.bigWig")
bw5 <- import("ATACseq/bigWig/ENCFF863PQY.bigWig")

print(bw1)

### single selection region of interest #################################################
bw_data = bw1 #or other bigwig file

library(GenomicRanges)
region_of_interest <- GRanges("chr16", IRanges(3054701, 3055248))
overlaps <- findOverlaps(region_of_interest, bw_data)
scores_in_region <- bw_data[subjectHits(overlaps)]
print(scores_in_region)

# Plot the scores
plot(scores_in_region$score, type = "l", xlab = "Region Index", ylab = "Score", main = "Score for each DNA region")


### selection region of interest #################################################

head(df)
dim(df)
library(GenomicRanges)

# Function to check if any score exceeds -log10(0.05)

bw_data = bw5

check_ATAC <- function(chromosome, start, end) {
  region_of_interest <- GRanges(chromosome, IRanges(start, end))
  overlaps <- findOverlaps(region_of_interest, bw_data)
  scores_in_region <- bw_data[subjectHits(overlaps)]
  if (any(scores_in_region$score > -log10(0.05))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Add new column df$ATAC
df$ATAC_ENCFF233TXT <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bw1
df$ATAC_ENCFF413OTJ <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bw2
df$ATAC_ENCFF552WZY <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bw3
df$ATAC_ENCFF665KOD <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bw4
df$ATAC_ENCFF863PQY <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bw5

# Print the updated dataframe
print(df)

sum(df$ATAC)
sum(df$ATAC==FALSE)

##############################################################################

#write.csv(df, file = "ATAC_TRUE_FALSE_pvalue05.csv", row.names = FALSE)



##########################################################################################
### ATAC fold change #####################################################################
##########################################################################################

library(rtracklayer)
bwfold1 <- import("ATACseq/bigWig/foldchange/1_ENCFF061KDF.bigWig")
bwfold2 <- import("ATACseq/bigWig/foldchange/2_ENCFF347FGY.bigWig")
bwfold3 <- import("ATACseq/bigWig/foldchange/3_ENCFF396WBB.bigWig")
bwfold4 <- import("ATACseq/bigWig/foldchange/4_ENCFF728DLQ.bigWig")
bwfold5 <- import("ATACseq/bigWig/foldchange/5_ENCFF938UXR.bigWig")



### single selection region of interest #################################################
bw_data = bwfold1

library(GenomicRanges)
region_of_interest <- GRanges("chr16", IRanges(3054701, 3055248))
overlaps <- findOverlaps(region_of_interest, bw_data)
scores_in_region <- bw_data[subjectHits(overlaps)]
print(scores_in_region)

# Plot the scores
plot(scores_in_region$score, type = "l", xlab = "Region Index", ylab = "Score", main = "Score for each DNA region")

### selection region of interest #################################################

head(df)
dim(df)
library(GenomicRanges)

# Function to check if any fold change exceeds 2

bw_data = bw5

check_ATAC <- function(chromosome, start, end) {
  region_of_interest <- GRanges(chromosome, IRanges(start, end))
  overlaps <- findOverlaps(region_of_interest, bw_data)
  scores_in_region <- bw_data[subjectHits(overlaps)]
  if (any(scores_in_region$score > 2)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# Add new column df$ATAC
df$ATAC_ENCFF061KDF <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bwfold1
df$ATAC_ENCFF347FGY <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bwfold2
df$ATAC_ENCFF396WBB <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bwfold3
df$ATAC_ENCFF728DLQ <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bwfold4
df$ATAC_ENCFF938UXR <- mapply(check_ATAC, df$chromosome, df$start, df$end) #bwfold5








