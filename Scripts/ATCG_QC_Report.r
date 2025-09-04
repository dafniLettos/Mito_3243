#!/usr/bin/env Rscript

# (SOS) - IF I RUN THE SCRIPT THROUGH NEXTFLOW, NOT SURE WHAT THE wd WOULD BE [WE SHALL SEE]
#setwd("../")
#Sys.setenv(HOME="/home/c0091724/Desktop/PostDoc/Mito3243/")
setwd("/home/c0091724/Desktop/PostDoc/Mito_3243/")
library("ggplot2")
library("ggrepel")

#   1. GET AMBIGUOUS SNPs
SNP_dat <- read.delim(file="./input_data/genotypes.bim", header=FALSE)
colnames(SNP_dat) <- c("CHR","ID","DIST","POS","ALLELE1","ALLELE2")
at <- subset(SNP_dat, (ALLELE1 == "A" & ALLELE2 == "T") | (ALLELE1 == "T" & ALLELE2 == "A"))            
gc <- subset(SNP_dat, (ALLELE1 == "C" & ALLELE2 == "G") | (ALLELE1 == "G" & ALLELE2 == "C"))
atgc <- rbind(at, gc)
print(paste0("Number of ambiguous (ATCG) markers: ", nrow(atgc)))
# 40,258
write.table(atgc[2], "./output_data/QC/ATCG_SNPs.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)

#   2. MARKER FILTERS
hwe <- read.table(file="./output_data/QC/QC.hwe", header=TRUE)
maf <- read.table(file="./output_data/QC/QC.frq", header=TRUE)
miss <- read.table(file="./output_data/QC/QC.lmiss", header=TRUE)

QC_markers <- cbind(hwe[c(1:2,4:5,9)], maf[5], miss[5])
colnames(QC_markers)[5] <- "P_hwe"
QC_markers$P_hwe_log10 <- log(QC_markers$P_hwe, base=10)
QC_markers <- QC_markers[,c(1:5,8,6:7)]
#write.table(QC_markers, file = "./output_data/QC/QC_Info_Plink_ALL.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

print(paste0("Total number of markers: ", nrow(QC_markers)))
# 643,848
# 1,931,544
print(paste0("Number of markers that don't pass HWE criterion: ", sum(QC_markers$P_hwe < 0.00000001)))
# 0
print(paste0("Number of markers with MAF < 0.01: ", sum(QC_markers$ALT_FREQS < 0.01)))
# 48,697
print(paste0("Number of markers with Missing Call Rate > 0.3: " , sum(QC_markers$F_MISS > 0.3)))
# 0

#   3. INDIVIDUAL FILTERS
imiss <- read.table(file="./output_data/QC/QC.imiss", header=TRUE)
het <- read.table(file="./output_data/QC/QC.het", header=TRUE)
QC_individuals <- cbind(het, imiss[6])
QC_individuals$HET = (QC_individuals$N.NM. - QC_individuals$O.HOM.)/QC_individuals$N.NM.
# (QC_individuals$OBS_CT - QC_individuals$O.HOM.)/QC_individuals$OBS_CT for PLINK2 output

#length(unique(QC_individuals$X.FID))  # 327
#length(unique(QC_individuals$IID))  # 408

print(paste0("Number of individuals with Missing Call Rate > 0.03: " , sum(QC_individuals$F_MISS > 0.03)))
# 395   Which are equal to 1, potentially "fake" family data individuals
# (SOS) after removing the "fake" individuals, we have 0 with F_MISS > 0.03
#sum(is.na(QC_individuals$HET))
# 395 [that's why no dots in heterozygosity plot]
# 0
#dim(subset(QC_individuals, subset=(F_MISS > 0.03), select=c(IID,HET)))
# There are some proper IDs there, maybe duplicates??? WTF
# ALL GOOD NOW

