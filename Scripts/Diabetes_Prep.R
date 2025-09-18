#!/usr/bin/env Rscript

# (SOS) - IF I RUN THE SCRIPT THROUGH NEXTFLOW, NOT SURE WHAT THE wd WOULD BE [WE SHALL SEE]
setwd("/home/c0091724/Desktop/PostDoc/Mito_3243/")
library("ggplot2")
library("ggrepel")
library("dplyr")

# Define NOT-IN function
`%!in%` <- Negate(`%in%`)

# GENOTYPE FILE - IDs
geno_fam <- read.table(file="./output_data/Geno_Imputed/merged_imputed_FINAL.fam", header=FALSE)

# COVARIATES [Fam, Age, Sex] & BINARY NMDAS
COVARIATES <- read.table(file="./input_data/FullCohortBinaryPedigree_170821.tab", header=TRUE)
COVARIATES <- COVARIATES[order(COVARIATES[,"famid"]),]
rownames(COVARIATES) <- 1:nrow(COVARIATES)
COVARIATES$ID <- paste(COVARIATES$famid, COVARIATES$pid, sep="_")

# WHICH IDs have GENOTYPES
index <- which(COVARIATES$ID %in% geno_fam$V1)
COVARIATES$Genotype = ifelse(COVARIATES$ID %in% geno_fam$V1, "Yes", "No")
which(geno_fam$V1 %!in% COVARIATES$ID) # 1 ID
geno_fam[258,1]
# [1] "EXLD24562-01_EXLD24562-01"
# (SOS) - THIS IS A MISTAKE, IT IS: 505165057D_EXLD24562-01, SO EDIT COVARIATES$pid to match it
which(COVARIATES$famid == "EXLD24562-01")
# [1] 567
COVARIATES[567,1] <- "EXLD24562-01"
COVARIATES$ID <- paste(COVARIATES$famid, COVARIATES$pid, sep="_")
COVARIATES$Genotype = ifelse(COVARIATES$ID %in% geno_fam$V1, "Yes", "No")
sum(geno_fam$V1 %in% COVARIATES$ID == TRUE)
# 384 GENOTYPES ALL CORRESPOND TO AN INDIVIDUAL IN THE COVARIATES FILE

# 553 INDIVIDUALS WITH NMDAS & HETEROPLASMY MEASUREMENTS:
COVARIATES2 <- subset(COVARIATES, subset=(pid %!in% c(1:12)))
sum(geno_fam$V1 %in% COVARIATES2$ID == TRUE)
# 384
rownames(COVARIATES2) <- 1:nrow(COVARIATES2)
# OUT OF THOSE
sum(COVARIATES2$Genotype=="Yes")
# [1] 384/553 HAVE GENOTYPE DATA
COVARIATES3 <- subset(COVARIATES2, subset=(Genotype == "Yes"))
# 384 INDIVIDUALS

# ASSOCIATION COHORT INFO (BINARY NMDAS)
colnames(COVARIATES3)
#sum(!is.na(COVARIATES3$encephalopathy))
## 253
#nrow(subset(COVARIATES3, subset=(!is.na(encephalopathy) & encephalopathy==1)))
## 49 CASES
#sum(!is.na(COVARIATES3$sle))
# 253
#nrow(subset(COVARIATES3, subset=(!is.na(sle) & sle==1)))
# 56 CASES OF Stroke-Like-Episodes
#sum(!is.na(COVARIATES3$hearing))
## 265
#nrow(subset(COVARIATES3, subset=(!is.na(hearing) & hearing==1)))
## 145 CASES
sum(!is.na(COVARIATES3$diabetes))
# 353 TOTAL
nrow(subset(COVARIATES3, subset=(!is.na(diabetes) & diabetes==1)))
# 205 CASES

# SLE SUBSET
covariates_Diabetes <- COVARIATES3[,c(1:2,47:48,5,16:17,46)]
sum(!is.na(covariates_Diabetes$diabetes))
# 353 from the individuals with genotypes [only EUR] data ALSO have diabetes data
rownames(covariates_Diabetes) <-1:nrow(covariates_Diabetes)

# INDIVIDUALS OF NON-EUR Ancestry (check)
non_eur_ids <- read.table(file="./PCA/Non_EUR_IDs_Full.txt", header=TRUE)
sum(non_eur_ids$Fam_PID %in% covariates_Diabetes$ID == TRUE)
# [1] 0
# CHECK IF DIABETES SUBSET INCLUDES NON-EUR (includes 16 non-EUR)
covariates_Diabetes <- subset(covariates_Diabetes, subset=(pid %!in% non_eur_ids$Individual.ID))
# 353 individuals

# CHECK FOR WHICH OF THESE INDIVIDUALS WE HAVE GENOTYPES
sum(covariates_Diabetes$Genotype=="Yes")
# [1] 353

# DIABETES SUBSET OVERVIEW
sum((covariates_Diabetes$diabetes ==1))
# [1] 205  CASES (148 CONTROLS)
which(is.na(covariates_Diabetes$diabetes))
# integer(0)
which(is.na(covariates_Diabetes$het))
# [1] 26 138 139 167 168 176 183 204 218 222 271 343 --- 12 missing values
which(is.na(covariates_Diabetes$diabetes_age))
# integer(0)
ids_missing_pheno <- covariates_Diabetes[c(which(is.na(covariates_Diabetes$sle_age)),which(is.na(covariates_Diabetes$het))),1]
# 12 individuals missing heteroplasmy
# SO REALLY, ALTHOUGH SAIGE CAN HANDLE NAs IN THE PHENOTYPE FILE, THESE SAMPLES WILL BE IGNORED.
# SO, THE GWAS WILL BE BASED ON 341 SAMPLES:
nrow(subset(covariates_Diabetes, subset=(pid%!in%ids_missing_pheno & diabetes==1)))
# 198 CASES
nrow(subset(covariates_Diabetes, subset=(pid%!in%ids_missing_pheno & diabetes==0)))
# 143 CONTROLS

# MAKE PHENO FILE FOR SAIGE: diabetes, age, het, sex; pID
covariates_diabetes_pheno_file <- covariates_Diabetes[,c(6:7,5,8,1:3)]
colnames(covariates_diabetes_pheno_file)[2] <- "age"
sum(is.na(covariates_diabetes_pheno_file$sex))
#[1] 0
for(i in 1:nrow(covariates_diabetes_pheno_file)){
  if(covariates_diabetes_pheno_file$sex[i] == "M"){covariates_diabetes_pheno_file$sex[i] <- 0}
  else{covariates_diabetes_pheno_file$sex[i] <- 1}}

write.table(covariates_diabetes_pheno_file, file="./output_data/SAIGE/Diabetes/diabetes_pheno.txt", row.names=F, col.names=T, quote=F)
