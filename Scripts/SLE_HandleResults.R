setwd("/home/c0091724/Desktop/PostDoc/Mito_3243/")
library(data.table)
library(qvalue)
library(qqman)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(gridExtra)
library(dplyr)

# Define "NOT in" function
`%!in%` <- Negate(`%in%`)

# Read-in Group Files
gf_chr5 <- read.table("./output_data/Geno_Imputed/Group_File_chr5", header=F, col.names = paste0("V",seq_len(max(count.fields("./output_data/Geno_Imputed/Group_File_chr5", sep = '\t')))), fill = TRUE)
gf_chr6 <- read.table("./output_data/Geno_Imputed/Group_File_chr6", header=F, col.names = paste0("V",seq_len(max(count.fields("./output_data/Geno_Imputed/Group_File_chr6", sep = '\t')))), fill = TRUE)

# Read-in SAIGE-GENE Results
res_files_ds <- paste0(rep("SAIGE_SBA_DS_chr",22),c(1:22))
res_GB_DS <- list()
for(file in res_files_ds){
  print(paste0("./output_data/SAIGE/Encephalopathy/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/SLE/Results/",file), header=T)
  res_GB_DS[[file]] <- data}

# Read-in SAIGE Results (GWAS)
res_files_res_ALL001_ds <- paste0(rep("SAIGE_SVA_maf001_chr",22),c(1:22),rep(".txt",22))
res_SVA001_DS <- list()
for(file in res_files_res_ALL001_ds){
  print(paste0("./output_data/SAIGE/SLE/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/SLE/Results/",file), header=T)
  res_SVA001_DS[[file]] <- data}

res_files_res_ALL005_ds <- paste0(rep("SAIGE_SVA_maf005_chr",22),c(1:22),rep(".txt",22))
res_SVA005_DS <- list()
for(file in res_files_res_ALL005_ds){
  print(paste0("./output_data/SAIGE/SLE/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/SLE/Results/",file), header=T)
  res_SVA005_DS[[file]] <- data}
#----------------------------------
#         QQ Plot Function
#----------------------------------
# Create a quantile-quantile plot with ggplot2.
# Assumptions:
#   - Expected P values are uniformly distributed.
#   - Confidence intervals assume independence between tests. We expect deviations past the confidence intervals
#     if the tests are not independent. For example, in a genome-wide association study, the genotype at any
#     position is correlated to nearby positions. Tests of nearby genotypes will result in similar test statistics.
# ps: Vector of p-values.
# ci: Size of the confidence interval, 95% by default.
qqplot_test <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps, na.last = TRUE)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1)))
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(mapping = aes(x = expected, ymin = clower, ymax = cupper), alpha = 0.1, fill = "orchid3") +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    #labs(x = "Expected -log<sub>10</sub>(p)", y = "Observed -log<sub>10</sub>(p)") +
    #ylab("Observed -log<sub>10</sub>(p)") +
    #xlab("Expected -log<sub>10</sub>(p)") +
    xlab(log10Pe) + ylab(log10Po) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.7, col = "orchid4", lwd=1)
}
#----------------------------------------------------------
#----------------------------------------------------------
#                     GENOME-WIDE TEST
#----------------------------------------------------------
#----------------------------------------------------------
# Merge results per chromosome into single dataframe [MAF 0.01]
res_ALL001_DS_ALL <- bind_rows(res_SVA001_DS)
# 7,586,961 SNPs
write.table(res_ALL001_DS_ALL, file="./output_data/SAIGE/SLE/Results/SAIGE_SVA_DS_ALL_MAF001.txt", sep="\t", col.names=T, row.names=F, quote=F)
hist(res_ALL001_DS_ALL$p.value)

# Merge results per chromosome into single dataframe
res_ALL005_DS_ALL <- bind_rows(res_SVA005_DS)
# 5,464,655 SNPs
write.table(res_ALL005_DS_ALL, file="./output_data/SAIGE/SLE/Results/SAIGE_SVA_DS_ALL_MAF005.txt", sep="\t", col.names=T, row.names=F, quote=F)
hist(res_ALL005_DS_ALL$p.value)
#--------------------------
#         CHR6  [DS]
#--------------------------
res_ALL001_DS_chr6 <- res_SVA001_DS$SAIGE_SVA_maf001_chr6.txt
res_ALL001_DS_chr6_reg <- subset(res_ALL001_DS_chr6, subset=((POS >= 143863601)&(POS <=173863601)))
hist(res_ALL001_DS_chr6_reg$p.value)
sum(res_ALL001_DS_chr6_reg$p.value < 5e-08)
#[1] 0
res_ALL001_DS_chr6_reg$p.value_log <- -log10(as.numeric(res_ALL001_DS_chr6_reg$p.value))
res_ALL001_DS_chr6_reg$SNP_ID <- paste(res_ALL001_DS_chr6_reg$MarkerID,res_ALL001_DS_chr6_reg$Allele1,res_ALL001_DS_chr6_reg$Allele2,sep=":")

res_ALL005_DS_chr6 <- res_SVA005_DS$SAIGE_SVA_maf005_chr6.txt
res_ALL005_DS_chr6_reg <- subset(res_ALL005_DS_chr6, subset=((POS >= 143863601)&(POS <=173863601)))
hist(res_ALL005_DS_chr6_reg$p.value)
sum(res_ALL005_DS_chr6_reg$p.value < 5e-08)
#[1] 0
res_ALL005_DS_chr6_reg$p.value_log <- -log10(as.numeric(res_ALL005_DS_chr6_reg$p.value))
res_ALL005_DS_chr6_reg$SNP_ID <- paste(res_ALL005_DS_chr6_reg$MarkerID,res_ALL005_DS_chr6_reg$Allele1,res_ALL005_DS_chr6_reg$Allele2,sep=":")
#--------------------------
#         CHR5  [DS]
#--------------------------
res_ALL001_DS_chr5 <- res_SVA001_DS$SAIGE_SVA_maf001_chr5.txt
res_ALL001_DS_chr5_reg <- subset(res_ALL001_DS_chr5, subset=((POS >= 140749194)&(POS <=188279842)))
hist(res_ALL001_DS_chr5_reg$p.value)
sum(res_ALL001_DS_chr5_reg$p.value < 5e-08)
#[1] 0
res_ALL001_DS_chr5_reg$p.value_log <- -log10(as.numeric(res_ALL001_DS_chr5_reg$p.value))
res_ALL001_DS_chr5_reg$SNP_ID <- paste(res_ALL001_DS_chr5_reg$MarkerID,res_ALL001_DS_chr5_reg$Allele1,res_ALL001_DS_chr5_reg$Allele2,sep=":")

res_ALL005_DS_chr5 <- res_SVA005_DS$SAIGE_SVA_maf005_chr5.txt
res_ALL005_DS_chr5_reg <- subset(res_ALL005_DS_chr5, subset=((POS >= 140749194)&(POS <=188279842)))
hist(res_ALL005_DS_chr5_reg$p.value)
sum(res_ALL005_DS_chr5_reg$p.value < 5e-08)
#[1] 0
res_ALL005_DS_chr5_reg$p.value_log <- -log10(as.numeric(res_ALL005_DS_chr5_reg$p.value))
res_ALL005_DS_chr5_reg$SNP_ID <- paste(res_ALL005_DS_chr5_reg$MarkerID,res_ALL005_DS_chr5_reg$Allele1,res_ALL005_DS_chr5_reg$Allele2,sep=":")
#---------------------------------------------------------------
#                       Manhattan Plots
#---------------------------------------------------------------
#---------------------------------------
#   CHR6 region of suggestive linkage
#---------------------------------------
## Suggestive significance level for the proportion of the genome tested:
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr6_reg)  # 4.37e-06  -log10(sig) --- 5.36
ggplot(res_ALL001_DS_chr6_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  #geom_label_repel(data = subset(res_ALL001_DS_chr6_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 20, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("midnightblue", "midnightblue"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(143000000,174000000,by=5000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
  theme_minimal() +
  labs(title="SAIGE SVA: Encephalopathy Associations",
       x ="",
       y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))

sig <- (nrow(res_ALL005_DS_ALL) * 5e-08) / nrow(res_ALL005_DS_chr6_reg)  # 4.23e-06  -log10(sig) --- 5.37
ggplot(res_ALL005_DS_chr6_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL005_DS_chr6_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 20, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("midnightblue", "midnightblue"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(143000000,174000000,by=5000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
  theme_minimal() +
  labs(title="SAIGE SVA: Encephalopathy Associations",
       x ="",
       y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))

#---------------------------------------
#   CHR5 region of suggestive linkage
#---------------------------------------
## Suggestive significance level for the proportion of the genome tested:
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr5_reg)  # 3.34e-06  -log10(sig) --- 5.48
ggplot(res_ALL001_DS_chr5_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL001_DS_chr5_reg, subset=(-log10(p.value) >= 4.1)), size = 3, max.overlaps = 20, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("midnightblue", "midnightblue"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(140000000,189000000,by=5000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
  theme_minimal() +
  labs(title="SAIGE SVA: Encephalopathy Associations",
       x ="",
       y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))

sig <- (nrow(res_ALL005_DS_ALL) * 5e-08) / nrow(res_ALL005_DS_chr5_reg)  # 3.34e-06  -log10(sig) --- 5.48
ggplot(res_ALL005_DS_chr5_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL005_DS_chr5_reg, subset=(-log10(p.value) >= 4.1)), size = 20, max.overlaps = 2, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("midnightblue", "midnightblue"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(140000000,189000000,by=5000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
  theme_minimal() +
  labs(title="SAIGE SVA: Encephalopathy Associations",
       x ="",
       y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))

#---------------------------------------------------------------------
#               Genomic inflation factor (lambda)
#---------------------------------------------------------------------
qchisq(1-median(res_ALL001_DS_ALL$p.value),1)/qchisq(0.5,1)
# [1] 1.018247
qchisq(1-median(res_ALL005_DS_ALL$p.value),1)/qchisq(0.5,1)
# [1] 1.010453
#-----------------------------------------------------------------
#                           QQ Plots
#-----------------------------------------------------------------
qqplot_test(res_ALL001_DS_ALL$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) Genome-Wide")
qqplot_test(res_ALL005_DS_ALL$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) Genome-Wide")
#----------------------------------------------------------
#----------------------------------------------------------
#                     GENE-BURDEN TEST
#----------------------------------------------------------
#----------------------------------------------------------
# How many genes were tested?
length(unique(res_GB_DS$SAIGE_SBA_DS_chr6$Region))
# 909 GENES on Chromosome 6
length(unique(res_GB_DS$SAIGE_SBA_DS_chr5$Region))
# 748 GENES on Chromosome 5

# Initialise SNPs dataframe that will contain variant information
SNPs_chr6 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(SNPs_chr6) <- c("Variant","Annotation","Gene")
SNPs_chr5 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(SNPs_chr5) <- c("Variant","Annotation","Gene")

# FOR LOOP to get some extra columns of info for variants per gene in the group file & create the SNPs dataframe
#--------------------------------------
#            Chromosome 6
#--------------------------------------
# These column in the group file 
gf_chr6$N_Var <- ""
gf_chr6$Variants <- ""

for(i in 1:nrow(gf_chr6)){
  if(gf_chr6[i,2] == "var"){
    empty <- sum(gf_chr6[i,] == "") - 2
    nVar <- ncol(gf_chr6) - empty - 4
    gene <- gf_chr6[i,1]
    
    gf_chr6$N_Var[i] <- nVar
    gf_chr6$Variants[i] <- paste(gf_chr6[i,3:(nVar+2)], collapse="|")
    SNPs_tmp <- data.frame(
            Variant = c(t(gf_chr6[i,3:(nVar+2)])),
            Annotation = c(t(gf_chr6[i+1,3:(nVar+2)])),
            Gene = rep(gene, times=nVar))
    SNPs_chr6 <- rbind(SNPs_chr6, SNPs_tmp)}
  else{
    gf_chr6$N_Var[i] <- ""
    gf_chr6$Variants[i] <- ""}
}
# Add POSITION column to SNPs dataframe
SNPs_chr6$SNP_POS <- sapply(strsplit(SNPs_chr6$Variant,":"), `[`, 2)
SNPs_chr6$Gene <- sapply(strsplit(SNPs_chr6$Gene,","), `[`, 1)
#--------------------------------------
#            Chromosome 5
#--------------------------------------
# These column in the group file 
gf_chr5$N_Var <- ""
gf_chr5$Variants <- ""

for(i in 1:nrow(gf_chr5)){
  if(gf_chr5[i,2] == "var"){
    empty <- sum(gf_chr5[i,] == "") - 2
    nVar <- ncol(gf_chr5) - empty - 4
    gene <- gf_chr5[i,1]
    
    gf_chr5$N_Var[i] <- nVar
    gf_chr5$Variants[i] <- paste(gf_chr5[i,3:(nVar+2)], collapse="|")
    SNPs_tmp <- data.frame(
      Variant = c(t(gf_chr5[i,3:(nVar+2)])),
      Annotation = c(t(gf_chr5[i+1,3:(nVar+2)])),
      Gene = rep(gene, times=nVar))
    SNPs_chr5 <- rbind(SNPs_chr5, SNPs_tmp)}
  else{
    gf_chr5$N_Var[i] <- ""
    gf_chr5$Variants[i] <- ""}
}
# Add POSITION column to SNPs dataframe
SNPs_chr5$SNP_POS <- sapply(strsplit(SNPs_chr5$Variant,":"), `[`, 2)
SNPs_chr5$Gene <- sapply(strsplit(SNPs_chr5$Gene,","), `[`, 1)
#------------------------------------------------------
#------------------------------------------------------
## Add GENE START for POSITION HERE
# Will be using the file from ANNOVAR
anno_genes <- fread("./annovar/humandb/hg19_refGene.txt", header=F)
anno_genes <- anno_genes[,c(3:6,13)]
colnames(anno_genes) <- c("Chr","Strand","Start","End","GeneName")
anno_genes$Chr <- sapply(strsplit(anno_genes$Chr,"_"), `[`, 1)
anno_genes$Chr2 <- as.numeric(sapply(strsplit(anno_genes$Chr,"chr"), `[`, 2))
anno_genes <- anno_genes[order(anno_genes$Chr),]
anno_genes <- anno_genes[order(anno_genes$Chr2),]

# Annotate SNPs dataframes with gene info
SNPs_chr6 <- SNPs_chr6[,c(1,4,2:3)]
for(i in 1:nrow(SNPs_chr6)){
  SNPs_chr6$Strand[i] <- with(anno_genes, Strand[GeneName == SNPs_chr6$Gene[i]])
  SNPs_chr6$Start[i] <-with(anno_genes, Start[GeneName == SNPs_chr6$Gene[i]])
  SNPs_chr6$End[i] <- with(anno_genes, End[GeneName == SNPs_chr6$Gene[i]])}

SNPs_chr5 <- SNPs_chr5[,c(1,4,2:3)]
for(i in 1:nrow(SNPs_chr5)){
  SNPs_chr5$Strand[i] <- with(anno_genes, Strand[GeneName == SNPs_chr5$Gene[i]])
  SNPs_chr5$Start[i] <-with(anno_genes, Start[GeneName == SNPs_chr5$Gene[i]])
  SNPs_chr5$End[i] <- with(anno_genes, End[GeneName == SNPs_chr5$Gene[i]])}

# Annotate results with gene info
res_chr6 <- res_GB_DS$SAIGE_SBA_DS_chr6
res_chr6$Region <- sapply(strsplit(res_chr6$Region,","), `[`, 1)
gf_chr6$V1 <- sapply(strsplit(gf_chr6$V1,","), `[`, 1)
for(i in 1:nrow(res_chr6)){
  res_chr6$Start[i] <-with(anno_genes, Start[GeneName == res_chr6$Region[i]])
  res_chr6$End[i] <- with(anno_genes, End[GeneName == res_chr6$Region[i]])
  print(i)
  res_chr6$N_Var[i] <- with(gf_chr6, N_Var[(V1 == res_chr6$Region[i])&(V2 == "var")])}
res_chr6 <- res_chr6[,c(1,14:16,2:13)]
res_chr6_reg <- subset(res_chr6, subset=((Start >= 143863601)&(End <=173863601)))
rownames(res_chr6_reg) <- 1:nrow(res_chr6_reg)

res_chr5 <- res_GB_DS$SAIGE_SBA_DS_chr5
res_chr5$Region <- sapply(strsplit(res_chr5$Region,","), `[`, 1)
gf_chr5$V1 <- sapply(strsplit(gf_chr5$V1,","), `[`, 1)
for(i in 1:nrow(res_chr5)){
  res_chr5$Start[i] <-with(anno_genes, Start[GeneName == res_chr5$Region[i]])
  res_chr5$End[i] <- with(anno_genes, End[GeneName == res_chr5$Region[i]])
  print(i)
  res_chr5$N_Var[i] <- with(gf_chr5, N_Var[(V1 == res_chr5$Region[i])&(V2 == "var")])}
res_chr5 <- res_chr5[,c(1,14:16,2:13)]
res_chr5_reg <- subset(res_chr5, subset=((Start >= 140749194)&(End <=188279842)))
rownames(res_chr5_reg) <- 1:nrow(res_chr5_reg)
#-------------------------------------------------------------------------------------
# CREATE GENES DATAFRAME, A SUBSET OF THE RESULTS TABLE WITH 1 METRICS ROW PER GENE
# Get gene annotation
GENES_chr6 <- data.frame(unique(res_chr6$Region))
colnames(GENES_chr6) <- "Gene"
for(i in 1:nrow(GENES_chr6)){
  GENES_chr6$Start[i] <-with(anno_genes, Start[GeneName == GENES_chr6$Gene[i]])
  GENES_chr6$End[i] <- with(anno_genes, End[GeneName == GENES_chr6$Gene[i]])
  GENES_chr6$N_Var[i] <- with(gf_chr6, N_Var[(V1 == GENES_chr6$Gene[i])&(V2 == "var")])}
# Get group with best pvalue
for(i in 1:nrow(GENES_chr6)){
 genename <- GENES_chr6[i,1]
 res_tmp <- subset(res_chr6, subset=(Region == genename))
 rownames(res_tmp) <- 1:nrow(res_tmp)
 if(nrow(res_tmp) == 1){GENES_chr6[i,5:16] <- res_tmp[1,5:16]}
 else{
   minP <- min(res_tmp$Pvalue, na.rm=T) 
   index_min <- which(res_tmp$Pvalue == minP)
   if(length(index_min)>1){print(index_min)}
   GENES_chr6[i,5:16] <- res_tmp[index_min[1],5:16]}}
sum(is.na(GENES_chr6$Pvalue))
# [1] 1
ind <- which(is.na(GENES_chr6$Pvalue))
GENES_chr6[ind,1:7]
#      Gene     Start       End N_Var                   Group max_MAF Pvalue
# 492 DEFB114 49927961 49931877     2 missense;lof;synonymous     0.5     NA
subset(res_chr6, subset=(Region == "DEFB114"))
subset(gf_chr6, subset=(V1 == "DEFB114"))
## Remove genes with NA pvalues
GENES_chr6 <- subset(GENES_chr6, subset=(!is.na(Pvalue)))
#-----------------------------------------------------------------------------------
# Get gene annotation
GENES_chr5 <- data.frame(unique(res_chr5$Region))
colnames(GENES_chr5) <- "Gene"
for(i in 1:nrow(GENES_chr5)){
  GENES_chr5$Start[i] <-with(anno_genes, Start[GeneName == GENES_chr5$Gene[i]])
  GENES_chr5$End[i] <- with(anno_genes, End[GeneName == GENES_chr5$Gene[i]])
  GENES_chr5$N_Var[i] <- with(gf_chr5, N_Var[(V1 == GENES_chr5$Gene[i])&(V2 == "var")])}
# Get group with best pvalue
for(i in 1:nrow(GENES_chr5)){
  genename <- GENES_chr5[i,1]
  res_tmp <- subset(res_chr5, subset=(Region == genename))
  rownames(res_tmp) <- 1:nrow(res_tmp)
  if(nrow(res_tmp) == 1){GENES_chr5[i,5:16] <- res_tmp[1,5:16]}
  else{
    minP <- min(res_tmp$Pvalue, na.rm=T) 
    index_min <- which(res_tmp$Pvalue == minP)
    if(length(index_min)>1){print(index_min)}
    GENES_chr5[i,5:16] <- res_tmp[index_min[1],5:16]}}
sum(is.na(GENES_chr5$Pvalue))
# [1] 3
ind <- which(is.na(GENES_chr5$Pvalue))
GENES_chr5[ind,1:7]
#      Gene     Start       End N_Var                   Group max_MAF Pvalue
# 365 ISOC1 128430441 128449719     2 missense;lof;synonymous     0.5     NA - 1 line, 2 synonymous variants
# 594 MFAP3 153418526 153437009     3 missense;lof;synonymous     0.5     NA - 1 line, 3 synonymous variants
# 744 OR2V2 180581942 180582890     2                    <NA>      NA     NA - 3 lines, 2 missense variants
subset(res_chr5, subset=(Region == "OR2V2"))
subset(gf_chr5, subset=(V1 == "OR2V2"))
## Remove genes with NA pvalues
GENES_chr5 <- subset(GENES_chr5, subset=(!is.na(Pvalue)))

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
# WE HAVE 905 GENES ON CHR6
# Region of interest on CHR6: 143,863,601 - 173,863,601
GENES_chr6_reg <- subset(GENES_chr6, subset=((Start >= 143863601)&(End <=173863601)))
rownames(GENES_chr6_reg) <- 1:nrow(GENES_chr6_reg)
# 116/905 in region
# Subset of variants in that region
SNPs_chr6_reg <- subset(SNPs_chr6, subset=(Gene %in% GENES_chr6_reg$Gene))
# 1,090/7,071 in region

# WE HAVE 745 GENES ON CHR5
# Region of interest on CHR5: 140,749,194 - 188,279,842
GENES_chr5_reg <- subset(GENES_chr5, subset=((Start >= 140749194)&(End <=188279842)))
rownames(GENES_chr5_reg) <- 1:nrow(GENES_chr5_reg)
# 246/745 in region
# Subset of variants in that region
SNPs_chr5_reg <- subset(SNPs_chr5, subset=(Gene %in% GENES_chr5_reg$Gene))
# 1,775/5,408 in region
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-----------------------------------------
#     Calculate GB lambda GENOME-WIDE
#-----------------------------------------
for(i in 1:22){
  res_GB_DS[[i]]$Chr <- rep(i,nrow(res_GB_DS[[i]]))
  res_GB_DS[[i]] <- res_GB_DS[[i]][,c(14,1:13)]}
# Merge results per chromosome into single dataframe
res_GB_ALL <- bind_rows(res_GB_DS)
res_GB_ALL_noCauchy <- subset(res_GB_ALL, subset=(Group != "Cauchy"))
write.table(res_GB_ALL, file="./output_data/SAIGE/SLE/Results/SAIGE_SBA_ALL_CHR", sep="\t", col.names=T, row.names=F, quote=F)

# SKAT-O
qchisq(1-median(res_GB_ALL_noCauchy$Pvalue, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 1.021603

# BURDEN
qchisq(1-median(res_GB_ALL_noCauchy$Pvalue_Burden, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 1.055074

# SKAT
qchisq(1-median(res_GB_ALL_noCauchy$Pvalue_SKAT, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 1.043386
#-----------------------------------------
#     QQ Plots - GENOME-WIDE
#-----------------------------------------
qqplot_test(res_GB_ALL_noCauchy$Pvalue) + theme_minimal() +
  labs(title="Q-Q Plot for gene burden analysis: SKAT-O P-values")
qqplot_test(res_GB_ALL_noCauchy$Pvalue_SKAT) + theme_minimal() +
  labs(title="Q-Q Plot for gene burden analysis: SKAT P-values")
qqplot_test(res_GB_ALL_noCauchy$Pvalue_Burden) + theme_minimal() +
  labs(title="Q-Q Plot for gene burden analysis: Burden P-values")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# MitoCarta 3.0 Genes
mitoGenes <- read.delim("./input_data/Human_MitoCarta_GeneList.txt", header=F)
GENES_chr6_reg$Mito <- ifelse(GENES_chr6_reg$Gene %in% mitoGenes$V1,"Yes","No")
sum(GENES_chr6_reg$Mito == "Yes")
# [1] 10 "Mito" genes in chr6 region
GENES_chr5_reg$Mito <- ifelse(GENES_chr5_reg$Gene %in% mitoGenes$V1,"Yes","No")
sum(GENES_chr5_reg$Mito == "Yes")
# [1] 7 "Mito" genes in chr5 region
#------------------------------------------------------------------------------------
ggplot(GENES_chr6_reg,aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  geom_hline(yintercept=c(-log10(0.05/10),-log10(0.05/116)), color = c("orange", "midnightblue"), linetype = c("dashed", "solid"), lwd = 1) +
  scale_color_manual(values=c("midnightblue", "orange")) +
  geom_label_repel(data = subset(GENES_chr6_reg, subset=(as.numeric(Pvalue) < 0.05|Mito == "Yes")), size = 4, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  theme_minimal() +
  labs(title="SAIGE-GENE: SLE Associations - SKAT-O Pvalues", x ="Chromosome 6 position", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 13), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.title.y = element_markdown(), axis.text.x = element_text(angle = 0, size = 11))

ggplot(GENES_chr5_reg,aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  geom_hline(yintercept=c(-log10(0.05/7),-log10(0.05/246)), color = c("orange", "midnightblue"), linetype = c("dashed", "solid"), lwd = 1) +
  scale_color_manual(values=c("midnightblue", "orange")) +
  geom_label_repel(data = subset(GENES_chr5_reg, subset=(as.numeric(Pvalue) < 0.05|Mito == "Yes")), size = 4, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  theme_minimal() +
  labs(title="SAIGE-GENE: SLE Associations - SKAT-O Pvalues", x ="Chromosome 5 position", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 13), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.title.y = element_markdown(), axis.text.x = element_text(angle = 0, size = 11))
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                        Manhattan & Gene Burden Plots (merged)
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------
#            Chromosome 6
#--------------------------------------
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr6_reg)  # 4.4e-06
MANHATTAN6 <- ggplot(res_ALL001_DS_chr6_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value))) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("gray","#fe7f2d"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  scale_color_manual(values = rep(c("#233d4d", "#233d4d"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=NULL) +
  theme_bw() +
  xlab(NULL) +
  ylab(bquote(-log[10](P))) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
        ggtitle("A. Single Variant Analysis")

GENES6 <- ggplot(GENES_chr6_reg, aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  geom_hline(yintercept=-log10(0.05/114), color = "#fe7f2d", linetype = "solid", lwd = 1) +
  scale_color_manual(values=c("#233d4d", "skyblue3")) +
  scale_x_continuous(breaks=seq(140000000,175000000,by=5000000), labels=scales::comma) +
  #guides(col = "none") +
  geom_label_repel(data = subset(GENES_chr6_reg, subset=(-log10(Pvalue) > 2|Mito == "Yes")),  fill = alpha(c("white"),0.8), size = 3, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  theme_bw() +
  xlab("Position (bp)") +
  ylab(bquote(-log[10](P))) +
  theme(legend.position = "none",
    axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.2),
    axis.title.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 12, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
    ggtitle("B. Gene Burden Analysis")

grid.arrange(MANHATTAN6, GENES6, nrow=2)

min(res_ALL001_DS_chr6_reg$p.value)
which(res_ALL001_DS_chr6_reg$p.value == 8.285175e-08)
res_ALL001_DS_chr6_reg[20939,]
# 6:151292437:G:A --- MTHFD1L gene
#--------------------------------------
#            Chromosome 5
#--------------------------------------
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr5_reg)  # 3.3e-06
MANHATTAN5 <- ggplot(res_ALL001_DS_chr5_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value))) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("gray","#fe7f2d"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  scale_color_manual(values = rep(c("#233d4d", "#233d4d"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=NULL) +
  theme_bw() +
  xlab(NULL) +
  ylab(bquote(-log[10](P))) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("A. Single Variant Analysis")

GENES5 <- ggplot(GENES_chr5_reg, aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  geom_hline(yintercept=-log10(0.05/247), color = "#fe7f2d", linetype = "solid", lwd = 1) +
  scale_color_manual(values=c("#233d4d", "skyblue3")) +
  scale_x_continuous(breaks=seq(140000000,190000000,by=10000000), labels=scales::comma) +
  #guides(col = "none") +
  geom_label_repel(data = subset(GENES_chr5_reg, subset=(-log10(Pvalue) > 2|Mito == "Yes")),  fill = alpha(c("white"),0.8), size = 3, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  theme_bw() +
  xlab("Position (bp)") +
  ylab(bquote(-log[10](P))) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ggtitle("B. Gene Burden Analysis")

grid.arrange(MANHATTAN5, GENES5, nrow=2)

min(res_ALL001_DS_chr5_reg$p.value)
which(res_ALL001_DS_chr5_reg$p.value == 4.104709e-05)
res_ALL001_DS_chr5_reg[75247,]
# 5:167983347:T:C
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                        COLOUR MANHATTAN BY GENE
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
TULP4_var <- as.matrix(subset(gf_chr6, select=Variants, subset=(V1=="TULP4" & V2=="var")))
TULP4_var <- unlist(strsplit(TULP4_var, split="\\|"))
TULP4_var <- sapply(strsplit(TULP4_var, ":"),`[`, 2)
## 22 variants
FRMD1_var <- as.matrix(subset(gf_chr6, select=Variants, subset=(V1=="FRMD1" & V2=="var")))
FRMD1_var <- unlist(strsplit(FRMD1_var, split="\\|"))
FRMD1_var <- sapply(strsplit(FRMD1_var, ":"),`[`, 2)
# 14 variants
MTHFD1L_var <- as.matrix(subset(gf_chr6, select=Variants, subset=(V1=="MTHFD1L" & V2=="var")))
MTHFD1L_var <- unlist(strsplit(MTHFD1L_var, split="\\|"))
MTHFD1L_var <- sapply(strsplit(MTHFD1L_var, ":"),`[`, 2)
# 5 variants

res_ALL001_DS_chr6_reg$GeneLabel <- "Other"
for(i in 1:nrow(res_ALL001_DS_chr6_reg)){
  if(res_ALL001_DS_chr6_reg$POS[i] %in% TULP4_var){res_ALL001_DS_chr6_reg$GeneLabel[i] <- "TULP4"}
  else if(res_ALL001_DS_chr6_reg$POS[i] %in% FRMD1_var){res_ALL001_DS_chr6_reg$GeneLabel[i] <- "FRMD1"}
  else if(res_ALL001_DS_chr6_reg$POS[i] %in% MTHFD1L_var){res_ALL001_DS_chr6_reg$GeneLabel[i] <- "MTHFD1L"}}

# GWAS was done with MAF > 0.01 cut-off so some variants were removed
sum(res_ALL001_DS_chr6_reg$GeneLabel=="TULP4")
#[1] 7
sum(res_ALL001_DS_chr6_reg$GeneLabel=="FRMD1")
##[1] 10
sum(res_ALL001_DS_chr6_reg$GeneLabel=="MTHFD1L")
##[1] 2

sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr6_reg)  # 4.4e-06
ggplot(data = res_ALL001_DS_chr6_reg, aes(x = POS, y = -log10(p.value))) +
  geom_hline(yintercept =-log10(sig), color = "gray", linetype = "dashed", lwd = 1) +
  geom_point(aes(color = GeneLabel, alpha = GeneLabel)) +
  scale_color_manual(values=setNames(c("#233d4d", "yellow", "cyan", "orange"),c("Other","TULP4","FRMD1", "MTHFD1L"))) +
  scale_alpha_manual(values=setNames(c(0.07, 1, 1, 1),c("Other","TULP4","FRMD1", "MTHFD1L"))) +
  scale_x_continuous(breaks=seq(140000000,175000000,by=5000000), labels=scales::comma) +
  theme_minimal() +
  labs(
    #title="SLE Associations",
    x ="",
    y ="-log<sub>10</sub>(P)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        #legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
SPINK5_var <- as.matrix(subset(gf_chr5, select=Variants, subset=(V1=="SPINK5" & V2=="var")))
SPINK5_var <- unlist(strsplit(SPINK5_var, split="\\|"))
SPINK5_var <- sapply(strsplit(SPINK5_var, ":"),`[`, 2)
## 25 variants

PANK3_var <- as.matrix(subset(gf_chr5, select=Variants, subset=(V1=="PANK3" & V2=="var")))
PANK3_var <- unlist(strsplit(PANK3_var, split="\\|"))
PANK3_var <- sapply(strsplit(PANK3_var, ":"),`[`, 2)
## 2 variants

res_ALL001_DS_chr5_reg$GeneLabel <- "Other"
for(i in 1:nrow(res_ALL001_DS_chr5_reg)){
  if(res_ALL001_DS_chr5_reg$POS[i] %in% SPINK5_var){res_ALL001_DS_chr5_reg$GeneLabel[i] <- "SPINK5"}
  else if(res_ALL001_DS_chr5_reg$POS[i] %in% PANK3_var){res_ALL001_DS_chr5_reg$GeneLabel[i] <- "PANK3"}}

# GWAS was done with MAF > 0.01 cut-off so some variants were removed
sum(res_ALL001_DS_chr5_reg$GeneLabel=="SPINK5")
#[1] 19
sum(res_ALL001_DS_chr5_reg$GeneLabel=="PANK3")
#[1] 0

sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr5_reg)  # 3.3e-06
ggplot(data = res_ALL001_DS_chr5_reg, aes(x = POS, y = -log10(p.value))) +
  geom_hline(yintercept =-log10(sig), color = "gray", linetype = "dashed", lwd = 1) +
  geom_point(aes(color = GeneLabel, alpha = GeneLabel)) +
  scale_color_manual(values=setNames(c("#233d4d", "orange"),c("Other","SPINK5"))) +
  scale_alpha_manual(values=setNames(c(0.07, 1),c("Other","SPINK5"))) +
  scale_x_continuous(breaks=seq(140000000,190000000,by=10000000), labels=scales::comma) +
  theme_minimal() +
  labs(
    #title="SLE Associations",
    x ="",
    y ="-log<sub>10</sub>(P)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        #legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))
