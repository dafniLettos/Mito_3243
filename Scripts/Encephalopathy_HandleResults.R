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

# Read-in Group File
gf_chr7 <- read.table("./output_data/Geno_Imputed/Group_File_chr7", header=F, col.names = paste0("V",seq_len(max(count.fields("./output_data/Geno_Imputed/Group_File_chr7", sep = '\t')))), fill = TRUE)
gf_chr11 <- read.table("./output_data/Geno_Imputed/Group_File_chr11", header=F, col.names = paste0("V",seq_len(max(count.fields("./output_data/Geno_Imputed/Group_File_chr7", sep = '\t')))), fill = TRUE)
gf_chr13 <- read.table("./output_data/Geno_Imputed/Group_File_chr13", header=F, col.names = paste0("V",seq_len(max(count.fields("./output_data/Geno_Imputed/Group_File_chr7", sep = '\t')))), fill = TRUE)

# Read-in SAIGE-GENE Results
res_files_ds <- paste0(rep("SAIGE_SBA_DS_chr",22),c(1:22))
res_GB_DS <- list()
for(file in res_files_ds){
  print(paste0("./output_data/SAIGE/Encephalopathy/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/Encephalopathy/Results/",file), header=T)
  res_GB_DS[[file]] <- data}

# Read-in SAIGE Results (GWAS)
res_files_res_ALL001_ds <- paste0(rep("SAIGE_SVA_maf001_chr",22),c(1:22),rep(".txt",22))
res_SVA001_DS <- list()
for(file in res_files_res_ALL001_ds){
  print(paste0("./output_data/SAIGE/Encephalopathy/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/Encephalopathy/Results/",file), header=T)
  res_SVA001_DS[[file]] <- data}

res_files_res_ALL005_ds <- paste0(rep("SAIGE_SVA_maf005_chr",22),c(1:22),rep(".txt",22))
res_SVA005_DS <- list()
for(file in res_files_res_ALL005_ds){
  print(paste0("./output_data/SAIGE/Encephalopathy/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/Encephalopathy/Results/",file), header=T)
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
# Merge results per chromosome into single dataframe
res_ALL001_DS_ALL <- bind_rows(res_SVA001_DS)
# 7,586,961 SNPs
write.table(res_ALL001_DS_ALL, file="./output_data/SAIGE/Encephalopathy/Results/SAIGE_SVA_DS_ALL_MAF001.txt", sep="\t", col.names=T, row.names=F, quote=F)
hist(res_ALL001_DS_ALL$p.value)
res_ALL001_DS_chr7 <- subset(res_ALL001_DS_ALL, subset=(CHR == 7))
res_ALL001_DS_chr7_reg <- subset(res_ALL001_DS_chr7, subset=((POS >= 91468550)&(POS <=121468550)))
hist(res_ALL001_DS_chr7_reg$p.value)
sum(res_ALL001_DS_chr7_reg$p.value < 5e-08)
#[1] 0
res_ALL001_DS_chr7_reg$p.value_log <- -log10(as.numeric(res_ALL001_DS_chr7_reg$p.value))
res_ALL001_DS_chr7_reg$SNP_ID <- paste(res_ALL001_DS_chr7_reg$MarkerID,res_ALL001_DS_chr7_reg$Allele1,res_ALL001_DS_chr7_reg$Allele2,sep=":")

# Merge results per chromosome into single dataframe
res_ALL005_DS_ALL <- bind_rows(res_SVA005_DS)
# 5,464,655 SNPs
write.table(res_ALL005_DS_ALL, file="./output_data/SAIGE/Encephalopathy/Results/SAIGE_SVA_DS_ALL_MAF005.txt", sep="\t", col.names=T, row.names=F, quote=F)
hist(res_ALL005_DS_ALL$p.value)
res_ALL005_DS_chr7 <- subset(res_ALL005_DS_ALL, subset=(CHR == 7))
res_ALL005_DS_chr7_reg <- subset(res_ALL005_DS_chr7, subset=((POS >= 91468550)&(POS <=121468550)))
hist(res_ALL005_DS_chr7_reg$p.value)
sum(res_ALL005_DS_chr7_reg$p.value < 5e-08)
#[1] 0
res_ALL005_DS_chr7_reg$p.value_log <- -log10(as.numeric(res_ALL005_DS_chr7_reg$p.value))
res_ALL005_DS_chr7_reg$SNP_ID <- paste(res_ALL005_DS_chr7_reg$MarkerID,res_ALL005_DS_chr7_reg$Allele1,res_ALL005_DS_chr7_reg$Allele2,sep=":")
#----------------------------------
#         Manhattan Plot
#----------------------------------
## Suggestive significance level for the proportion of the genome tested:
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr7_reg)  # 5.4e-06  -log10(sig) --- 5.26
ggplot(res_ALL001_DS_chr7_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL001_DS_chr7_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 2, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("midnightblue", "midnightblue"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(92000000,120000000,by=5000000), labels=scales::comma) +
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

sig <- (nrow(res_ALL005_DS_ALL) * 5e-08) / nrow(res_ALL005_DS_chr7_reg)  # 5.7e-06  -log10(sig) --- 5.25
ggplot(res_ALL005_DS_chr7_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL005_DS_chr7_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 2, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("midnightblue", "midnightblue"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(92000000,120000000,by=5000000), labels=scales::comma) +
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
#-----------------------------------------------------------------
#           1. Calculate genomic inflation factor (lambda)
#-----------------------------------------------------------------
# Lambda (λ) quantifies the deviation of observed test statistics from what's expected under the null hypothesis.
# A genomic inflation factor of 1 indicates no bias, while values greater than 1 suggest inflation,
# potentially due to factors like population structure or cryptic relatedness
# Convert p-values to z-scores: z_scores <- qnorm(Pvalue / 2, lower.tail = FALSE) 
# OR CALCULATE chi-squared stats STRAIGHT AWAY: chi_squared_stats <- z_scores^2
# Convert z-scores to chi-squared stats: chi_squared_stats <- qchisq(Pvalue, df = 1, lower.tail = FALSE)
# Calculate the genomic inflation factor (lambda): λ = OBSERVED MEDIAN/EXPECTED MEDIAN of chi-squared statistics
# lambda <- median(chi_squared_stats) / qchisq(0.5, df=1)
#-------------------------------------------------------------------------------------
qchisq(1-median(res_ALL001_DS_ALL$p.value),1)/qchisq(0.5,1)
#[1] 1.005991
qchisq(1-median(res_ALL005_DS_ALL$p.value),1)/qchisq(0.5,1)
#[1] 0.9968078
#-----------------------------------------------------------------
#                       2. Q-Q PLOTS
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
length(unique(res_GB_DS$SAIGE_SBA_DS_chr7$Region))
# 762 GENES

# Initialise SNPs dataframe that will contain variant information
SNPs <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(SNPs) <- c("Variant","Annotation","Gene")
# These column in the group file 
gf_chr7$N_Var <- ""
gf_chr7$Variants <- ""

# FOR LOOP to get some extra columns of info for variants per gene in the group file & create the SNPs dataframe
for(i in 1:nrow(gf_chr7)){
  if(gf_chr7[i,2] == "var"){
    empty <- sum(gf_chr7[i,] == "") - 2
    nVar <- ncol(gf_chr7) - empty - 4
    gene <- gf_chr7[i,1]
    #firstVar <- gf_chr7[i,3]
    #lastVar <- gf_chr7[i,(nVar+2)]
    
    gf_chr7$N_Var[i] <- nVar
    gf_chr7$Variants[i] <- paste(gf_chr7[i,3:(nVar+2)], collapse="|")
    #print(i)
    SNPs_tmp <- data.frame(
            Variant = c(t(gf_chr7[i,3:(nVar+2)])),
            Annotation = c(t(gf_chr7[i+1,3:(nVar+2)])),
            Gene = rep(gene, times=nVar))
    SNPs <- rbind(SNPs, SNPs_tmp)}
  else{
    gf_chr7$N_Var[i] <- ""
    gf_chr7$Variants[i] <- ""}
}
# Add POSITION column to SNPs dataframe
SNPs$SNP_POS <- sapply(strsplit(SNPs$Variant,":"), `[`, 2)
SNPs$Gene <- sapply(strsplit(SNPs$Gene,","), `[`, 1)

## SOS!! Add GENE START for POSITION HERE
# Will be using the file from ANNOVAR
anno_genes <- fread("./annovar/humandb/hg19_refGene.txt", header=F)
anno_genes <- anno_genes[,c(3:6,13)]
colnames(anno_genes) <- c("Chr","Strand","Start","End","GeneName")
anno_genes$Chr <- sapply(strsplit(anno_genes$Chr,"_"), `[`, 1)
anno_genes$Chr2 <- as.numeric(sapply(strsplit(anno_genes$Chr,"chr"), `[`, 2))
anno_genes <- anno_genes[order(anno_genes$Chr),]
anno_genes <- anno_genes[order(anno_genes$Chr2),]

# Annotate SNPs dataframe with gene info
SNPs <- SNPs[,c(1,4,2:3)]
for(i in 1:nrow(SNPs)){
  SNPs$Strand[i] <- with(anno_genes, Strand[GeneName == SNPs$Gene[i]])
  SNPs$Start[i] <-with(anno_genes, Start[GeneName == SNPs$Gene[i]])
  SNPs$End[i] <- with(anno_genes, End[GeneName == SNPs$Gene[i]])}

# Annotate results with gene info
res_chr7 <- res_GB_DS$SAIGE_SBA_DS_chr7
res_chr7$Region <- sapply(strsplit(res_chr7$Region,","), `[`, 1)
gf_chr7$V1 <- sapply(strsplit(gf_chr7$V1,","), `[`, 1)
for(i in 1:nrow(res_chr7)){
  res_chr7$Start[i] <-with(anno_genes, Start[GeneName == res_chr7$Region[i]])
  res_chr7$End[i] <- with(anno_genes, End[GeneName == res_chr7$Region[i]])
  #print(i)
  res_chr7$N_Var[i] <- with(gf_chr7, N_Var[(V1 == res_chr7$Region[i])&(V2 == "var")])}
res_chr7 <- res_chr7[,c(1,14:16,2:13)]

res_chr7_reg <- subset(res_chr7, subset=((Start >= 91468550)&(End <=121468550)))
rownames(res_chr7_reg) <- 1:nrow(res_chr7_reg)
#-------------------------------------------------------------------------------------
# CREATE GENES DATAFRAME, A SUBSET OF THE RESULTS TABLE WITH 1 METRICS ROW PER GENE
# Get gene annotation
GENES <- data.frame(unique(res_chr7$Region))
colnames(GENES) <- "Gene"
for(i in 1:nrow(GENES)){
  GENES$Start[i] <-with(anno_genes, Start[GeneName == GENES$Gene[i]])
  GENES$End[i] <- with(anno_genes, End[GeneName == GENES$Gene[i]])
  GENES$N_Var[i] <- with(gf_chr7, N_Var[(V1 == GENES$Gene[i])&(V2 == "var")])}
# Get group with best pvalue
for(i in 1:nrow(GENES)){
 print(i)
 genename <- GENES[i,1]
 res_tmp <- subset(res_chr7, subset=(Region == genename))
 rownames(res_tmp) <- 1:nrow(res_tmp)
 if(nrow(res_tmp) == 1){GENES[i,5:16] <- res_tmp[1,5:16]}
 else{
   minP <- min(res_tmp$Pvalue, na.rm=T) 
   index_min <- which(res_tmp$Pvalue == minP)
   if(length(index_min)>1){print(index_min)}
   GENES[i,5:16] <- res_tmp[index_min[1],5:16]}}

sum(is.na(GENES$Pvalue))
# [1] 2
ind <- which(is.na(GENES$Pvalue))
GENES[ind,1:7]
# Gene     Start       End N_Var Group max_MAF Pvalue
# 378 PON2  95034173  95064340     3  <NA>      NA     NA
# 564 CALU 128379412 128413454     2  <NA>      NA     NA

## Remove genes with NA pvalues
GENES <- subset(GENES, subset=(!is.na(Pvalue)))

# WE HAVE 756 GENES ON CHR7
# Region of interest on CHR7: 91,468,550 - 121,468,550
GENES_reg <- subset(GENES, subset=((Start >= 91468550)&(End <=121468550)))
rownames(GENES_reg) <- 1:nrow(GENES_reg)
# 178/756 in region

# Subset of variants in that region
SNPs_7q22 <- subset(SNPs, subset=(Gene %in% GENES_reg$Gene))
# 1,202/5,687 in region
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
res_GB_DS_ALL <- bind_rows(res_GB_DS)
res_GB_DS_ALL_noCauchy <- subset(res_GB_DS_ALL, subset=(Group != "Cauchy"))
write.table(res_GB_DS_ALL, file="./output_data/SAIGE/Encephalopathy/Results/SAIGE_SBA_DS_ALL_CHR", sep="\t", col.names=T, row.names=F, quote=F)

# SKAT-O
qchisq(1-median(res_GB_DS_ALL_noCauchy$Pvalue, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 0.9851934

# BURDEN
qchisq(1-median(res_GB_DS_ALL_noCauchy$Pvalue_Burden, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 1.016593

# SKAT
qchisq(1-median(res_GB_DS_ALL_noCauchy$Pvalue_SKAT, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 1.014835
#-----------------------------------------
#     QQ Plots - GENOME-WIDE
#-----------------------------------------
qqplot_test(res_GB_DS_ALL_noCauchy$Pvalue) + theme_minimal() +
  labs(title="Q-Q Plot for gene burden analysis: SKAT-O P-values")
qqplot_test(res_GB_DS_ALL_noCauchy$Pvalue_SKAT) + theme_minimal() +
  labs(title="Q-Q Plot for gene burden analysis: SKAT P-values")
qqplot_test(res_GB_DS_ALL_noCauchy$Pvalue_Burden) + theme_minimal() +
  labs(title="Q-Q Plot for gene burden analysis: Burden P-values")
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# MitoCarta 3.0 Genes
mitoGenes <- read.delim("./input_data/Human_MitoCarta_GeneList.txt", header=F)
GENES_reg$Mito <- ifelse(GENES_reg$Gene %in% mitoGenes$V1,"Yes","No")
sum(GENES_reg$Mito == "Yes")
# [1] 13
#------------------------------------------------------------------------------------
ggplot(GENES_reg,aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  #geom_hline(yintercept= -log10(0.05/149), color = "darkolivegreen", linetype = "solid", lwd = 1) +
  geom_hline(yintercept=c(-log10(0.05/13),-log10(0.05/178)), color = c("orange", "midnightblue"), linetype = c("dashed", "solid"), lwd = 1) +
  scale_color_manual(values=c("midnightblue", "orange")) +
  #guides(col = "none") +
  geom_label_repel(data = subset(GENES_reg, subset=(as.numeric(Pvalue) < 0.05|Mito == "Yes")), size = 4, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  theme_minimal() +
  labs(title="SAIGE-GENE: Encephalopathy Associations - SKAT-O Pvalues", x ="Chromosome 7 position", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 13), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.title.y = element_markdown(), axis.text.x = element_text(angle = 0, size = 11))
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                         LINKAGE PAPER PLOT
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr7_reg)  # 5.44e-06 -log10(sig) --- 5.3
MANHATTAN7 <- ggplot(res_ALL001_DS_chr7_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value))) +
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

GENES7 <- ggplot(GENES_reg, aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  geom_hline(yintercept=-log10(0.05/178), color = "#fe7f2d", linetype = "solid", lwd = 1) +
  scale_color_manual(values=c("#233d4d", "skyblue3")) +
  scale_x_continuous(breaks=seq(10000000,150000000,by=5000000), labels=scales::comma) +
  #guides(col = "none") +
  geom_label_repel(data = subset(GENES_reg, subset=(-log10(Pvalue) > 1.7|Mito == "Yes")),  fill = alpha(c("white"),0.8), size = 3, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
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

grid.arrange(MANHATTAN7, GENES7, nrow=2)

min(res_ALL001_DS_chr7_reg$p.value)
which(res_ALL001_DS_chr7_reg$p.value == 3.990941e-06)
res_ALL001_DS_chr7_reg[11755,]
# 7:97072414:G:T --- rs62497093

subset(GENES_reg, subset=(Gene == "PIK3CG"))# Pvalue = 0.005107425
subset(GENES_reg, subset=(Gene == "PLOD3"))# Pvalue = 0.01467675
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                        EXPORT RESULTS FOR SUPPLEMENTARY
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
res_ALL001_DS_chr7_reg <- res_ALL001_DS_chr7_reg[,c(1:3,25,4:24)]
colnames(res_ALL001_DS_chr7_reg)[4] <- "CHR:POS:A1:A2"
write.table(res_ALL001_DS_chr7_reg, file="./output_data/SAIGE/Encephalopathy/SupTable2_SAIGE_ds_chr7q22.tab", sep="\t", row.names=F, col.names=T, quote=F)

sum(is.na(GENES_reg[,7]))
which(is.na(GENES_reg[,6]))
# 168
GENES_reg[168,6] <- 0.5
which(GENES_reg$Group == "Cauchy")
GENES_reg[168,5] <- "missense;lof;synonymous"
write.table(GENES_reg, file="./output_data/SAIGE/Encephalopathy/SupTable3_SAIGEGENE_ds_chr7q22.tab", sep="\t", row.names=F, col.names=T, quote=F)
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                        COLOUR MANHATTAN BY GENE
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
PIK3CG_var <- as.matrix(subset(gf_chr7, select=Variants, subset=(V1=="PIK3CG" & V2=="var")))
PIK3CG_var <- unlist(strsplit(PIK3CG_var, split="\\|"))
PIK3CG_var <- sapply(strsplit(PIK3CG_var, ":"),`[`, 2)
# 10 variants
PLOD3_var <- as.matrix(subset(gf_chr7, select=Variants, subset=(V1=="PLOD3" & V2=="var")))
PLOD3_var <- unlist(strsplit(PLOD3_var, split="\\|"))
PLOD3_var <- sapply(strsplit(PLOD3_var, ":"),`[`, 2)
# 15 variants

res_ALL001_DS_chr7_reg$GeneLabel <- "Other"
for(i in 1:nrow(res_ALL001_DS_chr7_reg)){
  if(res_ALL001_DS_chr7_reg$POS[i] %in% PIK3CG_var){res_ALL001_DS_chr7_reg$GeneLabel[i] <- "PIK3CG"}
  else if(res_ALL001_DS_chr7_reg$POS[i] %in% PLOD3_var){res_ALL001_DS_chr7_reg$GeneLabel[i] <- "PLOD3"}}

# GWAS was done with MAF > 0.01 cut-off so some variants were removed
sum(res_ALL001_DS_chr7_reg$GeneLabel=="PIK3CG")
#[1] 8
sum(res_ALL001_DS_chr7_reg$GeneLabel=="PLOD3")
#[1] 2

sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr7_reg)  # 5.4e-06
ggplot(data = res_ALL001_DS_chr7_reg, aes(x = POS, y = -log10(p.value))) +
  geom_hline(yintercept =-log10(sig), color = "gray", linetype = "dashed", lwd = 1) +
  geom_point(aes(color = GeneLabel, alpha = GeneLabel)) +
  scale_color_manual(values=setNames(c("#233d4d", "yellow", "cyan"),c("Other","PIK3CG","PLOD3"))) +
  scale_alpha_manual(values=setNames(c(0.07, 1, 1),c("Other","PIK3CG","PLOD3"))) +
  scale_x_continuous(breaks=seq(92000000,120000000,by=5000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
  theme_minimal() +
  labs(
    #title="SAIGE SVA: Encephalopathy Associations",
    x ="",
    y ="-log<sub>10</sub>(P)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        #legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))
#-------------------------------------------------------------------------------------------------------------------------------
#                             Chr1: 1:10048151_G/T : Pval=4.69 × 10^-8 [rs139490379]
#                     Genome-wide significance (?), influencing gene NMNAT1
#                         (nicotinamide mononucleotide adenylyltransferase 1)
#                                         ~  NAD+ Biosynthesis!!!
#                       STANDALONE - NOT AN ACTUAL PEAK - SO, POSSIBLY NOT REAL
#-------------------------------------------------------------------------------------------------------------------------------
res_ALL001_DS_chr1 <- subset(res_ALL001_DS_ALL, subset=(CHR == 1))
res_ALL001_DS_chr1_reg <- subset(res_ALL001_DS_chr1, subset=((POS >= 693731)&(POS <=20800000)))
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr1_reg)  # 6.6e-06
ggplot(res_ALL001_DS_chr1, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=paste(MarkerID,Allele1,Allele2, sep=":"))) +
  geom_hline(yintercept = -log10(5e-08), color = "orange", linetype = "dashed", lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL001_DS_chr1, subset=(-log10(p.value) >= 6)), size = 3, max.overlaps = 10, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("#233d4d", "#233d4d"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(labels=scales::comma) +
  theme_minimal() +
  labs(title="Chromosome 1: Encephalopathy Associations", x ="", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                             Chr13:30976845:C/A
#-------------------------------------------------------------------------------------------------------------------------------
res_ALL001_DS_chr13 <- subset(res_ALL001_DS_ALL, subset=(CHR == 13))
res_ALL001_DS_chr13_reg <- subset(res_ALL001_DS_chr13, subset=((POS >= 10904797)&(POS <=40904797)))
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr13_reg)  # 5.5e-06
ggplot(res_ALL001_DS_chr13_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=paste(MarkerID,Allele1,Allele2,sep=":"))) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("#fe7f2d", "#233d4d"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL001_DS_chr13_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 30, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("#233d4d", "#233d4d"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(10904797,40904797,by=5000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
  theme_minimal() +
  labs(title="Chromosome 13: Encephalopathy Associations", x ="", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))
#----------------------------------------------------------------------------------------------------------------------
# How many genes were tested?
length(unique(res_GB_DS$SAIGE_SBA_DS_chr13$Region))
# 282 GENES
SNPs13 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(SNPs13) <- c("Variant","Annotation","Gene")
gf_chr13$N_Var <- ""
gf_chr13$Variants <- ""
for(i in 1:nrow(gf_chr13)){
  if(gf_chr13[i,2] == "var"){
    empty <- sum(gf_chr13[i,] == "") - 2
    nVar <- ncol(gf_chr13) - empty - 4
    gene <- gf_chr13[i,1]
    gf_chr13$N_Var[i] <- nVar
    gf_chr13$Variants[i] <- paste(gf_chr13[i,3:(nVar+2)], collapse="|")
    #print(i)
    SNPs_tmp <- data.frame(
      Variant = c(t(gf_chr7[i,3:(nVar+2)])),
      Annotation = c(t(gf_chr7[i+1,3:(nVar+2)])),
      Gene = rep(gene, times=nVar))
    SNPs13 <- rbind(SNPs13, SNPs_tmp)}
  else{
    gf_chr7$N_Var[i] <- ""
    gf_chr7$Variants[i] <- ""}
}

SNPs13$SNP_POS <- sapply(strsplit(SNPs13$Variant,":"), `[`, 2)
SNPs13$Gene <- sapply(strsplit(SNPs13$Gene,","), `[`, 1)
SNPs13 <- SNPs13[,c(1,4,2:3)]
for(i in 1:nrow(SNPs13)){
  SNPs13$Strand[i] <- with(anno_genes, Strand[GeneName == SNPs13$Gene[i]])
  SNPs13$Start[i] <-with(anno_genes, Start[GeneName == SNPs13$Gene[i]])
  SNPs13$End[i] <- with(anno_genes, End[GeneName == SNPs13$Gene[i]])}

res_chr13 <- res_GB_DS$SAIGE_SBA_DS_chr13
res_chr13$Region <- sapply(strsplit(res_chr13$Region,","), `[`, 1)
gf_chr13$V1 <- sapply(strsplit(gf_chr13$V1,","), `[`, 1)
for(i in 1:nrow(res_chr13)){
  res_chr13$Start[i] <-with(anno_genes, Start[GeneName == res_chr13$Region[i]])
  res_chr13$End[i] <- with(anno_genes, End[GeneName == res_chr13$Region[i]])
  res_chr13$N_Var[i] <- with(gf_chr13, N_Var[(V1 == res_chr13$Region[i])&(V2 == "var")])}
res_chr13 <- res_chr13[,c(1,14:16,2:13)]
res_chr13_reg <- subset(res_chr13, subset=((Start >= 10904797)&(End <=40904797)))
rownames(res_chr13_reg) <- 1:nrow(res_chr13_reg)
  
GENES13 <- data.frame(unique(res_chr13$Region))
colnames(GENES13) <- "Gene"
for(i in 1:nrow(GENES13)){
  GENES13$Start[i] <-with(anno_genes, Start[GeneName == GENES13$Gene[i]])
  GENES13$End[i] <- with(anno_genes, End[GeneName == GENES13$Gene[i]])
  GENES13$N_Var[i] <- with(gf_chr13, N_Var[(V1 == GENES13$Gene[i])&(V2 == "var")])}
# Get group with best pvalue
for(i in 1:nrow(GENES13)){
  print(i)
  genename <- GENES13[i,1]
  res_tmp <- subset(res_chr13, subset=(Region == genename))
  rownames(res_tmp) <- 1:nrow(res_tmp)
  if(nrow(res_tmp) == 1){GENES13[i,5:16] <- res_tmp[1,5:16]}
  else{
      minP <- min(res_tmp$Pvalue, na.rm=T) 
      index_min <- which(res_tmp$Pvalue == minP)
      if(length(index_min)>1){print(index_min)}
      GENES13[i,5:16] <- res_tmp[index_min[1],5:16]}}
sum(is.na(GENES13$Pvalue))
# [1] 2
ind <- which(is.na(GENES13$Pvalue))
GENES13[ind,1:7]
#       Gene    Start      End N_Var                   Group max_MAF Pvalue
# 205  SOX21 95361878 95364799     2 missense;lof;synonymous     0.5     NA
# 222 GPR183 99946793 99959653     3                    <NA>      NA     NA
GENES13 <- subset(GENES13, subset=(!is.na(Pvalue)))
# WE HAVE 279 GENES ON CHR13
# Region of interest on CHR13: 10,904,797 - 40,904,797
GENES13_reg <- subset(GENES13, subset=((Start >= 10904797)&(End <=40904797)))
rownames(GENES13_reg) <- 1:nrow(GENES13_reg)
# 94/279 in region
SNPs13_link <- subset(SNPs13, subset=(Gene %in% GENES13_reg$Gene))
# 755/2,122 in region

GENES13_reg$Mito <- ifelse(GENES13_reg$Gene %in% mitoGenes$V1,"Yes","No")
sum(GENES13_reg$Mito == "Yes")
# [1] 4
#----------------------------------------------------------------------------------------------------------------------
ggplot(GENES13_reg,aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  geom_hline(yintercept=c(-log10(0.05/4),-log10(0.05/94)), color = c("orange", "midnightblue"), linetype = c("dashed", "solid"), lwd = 1) +
  scale_color_manual(values=c("midnightblue", "orange")) +
  geom_label_repel(data = subset(GENES13_reg, subset=(as.numeric(Pvalue) < 0.05|Mito == "Yes")), size = 4, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  theme_minimal() +
  labs(title="SAIGE-GENE: Encephalopathy Associations - SKAT-O Pvalues", x ="Chromosome 13 position", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 13), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
      panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.title.y = element_markdown(), axis.text.x = element_text(angle = 0, size = 11))
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                             Chr11:93692022:A/G
#-------------------------------------------------------------------------------------------------------------------------------
res_ALL001_DS_chr11 <- subset(res_ALL001_DS_ALL, subset=(CHR == 11))
res_ALL001_DS_chr11_reg <- subset(res_ALL001_DS_chr11, subset=((POS >= 83831017)&(POS <=129080721)))
sig <- (nrow(res_ALL001_DS_ALL) * 5e-08) / nrow(res_ALL001_DS_chr11_reg)  # 2.9e-06
ggplot(res_ALL001_DS_chr11_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=paste(MarkerID,Allele1,Allele2,sep=":"))) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("#fe7f2d", "#233d4d"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  #geom_label_repel(data = subset(res_ALL001_DS_chr11_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 30, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("#233d4d", "#233d4d"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(83831017,129080721,by=10000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
  theme_minimal() +
  labs(title="Chromosome 11: Encephalopathy Associations", x ="", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))
#----------------------------------------------------------------------------------------------------------------------
# How many genes were tested?
length(unique(res_GB_DS$SAIGE_SBA_DS_chr11$Region))
# 1152 GENES
SNPs11 <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(SNPs11) <- c("Variant","Annotation","Gene")
gf_chr11$N_Var <- ""
gf_chr11$Variants <- ""
for(i in 1:nrow(gf_chr11)){
  if(gf_chr11[i,2] == "var"){
    empty <- sum(gf_chr11[i,] == "") - 2
    nVar <- ncol(gf_chr11) - empty - 4
    gene <- gf_chr11[i,1]
    gf_chr11$N_Var[i] <- nVar
    gf_chr11$Variants[i] <- paste(gf_chr11[i,3:(nVar+2)], collapse="|")
    #print(i)
    SNPs_tmp <- data.frame(
      Variant = c(t(gf_chr11[i,3:(nVar+2)])),
      Annotation = c(t(gf_chr11[i+1,3:(nVar+2)])),
      Gene = rep(gene, times=nVar))
    SNPs11 <- rbind(SNPs11, SNPs_tmp)}
  else{
    gf_chr11$N_Var[i] <- ""
    gf_chr11$Variants[i] <- ""}
}
SNPs11$SNP_POS <- sapply(strsplit(SNPs11$Variant,":"), `[`, 2)
SNPs11$Gene <- sapply(strsplit(SNPs11$Gene,","), `[`, 1)
SNPs11 <- SNPs11[,c(1,4,2:3)]
for(i in 1:nrow(SNPs11)){
  SNPs11$Strand[i] <- with(anno_genes, Strand[GeneName == SNPs11$Gene[i]])
  SNPs11$Start[i] <-with(anno_genes, Start[GeneName == SNPs11$Gene[i]])
  SNPs11$End[i] <- with(anno_genes, End[GeneName == SNPs11$Gene[i]])}    

res_chr11 <- res_GB_DS$SAIGE_SBA_DS_chr11
res_chr11$Region <- sapply(strsplit(res_chr11$Region,","), `[`, 1)
gf_chr7$V1 <- sapply(strsplit(gf_chr11$V1,","), `[`, 1)
for(i in 1:nrow(res_chr11)){   
    res_chr11$Start[i] <-with(anno_genes, Start[GeneName == res_chr11$Region[i]])
    res_chr11$End[i] <- with(anno_genes, End[GeneName == res_chr11$Region[i]])
    res_chr11$N_Var[i] <- with(gf_chr11, N_Var[(V1 == res_chr11$Region[i])&(V2 == "var")])}
res_chr11 <- res_chr11[,c(1,14:16,2:13)]
res_chr11_reg <- subset(res_chr11, subset=((Start >= 83831017)&(End <=129080721)))
rownames(res_chr11_reg) <- 1:nrow(res_chr11_reg)

GENES11 <- data.frame(unique(res_chr11$Region))
colnames(GENES11) <- "Gene"
for(i in 1:nrow(GENES11)){
  GENES11$Start[i] <-with(anno_genes, Start[GeneName == GENES11$Gene[i]])
  GENES11$End[i] <- with(anno_genes, End[GeneName == GENES11$Gene[i]])
  GENES11$N_Var[i] <- with(gf_chr11, N_Var[(V1 == GENES11$Gene[i])&(V2 == "var")])}
# Get group with best pvalue
for(i in 1:nrow(GENES11)){
  print(i)
  genename <- GENES11[i,1]
  res_tmp <- subset(res_chr11, subset=(Region == genename))
  rownames(res_tmp) <- 1:nrow(res_tmp)
  if(nrow(res_tmp) == 1){GENES11[i,5:16] <- res_tmp[1,5:16]}
  else{
    minP <- min(res_tmp$Pvalue, na.rm=T) 
    index_min <- which(res_tmp$Pvalue == minP)
    if(length(index_min)>1){print(index_min)}
    GENES11[i,5:16] <- res_tmp[index_min[1],5:16]}}
sum(is.na(GENES11$Pvalue))
# [1] 2
ind <- which(is.na(GENES11$Pvalue))
GENES11[ind,1:7]
# Gene    Start      End N_Var                   Group max_MAF Pvalue
# 195  LMO1  8245855  8285425     2 missense;lof;synonymous     0.5     NA
# 567 UQCC3 62439144 62441158     2                    <NA>      NA     NA
GENES11 <- subset(GENES11, subset=(!is.na(Pvalue)))
# WE HAVE 1146 GENES ON CHR11
# Region of interest on CHR11: 83,831,017 - 129,080,721
GENES11_reg <- subset(GENES11, subset=((Start >= 83831017)&(End <=129080721)))
rownames(GENES11_reg) <- 1:nrow(GENES11_reg)
# 280/1146 in region
SNPs11_link <- subset(SNPs11, subset=(Gene %in% GENES11_reg$Gene))
# 1,922/8,034 in region

GENES11_reg$Mito <- ifelse(GENES11_reg$Gene %in% mitoGenes$V1,"Yes","No")
sum(GENES11_reg$Mito == "Yes")
# [1] 11
#----------------------------------------------------------------------------------------------------------------------
ggplot(GENES11_reg,aes(x=as.numeric(Start),y=-log10(as.numeric(Pvalue)), color = Mito, label=Gene)) + geom_point() +
  geom_hline(yintercept=c(-log10(0.05/11),-log10(0.05/280)), color = c("orange", "midnightblue"), linetype = c("dashed", "solid"), lwd = 1) +
  scale_color_manual(values=c("midnightblue", "orange")) +
  geom_label_repel(data = subset(GENES11_reg, subset=(as.numeric(Pvalue) < 0.05|Mito == "Yes")), size = 4, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  theme_minimal() +
  labs(title="SAIGE-GENE: Encephalopathy Associations - SKAT-O Pvalues", x ="Chromosome 11 position", y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 13), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), axis.title.y = element_markdown(), axis.text.x = element_text(angle = 0, size = 11))
