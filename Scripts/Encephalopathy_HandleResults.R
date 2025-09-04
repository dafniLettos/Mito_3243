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

# Read-in SAIGE-GENE Results
#res_chr7 <- read.table(file="./output_data/SAIGE/Encephalopathy/Results/SAIGE_SBA_chr7", header=T)
res_files <- list.files(path="./output_data/SAIGE/Encephalopathy/Results/", pattern="SAIGE_SBA_chr*")
res_files <- paste0(rep("SAIGE_SBA_chr",22),c(1:22))
res_GB <- list()
for(file in res_files){
  print(paste0("./output_data/SAIGE/Encephalopathy/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/Encephalopathy/Results/",file), header=T)
  res_GB[[file]] <- data}

# Read-in SAIGE Results (GWAS)
res_ALL001 <- fread(file="./output_data/SAIGE/Encephalopathy/Results/SAIGE_SVA_sparseGRM_Firth_Enceph_RESULTS_maf001.txt", header=T)
res_ALL005 <- fread(file="./output_data/SAIGE/Encephalopathy/Results/SAIGE_SVA_sparseGRM_Firth_Enceph_RESULTS_maf005.txt", header=T)
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
hist(res_ALL001$p.value)
res_ALL001_chr7 <- subset(res_ALL001, subset=(CHR == 7))
res_ALL001_chr7_reg <- subset(res_ALL001_chr7, subset=((POS >= 91468550)&(POS <=121468550)))
hist(res_ALL001_chr7_reg$p.value)
sum(res_ALL001_chr7_reg$p.value < 5e-08)
#[1] 0
res_ALL001_chr7_reg$p.value_log <- -log10(as.numeric(res_ALL001_chr7_reg$p.value))
res_ALL001_chr7_reg$SNP_ID <- paste(res_ALL001_chr7_reg$MarkerID,res_ALL001_chr7_reg$Allele1,res_ALL001_chr7_reg$Allele2,sep=":")
#--------------------------
hist(res_ALL005$p.value)
res_ALL005_chr7 <- subset(res_ALL005, subset=(CHR == 7))
res_ALL005_chr7_reg <- subset(res_ALL005_chr7, subset=((POS >= 91468550)&(POS <=121468550)))
hist(res_ALL005_chr7_reg$p.value)
sum(res_ALL005_chr7_reg$p.value < 5e-08)
#[1] 0
res_ALL005_chr7_reg$p.value_log <- -log10(as.numeric(res_ALL005_chr7_reg$p.value))
res_ALL005_chr7_reg$SNP_ID <- paste(res_ALL005_chr7_reg$MarkerID,res_ALL005_chr7_reg$Allele1,res_ALL005_chr7_reg$Allele2,sep=":")

#----------------------------------
#         Manhattan Plot
#----------------------------------
#manhattan_enceph_chr7 <- res_ALL_chr7_reg[,c(26,1,2,25)]
#colnames(manhattan_enceph_chr7) <- c("SNP", "CHR", "BP", "P")
#manhattan(manhattan_enceph_chr7, main = "SAIGE SVA: Encephalopathy Associations",
#           col = "midnightblue", xlim = c(91468550,121468550), ylim = c(0,8),
#          suggestiveline = F, genomewideline = F)
#abline(h=-log10(5.7e-06), col = "hotpink4", lty=2, lwd = 3)
#abline(h=-log10(5e-08), col = "darkolivegreen", lwd = 3)

## Suggestive significance level for the proportion of the genome tested:
sig <- (nrow(res_ALL001) * 5e-08) / nrow(res_ALL001_chr7_reg)  # 5.5e-06  -log10(sig) --- 5.3
ggplot(res_ALL001_chr7_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
   geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
   geom_point(alpha = 0.70) +
   geom_label_repel(data = subset(res_ALL001_chr7_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 2, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
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

sig <- (nrow(res_ALL005) * 5e-08) / nrow(res_ALL005_chr7_reg)  # 5.7e-06
ggplot(res_ALL005_chr7_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=SNP_ID)) +
 geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("skyblue4", "hotpink4"), linetype = c("dashed", "solid"), lwd = 1) +
 geom_point(alpha = 0.70) +
 geom_label_repel(data = subset(res_ALL005_chr7_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 2, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
 scale_color_manual(values = rep(c("midnightblue", "midnightblue"))) +
 scale_size_continuous(range = c(0.5, 3)) +
 theme_minimal() +
 labs(title="SAIGE SVA: Encephalopathy Associations", x ="Chromosome 7 position", y ="-log<sub>10</sub>(p)") +
 theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
       legend.position = "none",
       panel.grid.major.x = element_blank(),
       panel.grid.minor.x = element_blank(),
       axis.title.y = element_markdown(),
       axis.text.x = element_text(angle = 0, size = 11))

#-----------------------------------------------------------------
#            QQ Plots & Genomic Inflation Factor
#-----------------------------------------------------------------
qqplot_test(res_ALL001$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) Genome-Wide")
qqplot_test(res_ALL001_chr7$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) on chr7")
qqplot_test(res_ALL001_chr7_reg$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) on chr7 (7q22)")
#-----------------------------------------------------
qqplot_test(res_ALL005$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) Genome-Wide")
qqplot_test(res_ALL005_chr7$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) on chr7")
qqplot_test(res_ALL005_chr7_reg$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) on chr7 (7q22)")
#-------------------------------------------------------------------------------------
#           Calculate genomic inflation factor (lambda)
#-------------------------------------------------------------------------------------
# Lambda (λ) quantifies the deviation of observed test statistics from what's expected under the null hypothesis.
# A genomic inflation factor of 1 indicates no bias, while values greater than 1 suggest inflation,
# potentially due to factors like population structure or cryptic relatedness
# Convert p-values to z-scores: z_scores <- qnorm(Pvalue / 2, lower.tail = FALSE) 
# OR CALCULATE chi-squared stats STRAIGHT AWAY: chi_squared_stats <- z_scores^2
# Convert z-scores to chi-squared stats: chi_squared_stats <- qchisq(Pvalue, df = 1, lower.tail = FALSE)
# Calculate the genomic inflation factor (lambda): λ = OBSERVED MEDIAN/EXPECTED MEDIAN of chi-squared statistics
# lambda <- median(chi_squared_stats) / qchisq(0.5, df=1)
#-------------------------------------------------------------------------------------
qchisq(1-median(res_ALL001$p.value),1)/qchisq(0.5,1)
#[1] 1.00265
qchisq(1-median(res_ALL005$p.value),1)/qchisq(0.5,1)
# [1] 0.9948877
qchisq(1-median(res_ALL001_chr7$p.value),1)/qchisq(0.5,1)
#[1] 0.9995107
qchisq(1-median(res_ALL005_chr7$p.value),1)/qchisq(0.5,1)
# [1] 0.9911946
qchisq(1-median(res_ALL001_chr7_reg$p.value),1)/qchisq(0.5,1)
#[1] 1.065139
qchisq(1-median(res_ALL005_chr7_reg$p.value),1)/qchisq(0.5,1)
# [1] 0.988943
#----------------------------------------------------------
#----------------------------------------------------------
#                     GENE-BURDEN TEST
#----------------------------------------------------------
#----------------------------------------------------------
# How many genes were tested?
length(unique(res_GB$SAIGE_SBA_chr7$Region))
# 760 GENES
gf_chr7 <- read.table("./output_data/Geno_Imputed/Group_File_chr7", header=F, col.names = paste0("V",seq_len(max(count.fields("./output_data/Geno_Imputed/Group_File_chr7", sep = '\t')))), fill = TRUE)
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
res_chr7 <- res_GB$SAIGE_SBA_chr7
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
# [1] 0
# ind <- which(is.na(GENES$Pvalue))
# GENES[ind,1:7]
## Remove genes with NA pvalues
#GENES <- subset(GENES, subset=(!is.na(Pvalue)))

# WE HAVE 756 GENES ON CHR7
# Region of interest on CHR7: 91,468,550 - 121,468,550
GENES_reg <- subset(GENES, subset=((Start >= 91468550)&(End <=121468550)))
rownames(GENES_reg) <- 1:nrow(GENES_reg)
# 177/757 in region

# Subset of variants in that region
SNPs_7q22 <- subset(SNPs, subset=(Gene %in% GENES_reg$Gene))
# 1,205/5,690 in region
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#-----------------------------------------
#     Calculate GB lambda GENOME-WIDE
#-----------------------------------------
for(i in 1:22){
  res_GB[[i]]$Chr <- rep(i,nrow(res_GB[[i]]))
  res_GB[[i]] <- res_GB[[i]][,c(14,1:13)]}
# Merge results per chromosome into single dataframe
res_GB_ALL <- bind_rows(res_GB)
res_GB_ALL_noCauchy <- subset(res_GB_ALL, subset=(Group != "Cauchy"))
write.table(res_GB_ALL, file="./output_data/SAIGE/Encephalopathy/Results/SAIGE_SBA_ALL_CHR", sep="\t", col.names=T, row.names=F, quote=F)

# SKAT-O
qchisq(1-median(res_GB_ALL_noCauchy$Pvalue, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 0.9769728

# BURDEN
qchisq(1-median(res_GB_ALL_noCauchy$Pvalue_Burden, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 0.997195

# SKAT
qchisq(1-median(res_GB_ALL_noCauchy$Pvalue_SKAT, na.rm=TRUE),1)/qchisq(0.5,1)
# [1] 1.005488
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
sig <- (nrow(res_ALL001) * 5e-08) / nrow(res_ALL001_chr7_reg)  # 5.5e-06
# -log10(sig) --- 5.3
MANHATTAN7 <- ggplot(res_ALL001_chr7_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value))) +
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
  geom_hline(yintercept=-log10(0.05/177), color = "#fe7f2d", linetype = "solid", lwd = 1) +
  scale_color_manual(values=c("#233d4d", "skyblue3")) +
  scale_x_continuous(breaks=seq(10000000,150000000,by=5000000), labels=scales::comma) +
  #guides(col = "none") +
  geom_label_repel(data = subset(GENES_reg, subset=(-log10(Pvalue) > 2|Mito == "Yes")),  fill = alpha(c("white"),0.8), size = 3, max.overlaps = 100, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
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

min(res_ALL001_chr7_reg$p.value)
which(res_ALL001_chr7_reg$p.value == 5.737226e-06)
res_ALL001_chr7_reg[11615,]
# 7:97072414:G:T --- rs62497093

subset(GENES_reg, subset=(Gene == "PIK3CG"))# Pvalue = 0.001788563
subset(GENES_reg, subset=(Gene == "PLOD3"))# Pvalue = 0.007459452
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#                                        EXPORT RESULTS FOR SUPPLEMENTARY
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
res_ALL001_chr7_reg <- res_ALL001_chr7_reg[,c(1:3,25,4:24)]
colnames(res_ALL001_chr7_reg)[4] <- "CHR:POS:A1:A2"
write.table(res_ALL001_chr7_reg, file="./output_data/SAIGE/Encephalopathy/SupTable2_SAIGE_chr7q22.tab", sep="\t", row.names=F, col.names=T, quote=F)

GENES_reg[167,5] <- "missense;lof;synonymous"
GENES_reg[167,6] <- 0.5
write.table(GENES_reg, file="./output_data/SAIGE/Encephalopathy/SupTable3_SAIGEGENE_chr7q22.tab", sep="\t", row.names=F, col.names=T, quote=F)
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

res_ALL001_chr7_reg$GeneLabel <- "Other"
for(i in 1:nrow(res_ALL001_chr7_reg)){
  if(res_ALL001_chr7_reg$POS[i] %in% PIK3CG_var){res_ALL001_chr7_reg$GeneLabel[i] <- "PIK3CG"}
  else if(res_ALL001_chr7_reg$POS[i] %in% PLOD3_var){res_ALL001_chr7_reg$GeneLabel[i] <- "PLOD3"}}

# GWAS was done with MAF > 0.01 cut-off so some variants were removed
sum(res_ALL001_chr7_reg$GeneLabel=="PIK3CG")
#[1] 8
sum(res_ALL001_chr7_reg$GeneLabel=="PLOD3")
#[1] 2

sig <- (nrow(res_ALL001) * 5e-08) / nrow(res_ALL001_chr7_reg)  # 5.5e-06
ggplot(data = res_ALL001_chr7_reg, aes(x = POS, y = -log10(p.value))) +
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
#                     That is an actual hit, influencing gene NMNAT1 (nicotinamide mononucleotide adenylyltransferase 1)
#                                                 ~  NAD+ Biosynthesis!!!
#-------------------------------------------------------------------------------------------------------------------------------
res_ALL001_chr1 <- subset(res_ALL001, subset=(CHR == 1))
res_ALL001_chr1_reg <- subset(res_ALL001_chr1, subset=((POS >= 693731)&(POS <=20800000)))

sig <- (nrow(res_ALL001) * 5e-08) / nrow(res_ALL001_chr1_reg)  # 6.7e-06
ggplot(res_ALL001_chr1_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value), label=paste(MarkerID,Allele1,Allele2, sep=":"))) +
  geom_hline(yintercept = -log10(5e-08), color = "orange", linetype = "dashed", lwd = 1) +
  geom_point(alpha = 0.70) +
  geom_label_repel(data = subset(res_ALL001_chr1_reg, subset=(-log10(p.value) >= 6)), size = 3, max.overlaps = 10, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_color_manual(values = rep(c("#233d4d", "#233d4d"))) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_x_continuous(breaks=seq(693731,20800000,by=5000000), labels=scales::comma) +
  #scale_x_continuous(breaks=NULL) +
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
#                             Chr13: 1:10048151_G/T : Pval=4.69 × 10^-8 [rs139490379]
#                     That is an actual hit, influencing gene NMNAT1 (nicotinamide mononucleotide adenylyltransferase 1)
#                                                 ~  NAD+ Biosynthesis!!!
#-------------------------------------------------------------------------------------------------------------------------------
res_ALL001_chr13 <- subset(res_ALL001, subset=(CHR == 13))
res_ALL001_chr13_reg <- subset(res_ALL001_chr13, subset=((POS >= 10904797)&(POS <=40904797)))

sig <- (nrow(res_ALL001) * 5e-08) / nrow(res_ALL001_chr13_reg)  # 5.6e-06
ggplot(res_ALL001_chr13_reg, aes(x = POS, y = -log10(p.value), color = as.factor(CHR), size = -log10(p.value))) +
  geom_hline(yintercept = c(-log10(sig), -log10(5e-08)), color = c("#fe7f2d", "#233d4d"), linetype = c("dashed", "solid"), lwd = 1) +
  geom_point(alpha = 0.70) +
  #geom_label_repel(data = subset(res_ALL001_chr13_reg, subset=(-log10(p.value) >= 4)), size = 3, max.overlaps = 2, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
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

