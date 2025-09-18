#!/usr/bin/env Rscript

# (SOS) - IF I RUN THE SCRIPT THROUGH NEXTFLOW, NOT SURE WHAT THE wd WOULD BE [WE SHALL SEE]
setwd("/home/c0091724/Desktop/PostDoc/Mito_3243/")
list.files()
library("ggplot2")
library("ggrepel")
library("dplyr")

# Define NOT-IN function
`%!in%` <- Negate(`%in%`)

# MISSINGNESS VS HETEROZYGOSITY PLOT
imiss <- read.delim(file="./output_data/QC_postImp/QC_postImp.smiss", header=TRUE)
het <- read.delim(file="./output_data/QC_postImp/QC_postImp.het", header=TRUE)
QC_individuals <- cbind(het, imiss[5])
QC_individuals$HET = (QC_individuals$E.HOM. - QC_individuals$O.HOM.)/QC_individuals$E.HOM.
# Calculate Heterozygosity Mean & SD
het_mean <- mean(QC_individuals$HET, na.rm=TRUE)
het_sd = sd(QC_individuals$HET, na.rm=TRUE)
QC_individuals$label <- ifelse(QC_individuals$HET>=(het_mean+(3*het_sd))|QC_individuals$HET<=(het_mean-(3*het_sd))|QC_individuals$F>0.03,QC_individuals$IID,"")

pdf(file = "./output_data/QC_postImp/Heterozygosity_plot.pdf",  width = 8.3, height = 5.8)
ggplot(QC_individuals, aes(x=F_MISS, y=HET)) + theme_bw() +
  geom_point() + xlim(0,0.1) +
  xlab("Call Rate") +
  ylab("Heterozygosity") +
  geom_vline(xintercept=0.03, color="lightskyblue3") +
  #geom_hline(yintercept=het_mean, color="lightskyblue3") +
  geom_hline(yintercept=(het_mean+(3*het_sd)), linetype="dashed", color="maroon4") +
  geom_hline(yintercept=(het_mean-(3*het_sd)), linetype="dashed", color="maroon4") +
  geom_text_repel(label=QC_individuals$label, nudge_y=0.001, nudge_x=0.001, size=4)
dev.off()
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
# PRINCIPAL COMPONENT ANALYSIS [vs 1000 Genomes]

# Produce file that includes Affy SNP IDs and snpIDs [Format of CHR:POS:REF:ALL]
probe_info = read.table(file="./input_data/probe_info.txt",header=TRUE,sep="\t")
probe_info$snpID = paste(probe_info$Chromosome, probe_info$Chromosome.Start, probe_info$Allele.A, probe_info$Allele.B, sep=":")
write.table(probe_info[,c("affy_snp_id", "snpID")], file="./input_data/affyID_to_snpID.txt", row.names=F, col.names=F, quote=F)

# Eigenvectors of merged data (1000 Genomes and m.3243 cohort)
eigenvec <- read.table(file="./PCA/1000G_M3243.eigenvec",header=FALSE)
colnames(eigenvec) <- c("FamID","PID",1:20)

# 1000 Genomes PED file
ped_1kGenomes <- read.table(file="./PCA/1000G_20130606_g1k.ped",header=TRUE,sep="\t")

# Create an m.3243 PED file
m3243_PIDs <- subset(eigenvec, subset=(FamID!=0) , select=(PID))[,1]
ped_m3243 <- data.frame("Family.ID"=NA, "Individual.ID"= m3243_PIDs, "Paternal.ID" =NA, "Maternal.ID" =NA, "Gender" =NA, "Phenotype" =NA, "Population" ="M3243", "Relationship" =NA, "Siblings"  =NA,
                     "Second.Order" =NA,  "Third.Order"  =NA, "Other.Comments"=NA)

# Create a merged PED
ped_merged <- rbind(ped_1kGenomes, ped_m3243)

# Create PCA object (PIDs, Population, PCs)
PCA_full <- merge(ped_merged, eigenvec, by.x = "Individual.ID", by.y = "PID")
# Create an extra column that summarizes ancestry
for(i in 1:nrow(PCA_full)){
  if(PCA_full$Population[i] %in% c("ACB","ASW","ESN","GWD","LWK","MSL","YRI")){PCA_full$PopSum[i] <- "African"}
  else if(PCA_full$Population[i] %in% c("CLM","MXL","PEL","PUR")){PCA_full$PopSum[i] <- "Ad Mixed American"}
  else if(PCA_full$Population[i] %in% c("CDX","CHB","CHS","JPT","KHV")){PCA_full$PopSum[i] <- "East Asian"}
  else if(PCA_full$Population[i] %in% c("CEU","FIN","GBR","IBS","TSI")){PCA_full$PopSum[i] <- "European"}
  else if(PCA_full$Population[i] %in% c("BEB","GIH","ITU","PJL","STU")){PCA_full$PopSum[i] <- "South Asian"}
  else{PCA_full$PopSum[i] <- "Mito"}}

PCA_full <- PCA_full[,c(1,7,34,14:33)]
for(i in 1:20){colnames(PCA_full)[3+i] <- paste0("PC",i)}

write.table(PCA_full, file="./PCA/PCA_DATA_all.txt", sep="\t", col.names=T, row.names=F, quote=F)
PCA_full <- read.table(file="./PCA/PCA_DATA_all.txt",header=TRUE,sep="\t")

# PCA_dev <- as.data.frame(matrix(0, nrow=2912, ncol=20))
# # For each element calculate (PC1[i,j]-mean(PC1))^2
# for(j in 1:20){for(i in 1:nrow(PCA_dev)){PCA_dev[i,j] <- (PCA_full[i,3+j]-mean(PCA_full[,3+j]))^2}}
# # Sum the squared deviations for each PC
# PCA_sdev <- rep(0,20)
# for(i in 1:20){PCA_sdev[i] <- sum(PCA_dev[,i])}
# # Divide the sum by (N-1), where N is the number of elements in the vector,
# # to get the sample variance (for population variance, you would divide by N). 
# for(i in 1:20){PCA_sdev[i] <- PCA_sdev[i]/(nrow(PCA_dev)-1)}
# # Take the square root
# for(i in 1:20){PCA_sdev[i] <- PCA_sdev[i]^2}
 
#eigenvalues <- pca$sdev^2
eigenvalues <- read.table(file="./PCA/1000G_M3243.eigenval",header=FALSE)

pdf(file = "./PCA/PCA_VarianceExplained.pdf",  width = 8.3, height = 5.8)
#barplot((pca$sdev^2 / sum(pca$sdev^2))[1:20], 
barplot(eigenvalues$V1 / sum(eigenvalues$V1), 
        main = "Scree Plot - Variance Explained by Genotype PCs",
        xlab = "PCs",
        ylab = "Proportion of Variance",
        names.arg = 1:20,
        xlim = c(0,23),
        ylim = c(0,0.5),
        col = "#233d4d")
dev.off()
--------------------------------
#   PCA Plots (without labels)
--------------------------------
pca_colours <- c("violet", "goldenrod2", "darkolivegreen", "slategray2", "firebrick", "yellowgreen")

ggplot(data = PCA_full) + theme_bw() +
  geom_point(aes(x = PC1, y = PC2, colour = PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC1, y = PC2, colour = PopSum)) +
  labs(x = "Principal Component 1", y = "Principal Component 2", title = "PCA Plot: PCs 1 VS 2") 

ggplot(data = PCA_full) + theme_bw() +
  geom_point(aes(x = PC2, y = PC3, colour = PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC1, y = PC2, colour = PopSum)) +
  labs(x = "Principal Component 2", y = "Principal Component 3", title = "PCA Plot: PCs 2 VS 3")

ggplot(data = PCA_full) + theme_bw() +
  geom_point(aes(x = PC3, y = PC4, colour = PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC3, y = PC4, colour = PopSum)) +
  labs(x = "Principal Component 3", y = "Principal Component 4", title = "PCA Plot: PCs 3 VS 4")

PCA_full$label <- ifelse((PCA_full$PopSum=="Mito" & PCA_full$PC1<0) | (PCA_full$PopSum=="Mito" & PCA_full$PC2>-0.0175), as.character(PCA_full$Individual.ID), "")
#PCA_full$label <- ifelse((PCA_full$PopSum=="Mito" & PCA_full$PC1<0) | (PCA_full$PopSum=="Mito" & PCA_full$PC2>-0.0175) | PCA_full$PopSum=="Mito" & PCA_full$PC3>0.004, as.character(PCA_full$Individual.ID), "")
sum(PCA_full$label!="")
#[1] 24 INDIVIDUALS FALLING OUTSIDE "EUR-BORDERS" (???)
--------------------------------
#   PCA Plots (with labels)
--------------------------------
pdf(file="./PCA/PCA_1_2.pdf")
ggplot(data = PCA_full, aes(x = PC1, y = PC2, label=label)) + theme_bw() +
  geom_point(aes(colour=PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC1, y = PC2, colour = PopSum)) +
  labs(x = "Principal Component 1", y = "Principal Component 2", title = "PCA Plot: PCs 1 VS 2") +
  geom_text_repel(data=subset(PCA_full, subset=(label!="")),max.overlaps = 20)
dev.off()

pdf(file="./PCA/PCA_1_2_noIDs.pdf")
ggplot(data = PCA_full, aes(x = PC1, y = PC2, label=label)) + theme_bw() +
  geom_point(aes(colour=PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC1, y = PC2, colour = PopSum)) +
  labs(x = "Principal Component 1", y = "Principal Component 2", title = "PCA Plot: PCs 1 VS 2")
dev.off()

pdf(file="./PCA/PCA_2_3.pdf")
ggplot(data = PCA_full, aes(x = PC2, y = PC3, label=label)) + theme_bw() +
  geom_point(aes(colour=PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC2, y = PC3, colour = PopSum)) +
  labs(x = "Principal Component 2", y = "Principal Component 3", title = "PCA Plot: PCs 2 VS 3") +
  geom_text_repel(data=subset(PCA_full, subset=(label!="")),max.overlaps = 20)
dev.off()

pdf(file="./PCA/PCA_2_3_noIDs.pdf")
ggplot(data = PCA_full, aes(x = PC2, y = PC3, label=label)) + theme_bw() +
  geom_point(aes(colour=PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC2, y = PC3, colour = PopSum)) +
  labs(x = "Principal Component 2", y = "Principal Component 3", title = "PCA Plot: PCs 2 VS 3")
dev.off()

pdf(file="./PCA/PCA_3_4.pdf")
ggplot(data = PCA_full, aes(x = PC3, y = PC3, label=label)) + theme_bw() +
  geom_point(aes(colour=PopSum)) +
  scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
  geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC1, y = PC2, colour = PopSum)) +
  labs(x = "Principal Component 3", y = "Principal Component 4", title = "PCA Plot: PCs 3 VS 4") +
  geom_text_repel(data=subset(PCA_full, subset=(label!="")),max.overlaps = 20)
dev.off()


plotA <-ggplot(eigenvalues, aes(x=1:20, y=(V1/sum(V1))))+
        geom_bar(stat="identity", fill="#233d4d")+
        scale_x_continuous(breaks=seq(1,20,by=1)) +
        theme_bw() +
        xlab("Principal Components") +
        ylab("Proportion of Variance") +
        theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          axis.title.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
          ggtitle("A. % Variance Explained by Principal Components")
  
plotB <-ggplot(data = PCA_full, aes(x = PC1, y = PC2, label=label)) + theme_bw() +
        geom_point(aes(colour=PopSum)) +
        scale_colour_manual(values=pca_colours, name="Population",
                      breaks = c("African", "European", "East Asian", "Ad Mixed American", "Mito", "South Asian"), 
                      labels = c("African", "European", "East Asian", "Ad Mixed American", "m3243A>G", "South Asian")) +
        geom_point(data = (subset(PCA_full, Population == "Mito")), aes(x = PC1, y = PC2, colour = PopSum)) +
        xlab("Principal Component 1") +
        ylab("Principal Component 2") +
        theme(axis.text.x = element_text(angle = 0, size = 12, vjust = 0.5, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.ticks = element_line(colour = "black", linewidth = 0.2),
          axis.title.y = element_text(size = 12, colour = "black"),
          axis.title.x = element_text(size = 12, colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
          ggtitle("B.PCA Plot")

grid.arrange(plotA, plotB, nrow=2)
#--------------------------------------------------------------------------------------------
# CREATE FILE TO BE USED IN PLINK TO --remove THESE SAMPLES FROM THE IMPUTED GENOTYPE FILES
#--------------------------------------------------------------------------------------------
# Read-in the COVARIATES [Fam, Age, Sex] & BINARY NMDAS
# We need the famID
COVARIATES <- read.table(file="./input_data/FullCohortBinaryPedigree_170821.tab", header=TRUE)

for(i in 1:nrow(PCA_full)){
  if(PCA_full$Individual.ID[i] %in% COVARIATES$pid){
  PCA_full$Family.ID[i] <- with(COVARIATES, famid[pid == PCA_full$Individual.ID[i]])}
  else{PCA_full$Family.ID[i] <- NA}
}
non_eur_ids <- subset(PCA_full, subset=(label!=""), select=c(Family.ID,Individual.ID))
rownames(non_eur_ids) <- 1:nrow(non_eur_ids)
non_eur_ids$Fam_PID <- paste(non_eur_ids$Family.ID,non_eur_ids$Individual.ID,sep="_")
non_eur_ids$Fam_PID2 <- non_eur_ids$Fam_PID
# 24 IDs to export
write.table(non_eur_ids[,3:4], file="./PCA/Non_EUR_IDs_PLINK.txt", row.names=F, col.names=F, quote=F)
write.table(non_eur_ids[,1:3], file="./PCA/Non_EUR_IDs_Full.txt", row.names=F, col.names=T, quote=F)

# NEED TO MAKE VCF PER CHROMOSOME, BUT HAVE ISSUE WITH IDs AND RECODING TO VCF
# I AM GOING TO USE PLINK --UPDATE-IDS ON ./output_data/Geno_Imputed/merged_imputed_FINAL BEFORE SPLITING INTO VCFs
# FOR THAT I NEED A FILE THAT CONTAINS 4 COLUMNS: Old FamID; Old pID; New FamID; New pID
# USE THE INFO FROM ID_hash...

# First, make a hash:
geno_fam <- read.table(file="./output_data/Geno_Imputed/merged_imputed_FINAL.fam", header=FALSE)
COVARIATES <- read.table(file="./input_data/FullCohortBinaryPedigree_170821.tab", header=TRUE)
COVARIATES <- COVARIATES[order(COVARIATES[,"famid"]),]
rownames(COVARIATES) <- 1:nrow(COVARIATES)
COVARIATES$ID <- paste(COVARIATES$famid, COVARIATES$pid, sep="_")
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
COVARIATES2 <- subset(COVARIATES, subset=(pid %!in% c(1:12)))
COVARIATES2 <- subset(COVARIATES2, subset=(Genotype=="Yes"))
ID_hash <- COVARIATES2[,c("famid","pid","ID")]
write.table(ID_hash, file="./input_data/FamID_pID_ID.txt", row.names=F, col.names=F, quote=F)

swap_ids <- data.frame(geno_fam[1:2])
swap_ids$V3 <- ID_hash$famid
swap_ids$V4 <- ID_hash$pid
write.table(swap_ids, file="./input_data/IDs_updateIDs.txt", row.names=F, col.names=F, quote=F)


