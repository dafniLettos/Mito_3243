setwd("/home/c0091724/Desktop/PostDoc/Mito_3243/")
library(optparse)
library(data.table)
library(methods)

missenseList = c("nonsynonymous SNV", "nonsynonymous")
synonymousList = c("synonymous SNV", "synonymous")
lofList = c("frameshift deletion", "frameshift insertion",
            "nonframeshift deletion", "nonframeshift insertion",
            "splicing", "stopgain", "stoploss", "startloss")
combineList = c(missenseList, synonymousList, lofList)
combineList = unique(combineList)

# MY IMPUTED GENOTYPES + FILTERS
Data1 = data.table::fread("./output_data/Geno_Imputed/merged_imputed_FINAL2.hg19_multianno.csv",sep=",")
# 14,391,036 markers (ds)
# 14,706,768 markers


# SARAH'S IMPUTED GENOTYPES + FILTERS
#Data1 = data.table::fread("./output_data/Encephalopathy_Sarah_Rep/merged_imputed.hg19_multianno.csv",sep=",")
## 7,770,220 markers

# SARAH'S IMPUTED GENOTYPES WITHOUT FILTERS
#Data1 = data.table::fread("./output_data/Encephalopathy_Sarah_Rep/merged_imputed_NoFilters.hg19_multianno.csv",sep=",")
# 15,669,496 markers

###colnames(Data1)[4:5] <- c("Alt","Ref")
#Data1$Gene.refGene <- gsub("\"", "", Data1$Gene.refGene)
#Data1$ExonicFunc.refGene <- gsub("\"", "", Data1$ExonicFunc.refGene)
#Data1$Func.refGene <- gsub("\"", "", Data1$Func.refGene)
Data1$Other = paste0(Data1$Chr,":", Data1$Start, ":", Data1$Ref, ":", Data1$Alt)
Data1 = Data1[,c("Gene.refGene", "ExonicFunc.refGene", "Other", "Func.refGene", "Start", "Chr")]
Data1$Chr = paste0("chr",Data1$Chr)
Data1 = Data1[which((Data1$ExonicFunc.refGene %in% combineList) | (Data1$Func.refGene  %in% combineList )), , drop=F]
# 119,695 [ds]
# [52,931 rows OLD]  119,726 rows
# 53,883 from Sarah's
# 127,644 from Sarah's (no filters)
rowswithDupGeneNames =  grep(";", Data1$Gene.refGene)
# 45 [ds]
# [22 OLD] 45
# 23 from Sarah's
# 48 from Sarah's (no filters)

#if(length(rowswithDupGeneNames) > 0){
  Data2 = Data1[-rowswithDupGeneNames,]  
  for(i in rowswithDupGeneNames){
    print(i)
    a = Data1[i,, drop=F]
    newgeneList = strsplit(a$Gene.refGene, split=";")[[1]]
    newgeneList = unique(newgeneList)
    aDF = NULL
    for(j in newgeneList){
      a$Gene.refGene = j
      aDF = rbind(aDF, a)}
  Data2 = rbind(Data2, aDF)	
  }
#  }else{Data2 = Data1}

Data2 = Data2[order(Data2$Chr,Data2$Start), ]
Data2$Group = rep("missense", nrow(Data2))
Data2$Group[which((Data2$ExonicFunc.refGene %in% lofList) | (Data2$Func.refGene  %in% lofList ))] = "lof"
Data2$Group[which(Data2$ExonicFunc.refGene %in% missenseList)] = "missense"
Data2$Group[which(Data2$ExonicFunc.refGene %in% synonymousList)] = "synonymous"
geneList = unique(Data2$Gene.refGene)
# 16,873 [ds]
# [14,082 genes OLD] 16,875
# 14,169 genes from Sarah's
# 17,047 genes from Sarah's (no filters)

write.table(Data2, file="./output_data/Geno_Imputed/GroupFile_ALL_temp.txt", row.names=F, col.names=T, quote=F, sep="\t")
# Data2 = data.table::fread("./output_data/Geno_Imputed/GroupFile_ALL_temp.txt", sep="\t")
# 119,730 rows

#write.table(Data2, file="./output_data/Encephalopathy_Sarah_Rep/GroupFiles/GroupFile_ALL_temp.txt", row.names=F, col.names=T, quote=F, sep="\t")
#Data2 = data.table::fread("./output_data/Encephalopathy_Sarah_Rep/GroupFiles/GroupFile_ALL_temp.txt", sep="\t")

#write.table(Data2, file="./output_data/Encephalopathy_Sarah_Rep/GroupFiles_NoFilters/GroupFile_ALL_temp.txt", row.names=F, col.names=T, quote=F, sep="\t")
#Data2 = data.table::fread("./output_data/Encephalopathy_Sarah_Rep/GroupFiles_NoFilters/GroupFile_ALL_temp.txt", sep="\t")


# Create a list of data frames corresponding to the equivalent chromosome subset of Data2
DATA <- list()
for(i in 1:22){DATA[[i]] <- subset(Data2, subset=(Chr==paste0("chr",i)))}

for(i in 1:22){
for(gene in unique(DATA[[i]]$Gene.refGene)){
    rowHeaderData = cbind(rep(gene,2), c("var", "anno"))
    geneTemp = DATA[[i]][which(DATA[[i]]$Gene.refGene == gene), c("Start", "Other", "Group")]
    geneTemp$rank = 1
    geneTemp$rank[which(geneTemp$Group == "missense")] = 2
    geneTemp$rank[which(geneTemp$Group == "synonymous")] = 3
    geneTemp$Start = as.numeric(geneTemp$Start)
    geneTemp = geneTemp[order(geneTemp$rank),]
    geneTemp = geneTemp[!duplicated(geneTemp$Other), ]
    geneTemp2 = geneTemp[order(geneTemp$Start),]
    geneTempt = t(geneTemp2[,2:3])
    geneTempt_new = cbind(rowHeaderData, geneTempt)
    write.table(geneTempt_new, paste0("./output_data/Geno_Imputed/Group_File_chr",i), sep="\t", col.names=F, row.names=F, quote=F, append=T)
    #write.table(geneTempt_new, paste0("./output_data/Encephalopathy_Sarah_Rep/GroupFiles/Group_File_chr",i), sep="\t", col.names=F, row.names=F, quote=F, append=T)
    #write.table(geneTempt_new, paste0("./output_data/Encephalopathy_Sarah_Rep/GroupFiles_NoFilters/Group_File_chr",i), sep="\t", col.names=F, row.names=F, quote=F, append=T)
}}



