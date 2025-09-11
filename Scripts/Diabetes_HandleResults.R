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

# Read-in SAIGE Results (GWAS)
res_files_res_ALL001_ds <- paste0(rep("SAIGE_SVA_maf001_chr",22),c(1:22),rep(".txt",22))
res_SVA001_DS <- list()
for(file in res_files_res_ALL001_ds){
  print(paste0("./output_data/SAIGE/Diabetes/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/Diabetes/Results/",file), header=T)
  res_SVA001_DS[[file]] <- data}

res_files_res_ALL005_ds <- paste0(rep("SAIGE_SVA_maf005_chr",22),c(1:22),rep(".txt",22))
res_SVA005_DS <- list()
for(file in res_files_res_ALL005_ds){
  print(paste0("./output_data/SAIGE/Diabetes/Results/",file))
  data <- read.table(file=paste0("./output_data/SAIGE/Diabetes/Results/",file), header=T)
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
#-------------------------------------------------------------------
#-------------------------------------------------------------------
#                     GENOME-WIDE TEST Manhattan Plots 
#-------------------------------------------------------------------
#-------------------------------------------------------------------
# Merge results per chromosome into single dataframe [MAF 0.01]
res_ALL001 <- bind_rows(res_SVA001_DS)
# 7,586,961 SNPs
write.table(res_ALL001, file="./output_data/SAIGE/Diabetes/Results/SAIGE_SVA_DS_ALL_MAF001.txt", sep="\t", col.names=T, row.names=F, quote=F)
hist(res_ALL001$p.value)

# Merge results per chromosome into single dataframe
res_ALL005 <- bind_rows(res_SVA005_DS)
# 5,464,655 SNPs
write.table(res_ALL005, file="./output_data/SAIGE/Diabetes/Results/SAIGE_SVA_DS_ALL_MAF005.txt", sep="\t", col.names=T, row.names=F, quote=F)
hist(res_ALL005$p.value)

data_ALL001 <- res_ALL001 %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(POS)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(res_ALL001, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, POS) %>%
  mutate(POScum=POS+tot)
axisdf = data_ALL001 %>%
  group_by(CHR) %>%
  summarize(center=(max(POScum) + min(POScum))/2)
ggplot(data_ALL001, aes(x = POScum, y = -log10(p.value), label=paste(MarkerID, Allele1, Allele2,sep=":"))) +
  geom_hline(yintercept = -log10(5e-08), color = "#fe7f2d", linetype = "solid", lwd = 1) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#233d4d", "darkslategray4"), 23)) +
  geom_label_repel(data = subset(data_ALL001, subset=(-log10(p.value) >= 10)), size = 3, max.overlaps = 20, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  theme_minimal() +
  labs(title="SAIGE SVA: Diabetes Associations",
       x ="Chromosome",
       y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))

data_ALL005 <- res_ALL005 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(POS)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(res_ALL001, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, POS) %>%
  mutate(POScum=POS+tot)
ggplot(data_ALL005, aes(x = POScum, y = -log10(p.value), label=paste(MarkerID, Allele1, Allele2,sep=":"))) +
  geom_hline(yintercept = -log10(5e-08), color = "#fe7f2d", linetype = "solid", lwd = 1) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#233d4d", "darkslategray4"), 23)) +
  geom_label_repel(data = subset(data_ALL005, subset=(-log10(p.value) >= 10)), size = 3, max.overlaps = 20, box.padding = 0.01, point.padding = 0.01, segment.color = 'grey50') +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  theme_minimal() +
  labs(title="SAIGE SVA: Diabetes Associations",
       x ="Chromosome",
       y ="-log<sub>10</sub>(p)") +
  theme(title = element_text( size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_markdown(),
        axis.text.x = element_text(angle = 0, size = 11))
#-----------------------------------------------------------------
#                           QQ Plots
#-----------------------------------------------------------------
qqplot_test(res_ALL001$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) Genome-Wide")
qqplot_test(res_ALL005$p.value) + theme_minimal() +
  labs(title="Q-Q Plot for single variant analysis (SVA) Genome-Wide")
#---------------------------------------------------------------------
#               Genomic inflation factor (lambda)
#---------------------------------------------------------------------
qchisq(1-median(res_ALL001$p.value),1)/qchisq(0.5,1)
# [1] 1.022573
qchisq(1-median(res_ALL005$p.value),1)/qchisq(0.5,1)
# [1] 1.017232