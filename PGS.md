# Polygenic Score (PGS) Calculation Workflow

```
plink2 --vcf ./output_data/Geno_Imputed/merged_imputed_FINAL2.vcf 'dosage=DS' --not-chr X --make-pgen --out ./output_data/Geno_Imputed/merged_imputed_FINAL2
# To make into pgen (primary binary genotype file format; excluding chromosome X)

plink2 -pfile ./output_data/Geno_Imputed/merged_imputed_FINAL2 --rm-dup 'exclude-all' --make-pgen --out ./output_data/Geno_Imputed/merged_imputed_FINAL2_noDups
# To remove duplicates
# 14,394,148 variants loaded
# 17,662 duplicated IDs - 35,348 variants removed

plink2 -pfile ./output_data/Geno_Imputed/merged_imputed_FINAL2_noDups --export vcf vcf-dosage="DS" id-delim=/ --out ./output_data/PGS/Raw_data
14,358,800 variants loaded
```