
# Annotate variants with ANNOVAR for Gene Burden Analysis (SAIGE-GENE+)
Visit https://www.openbioinformatics.org/annovar/annovar_download_form.php to register and receive the compressed file via email.

Open link with wget or using mozilla firefox.
```
tar -xvzf annovar.latest.tar.gz
rm annovar.latest.tar.gz
```
## Download additional databases (just for reference)
```
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGeneWithVer ./annovar/humandb/
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cytoBand ./annovar/humandb/
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome ./annovar/humandb/
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp151 ./annovar/humandb/
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp47a ./annovar/humandb/
```
## Get filtered genotypes and create vcf to use for annotation 
```
plink --bfile ./output_data/Geno_Imputed/merged_imputed_FINAL --keep-allele-order --recode vcf --out ./output_data/Geno_Imputed/merged_imputed_anno
```
## Transform vcf to ANNOVAR format
```
perl ./annovar/convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 ./output_data/Geno_Imputed/merged_imputed_anno.vcf > ./output_data/Geno_Imputed/merged_imputed_anno.avinput

NOTICE: Finished reading 7649412 lines from VCF file
NOTICE: A total of 7649383 locus in VCF file passed QC threshold, representing 7649383 SNPs (5248677 transitions and 2400706 transversions) and 0 indels/substitutions
NOTICE: Finished writing allele frequencies based on 3120948264 SNP genotypes (2141460216 transitions and 979488048 transversions) and 0 indels/substitutions for 408 samples
```
## Create multi-anno CSV
```
cd annovar \
perl table_annovar.pl ../output_data/Geno_Imputed/merged_imputed_anno.avinput \
  humandb/ \
  -buildver hg19 \
  -out ../output_data/Geno_Imputed/merged_imputed \
  -remove \
  --protocol refGene \
  -operation gx \
  -xref ./humandb/hg19_refGene.txt \
  -nastring . \
  -polish \
  -csvout \

cd ..  
```
```
# If you want to annotate with >1 databases
cd annovar
perl table_annovar.pl ../output_data/Geno_Imputed/merged_imputed_anno.avinput \
  humandb/ \
  -buildver hg19 \
  -out ../output_data/Geno_Imputed/merged_imputed \
  -remove \
  --protocol refGene,cytoBand,gnomad211_exome,avsnp151,dbnsfp47a \
  -operation gx,r,f,f,f \
  -xref ./humandb/hg19_refGene.txt \
  -nastring . \
  -polish \
  -csvout

cd ..  
```

## Create GROUP FILE per chromosome
```
*SHOULD MAKE IT WORK VIA TERMINAL*

Scripts/CreateGroupFile_ANNOVAR.R

# Rscript ./Scripts/GroupFile_for_SaigeGene_with_ANNOVAR.r \
# --inputfile=./output_data/Geno_Imputed/merged_imputed.hg19_multianno.csv \
# --outputfile=./output_data/Geno_Imputed/Group_File_ALL.txt
```
