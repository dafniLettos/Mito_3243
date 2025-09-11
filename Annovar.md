
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
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb cytoBand ./annovar/humandb/
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_exome ./annovar/humandb/
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp151 ./annovar/humandb/
perl ./annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp47a ./annovar/humandb/
```
## Transform vcf to ANNOVAR format
```
bgzip -dc ./output_data/Geno_Imputed/merged_imputed_FINAL2.vcf.gz > ./output_data/Geno_Imputed/merged_imputed_FINAL2.vcf

perl ./annovar/convert2annovar.pl -includeinfo -allsample -withfreq -format vcf4 ./output_data/Geno_Imputed/merged_imputed_FINAL2.vcf > ./output_data/Geno_Imputed/merged_imputed_FINAL2.avinput

merged_imputed_FINAL2:
NOTICE: Finished reading 14391080 lines from VCF file
NOTICE: A total of 14391036 locus in VCF file passed QC threshold, representing 14391036 SNPs (9914380 transitions and 4476656 transversions) and 0 indels/substitutions
NOTICE: Finished writing allele frequencies based on 5526157824 SNP genotypes (3807121920 transitions and 1719035904 transversions) and 0 indels/substitutions for 384 samples

rm ./output_data/Geno_Imputed/merged_imputed_FINAL2.vcf
```

## Create multi-anno CSV
```
cd annovar
perl table_annovar.pl ../output_data/Geno_Imputed/merged_imputed_FINAL2.avinput \
  humandb/ \
  -buildver hg19 \
  -out ../output_data/Geno_Imputed/merged_imputed_FINAL2 \
  -remove \
  --protocol refGene \
  -operation gx \
  -xref ./humandb/hg19_refGene.txt \
  -nastring . \
  -polish \
  -csvout

cd .. 

rm rm ./output_data/Geno_Imputed/merged_imputed_FINAL2.avinput
```
### Annotate with >1 database (not necessary)
```
cd annovar
perl table_annovar.pl ../output_data/Geno_Imputed/merged_imputed_FINAL2.avinput \
  humandb/ \
  -buildver hg19 \
  -out ../output_data/Geno_Imputed/merged_imputed_FINAL2 \
  -remove \
  --protocol refGene,cytoBand,gnomad211_exome,avsnp151 \
  -operation gx,r,f,f \
  -xref ./humandb/hg19_refGene.txt \
  -nastring . \
  -polish \
  -csvout

cd .. 

rm rm ./output_data/Geno_Imputed/merged_imputed_FINAL2.avinput
```

## Create GROUP FILE per chromosome
```
*SHOULD MAKE IT WORK VIA TERMINAL*

Scripts/CreateGroupFile_ANNOVAR.R

# Rscript ./Scripts/CreateGroupFile_ANNOVAR.R \
# --inputfile=./output_data/Geno_Imputed/merged_imputed_FINAL2.hg19_multianno.csv \
# --outputfile=./output_data/Geno_Imputed/Group_File_ALL.txt
```
