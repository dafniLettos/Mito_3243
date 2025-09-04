# GenoIMpPGS WORKFLOW (sans Nextflow)

## Steps Summary
0. Filter out samples with >10% missing call rates --prefix:genotypes
1. Generate MAF, HWE P-values and Missingness using PLINK
2. Get ATCG SNP list (R script)
3. Apply filters (PLINK): keep markers with P_hwe > 0.00000001, MAF > 0.01, Missingness < 0.03 and exclude ATCG (ambiguous) markers
4. Convert to a single compressed vcf -- genotypes_filtered.vcf.gz
5. Check alignment of markers to the positive strand of GRChb37
6. Pre-imputation data preparation (Will Rayner's Toolbox)
7. Michigan Imputation Server (online/interactive)
8. Post-imputation quality control

Convert from PLINK to VCF format (per chromosome)
Impute genotypes using the Michigan Imputation Server
PCA on imputed genotypes using PLINK (?)
Calculate Polygenic Scores (pgsc_calc)

## Filter out samples with >50% missing call rates ("dummy" samples used for LD analysis)
#### Outside project directory:
```
plink --bfile merged_cleaner --mind 0.5 --make-bed --out genotypes
```
* -- mind filters out all samples with missing call rates > 50% (default 0.1 - 10%)
* Individuals: 803 - 395 = 408
* Markers: 643,848

## 1.Generate QC Files
```
plink --bfile ./input_data/genotypes --freq --hardy --missing --het --out ./output_data/QC/QC
```
* Number of markers that don't pass HWE criterion: 0
* Number of markers with missing call rate > 0.3: 0
* Number of markers with MAF < 0.01: 48,697

## 2. Get ATCG SNPs (ambiguous markers)
```
Rscript ./Scripts/ATCG_QC_Report.r
```
* Number of ambiguous (ATCG) markers: 40,258
* Generate ATCG_SNPs.txt file

## 3. Apply Filters
```
plink --bfile ./input_data/genotypes --exclude ./output_data/QC/ATCG_SNPs.txt --make-bed --hwe 0.00000001 --maf 0.01 --geno 0.3 --out ./input_data/genotypes_filtered
```
* Individuals: 408
* Markers: 558,579

## 4. Convert to a single compressed vcf
```
plink --bfile ./input_data/genotypes_filtered --keep-allele-order --not-chr XY --recode vcf bgz --out ./input_data/genotypes_filtered
```
------------
------------
------------
## 5 (optional). Align markers to the positive strand of Genome Reference Consortium Human genome build 37 co-ordinates (GRChb37)
Get reference fasta from GENCODE and correct chr sequence names:
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz`
bgzip -d GRCh37.p13.genome.fa.gz
```
Note: In this FASTA file each chromosome is a sequence (obviously) but the sequence names are in the format > chr1 1. This causes an issue when parsing the file, given that our chromosome names are numeric (hg). To fix this, I edit the FASTA file as follows:
```
for i in {1..22} X; do sed -i -e 's/>chr\([$i]\) />/g' GRCh37.p13.genome.fa; done
```
Run the stats to learn the number of REF allele mismatches and the number of non-biallelic sites:
```
bcftools +fixref ./input_data/genotypes_filtered.vcf.gz -- -f GRCh37.p13.genome.fa
```
------------
------------
------------
## 6. Pre-imputation data preparation
Follow M.I.S instructions
https://imputationserver.readthedocs.io/en/latest/prepare-your-data/

Download  the latest version of the Will Rayner's Toolbox tool from
https://www.chg.ox.ac.uk/~wrayner/tools/

Here we used v4.2.13. All scripts and files are in "Rayners_Toolbox directory.
A frequency file is necessary for the tool's PERL script to run.
```
plink --freq --bfile ./input_data/genotypes_filtered --out ./output_data/QC/Post_filtering

gunzip ./Rayners_Toolbox/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

perl ./Rayners_Toolbox/HRC-1000G-check-bim.pl -b ./input_data/genotypes_filtered.bim -f ./output_data/QC/Post_filtering.frq -r ./Rayners_Toolbox/HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

sh ./input_data/Run-plink.sh
```
This PERL script created files:
* Chromosome-genotypes_filtered-HRC.txt
* Exclude-genotypes_filtered-HRC.txt
* Force-Allele1-genotypes_filtered-HRC.txt
* FreqPlot-genotypes_filtered-HRC.txt
* ID-genotypes_filtered-HRC.txt
* LOG-genotypes_filtered-HRC.txt
* Position-genotypes_filtered-HRC.txt
* Strand-Flip-genotypes_filtered-HRC.txt
* LOG-genotypes_filtered-HRC.txt
* Run-plink.sh

This bash script runs plink multiple times, using files generated above AND splits into chromosomes BUT ignores chrX [excludes all markers on chrX]. In r1 HRC there are only autosomes, so 1-22 only considered by script. HRC r1.1 does include data for chromosome X and still they are getting excluded.
To deal with the chrX issue I subset my genotype file (using plink) to only include chrX. Then, using awk I subset the reference to only include chrX sites as well as the Post_filtering.afreq (for the MAFs) and re-run the PERL script using the flag -c 23:

```
plink --bfile ./input_data/genotypes_filtered --chr X --make-bed --out ./input_data/chrX_filtered

awk '$1 == "chrX" || $1 == "X" {print}' ./Rayners_Toolbox/HRC.r1-1.GRCh37.wgs.mac5.sites.tab > ./Rayners_Toolbox/HRC.r1-1.GRCh37.wgs.mac5.sites_chrX.tab

awk '$1 == "chrX" || $1 == "X" {print}' ./output_data/QC/Post_filtering.frq > ./output_data/QC/Post_filtering_chrX.frq

perl ./Rayners_Toolbox/HRC-1000G-check-bim.pl -b ./input_data/chrX.bim -f ./output_data/QC/Post_filtering_chrX.frq -r ./Rayners_Toolbox/HRC.r1-1.GRCh37.wgs.mac5.sites_chrX.tab -h -c 23

mv ./input_data/Run-plink.sh ./input_data/Run-plink_chrX.sh
sh ./input_data/Run-plink_chrX.sh
```
### Organise files in "input_data" directory
```
mkdir ./input_data/Filter_files
mv ./input_data/*txt ./input_data/Filter_files
mkdir ./input_data/MIS_input
mv ./input_data/genotypes_filtered-updated-chr*.vcf ./input_data/MIS_input
mv ./input_data/chrX_filtered-updated-chr23.vcf ./input_data/MIS_input
tar czf ./input_data/genotypes_filtered-updated_PLINK.tar.gz ./input_data/genotypes_filtered-updated-chr* ./input_data/chrX_filtered-updated-chr23.*
rm ./input_data/genotypes_filtered-updated-chr*
rm ./input_data/chrX_filtered-updated-chr23.*
```
### Compress VCFs to upload on the Michigan Imputation Server 2
```
for i in {1..23}; do bgzip -c ./input_data/MIS_input/genotypes_filtered-updated-chr"$i".vcf > ./input_data/MIS_input/chr"$i".vcf.gz; done

rm ./input_data/MIS_input/genotypes_filtered-updated-chr*
```

## 7. Michigan Imputation Server

* Job:	PGS3243_CHR1-22
* Pipeline Version	v2.0.6
* Date	Thu Jun 12 09:44:25 EDT 2025
* Samples:	408
* Chromosomes:	1 10 11 12 13 14 15 16 17 18 19 2 20 21 22 3 4 5 6 7 8 9
* SNPs:	378,805
* Chunks:       154
* Datatype:	unphased
* Build:	hg19
* Reference Panel:	hrc-r1.1 (hg19)
* rsq Filter: 0.3
* Population:	eur
* Phasing:	beagle
* Mode:	QC & imputation
* Alternative allele frequency > 0.5 sites:	27,763
* Reference Overlap:	100.00 %
* Match:	378,805
* Password: [YSha#3WfUyE6j
-------
* Job:	PGS3243_CHRX
* Pipeline Version	v2.0.6
* Date	Thu Jun 12 10:38:28 EDT 2025
* Samples:	408
* Chromosomes:	23
* SNPs:	7,113
* Chunks:       8
* Datatype:	unphased
* Build:	hg19
* Reference Panel:	hrc-r1.1 (hg19)
* rsq Filter: 0.3
* Population:	eur
* Phasing:	beagle
* Mode:	QC & imputation
* Alternative allele frequency > 0.5 sites:	631
* Reference Overlap:	100.00 %
* Match:	7,113
* Password: J>Cb5fL9MoBvk(
-------
550862
### Download & Unzip Files
```
curl -sL https://imputationserver.sph.umich.edu/get/AZzL1TBzj0nz4CjSLq0p4c9wXYelk4dSF5du5If2 | bash

wget https://imputationserver.sph.umich.edu/share/results/5a6c0be43cae00582735c44dfdc674e73477c6c627bed6ab6ab1d4ce79fcd561/chr_X.zip
wget https://imputationserver.sph.umich.edu/share/results/886a3f1f9bc42fcedea278109640387df50c919a41ac201acb7b7b373b7b211d/qc_report.txt
wget https://imputationserver.sph.umich.edu/share/results/db02f7e6479e393b7cb57ca90085e0f52f73e625e02cfb674bc585b0b90e45d0/quality-control.html

for i in {1..22} X; do unzip ./output_data/Geno_Imputed/chr_"$i".zip -d ./output_data/Geno_Imputed/; done

It will ask for the password for every chromosome...
```
### Merge into single vcf file
```
cd ./output_data/Geno_Imputed/

bcftools concat -Oz -o merged_imputed.vcf.gz chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz chrX.dose.vcf.gz
```
### Compress files that will not be used now
```
tar czf ./output_data/Geno_Imputed/EmpiricalDose_and_info_tar.gz ./output_data/Geno_Imputed/chr*.empiricalDose.* ./output_data/Geno_Imputed/chr*.info.gz

rm ./output_data/Geno_Imputed/chr*.empiricalDose.*
rm ./output_data/Geno_Imputed/chr*.info.gz
```

## 8. Post-imputation QC
#### Generate PLINK files from merged VCF
```
plink --vcf ./output_data/Geno_Imputed/merged_imputed.vcf.gz --make-bed --double-id --out ./output_data/Geno_Imputed/merged_imputed

14,844,986 variants
408 people (0 males, 0 females, 408 ambiguous)
14,844,986 variants and 408 people pass filters and QC
```
#### Generate QC Files
```
plink --bfile ./output_data/Geno_Imputed/merged_imputed --freq --hardy --missing --het --ibc --out ./output_data/QC_postImp/QC_postImp
```
--------------------------
### I. Variant-Level QC
--------------------------
We have already applied the Rsq > 0.3 filter. This is the "imputation quality score" calculated by minimac which is the software that the Michigan Imputation Server uses for the imputation step.
In addition, there were no inconsistencies with the reference panel to consider, as we had a 100% overlap as seen in the MIS QC report.

#### Apply HWE filter of 0.00000001
```
plink --bfile ./output_data/Geno_Imputed/merged_imputed --hwe 0.00000001 --make-bed --out ./output_data/Geno_Imputed/merged_imputed_hwe

--hwe: 138,218 variants removed due to Hardy-Weinberg exact test.
14,706,768 variants and 408 people pass filters and QC.
```
--------------------------
### II. Sample-Level QC
--------------------------
#### Heterozygosity (??)
```
Post_Imp_QC.r
```
#### Ancestry: PCA on hard-called genotypes using 1000 Genomes [23 individuals to be removed]
##### Recode families in our dataset (don't know why)
```
plink --bfile ./input_data/genotypes_filtered --make-bed --update-ids ./input_data/recode_fams.txt -out ./input_data/genotypes_filtered_recodeFams
```
##### Replace Affy SNP ID with snpID of format CHR:POS:REF:ALT and extract our SNP list!!!!!!!!!!!!!!!!!!!!!!!!!!!
```
plink --bfile ./input_data/genotypes_filtered_recodeFams --update-map ./input_data/affyID_to_snpID.txt --update-name --make-bed --out ./input_data/genotypes_filtered_recodeFams_snpID
plink --bfile ./input_data/genotypes_filtered_recodeFams_snpID --write-snplist --make-bed --out ./input_data/genotypes_filtered_recodeFams_snpID
```
##### Extract SNP from the 1000G genotypes that are in our list where MAF > 0.05 and Call rate ≥ 95% 
```
plink --bfile ./PCA/1000G --extract ./input_data/genotypes_filtered_recodeFams_snpID.snplist --maf 0.05 --geno 0.05 --write-snplist --make-bed --out ./PCA/1000G_SNPsub
```
##### Keep from our genotypes the overlapping SNPs and place genotype files in the PCA directory
```
plink --bfile ./input_data/genotypes_filtered_recodeFams_snpID --extract ./PCA/1000G_SNPsub.snplist --make-bed --out ./PCA/M3243
```
##### Merge the two datasets [1000G and M3243] 
```
plink --bfile ./PCA/1000G_SNPsub --bmerge ./PCA/M3243 --make-bed --out ./PCA/1000G_M3243
13,288 variants and 2,912 people pass filters and QC
```
##### Run PCA on the merged data
```
plink --bfile ./PCA/1000G_M3243 --pca --out ./PCA/1000G_M3243
```
RRRRRRRRRRRRRRRRRRRRRRR script
which will also produce a txt file with the IDs to be excluded [format: famId_pID\tfamId_pID; Non_EUR_IDs.txt]
For a reason, famID==pID==famId_pID. PLINK needs 2 columns in the --remove file, so I just export the same column twice (maybe fix later)
Also, Ro says 23 individuals excluded in the paper (?).
##### Exclude non-EUR individuals from genotype file 
```
plink --bfile ./output_data/Geno_Imputed/merged_imputed_hwe --remove ./PCA/Non_EUR_IDs_PLINK.txt --make-bed --out ./output_data/Geno_Imputed/merged_imputed_FINAL

14,706,768 variants and 384 people pass filters and QC
```
##### Fix Pid and Famid on genotype file
```
plink --bfile ./output_data/Geno_Imputed/merged_imputed_FINAL --make-bed --update-ids ./input_data/IDs_updateIDs.txt -out ./output_data/Geno_Imputed/merged_imputed_FINAL2
```
Note: The 'vcf', 'vcf-fid', and 'vcf-iid' modifiers result in production of a VCFv4.2 file. 'vcf-fid' and 'vcf-iid' cause family IDs and within-family IDs respectively to be used for the sample IDs in the last header row, while 'vcf' merges both IDs and puts an underscore between them (in this case, a warning will be given if an ID already contains an underscore).
If the 'bgz' modifier is added, the VCF file is block-gzipped. (Gzipping of other --recode output files is not currently supported.)
The A2 allele is saved as the reference and normally flagged as not based on a real reference genome ('PR' INFO field value). When it is important for reference alleles to be correct, you'll usually also want to include --a2-allele and --real-ref-alleles in your command.

-------------------------------------------------------------------------
### III. Create extra files with MAF filters to be used where necessary
-------------------------------------------------------------------------
```
plink --bfile ./output_data/Geno_Imputed/merged_imputed_FINAL2 --make-bed --maf 0.01 --out ./output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.01

7,708,725 (out of 14,706,768) variants and 384 people pass filters and QC
```
```
plink --bfile ./output_data/Geno_Imputed/merged_imputed_FINAL2 --make-bed --maf 0.05 --out ./output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.05

5,411,913 (out of 14,706,768) variants and 384 people pass filters and QC
```
---------------------------------------------------------------------------------------------------------------------
### Extract variants present in chr7q22 and chr5,6 significant|suggestive linkage regions for encephalopathy and SLE
---------------------------------------------------------------------------------------------------------------------
```
bcftools view -r 7:91468550-121468550 ./output_data/Geno_Imputed/chr7.vcf.gz > ./output_data/Geno_Imputed/CHR7q22_VARIANTS_ALL.vcf
140,722 (wc -l) - 10 (header) = 140,712 variants in 7q22 without any MAF filter

bcftools view -r 6:143863601-173863601 ./output_data/Geno_Imputed/chr6.vcf.gz > ./output_data/Geno_Imputed/CHR6reg_VARIANTS_ALL.vcf
159,176 (wc -l) - 10 (header) = 159,166 variants in chr6 region without any MAF filter

bcftools view -r 5:140749194-188279842 ./output_data/Geno_Imputed/chr5.vcf.gz > ./output_data/Geno_Imputed/CHR5reg_VARIANTS_ALL.vcf
215,185 (wc -l) - 10 (header) = 215,175 variants in chr5 region without any MAF filter
```