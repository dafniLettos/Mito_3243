# Association Analyses Using SAIGE

```
git clone git@github.com:weizhouUMICH/SAIGE.git
```
```
docker pull wzhou88/saige:1.3.0
```
## To mount data in the container
```
docker container run --rm -it -v .:/app/mydata 'wzhou88/saige:1.3.0' /bin/bash
```
Interactive session...
-------------------
### CREATING A SPARSE Genomic Relationship Matrix (GRM):
#### Use hard-called filtered genotypes, after remoning individuals of non European ancestry
```
Starting with: 408 Individuals and 558,579 Markers 

plink --bfile ./input_data/genotypes_filtered --remove ./PCA/Non_EUR_IDs_PLINK_3.txt --make-bed --out ./input_data/genotypes_filtered_EUR

558,579 variants and 384 people pass filters and QC
```
#### Filter out variants with MAF<0.05
```
plink --bfile ./input_data/genotypes_filtered_EUR --make-bed --maf 0.05 --out ./input_data/genotypes_filtered_EUR_maf0.05

348,330 variants and 384 people pass filters and QC
```
#### Prune variants for LD
```
plink --bfile ./input_data/genotypes_filtered_EUR_maf0.05 --indep-pairwise 50 5 0.3 --out ./input_data/LD_markers

182,351 of 348,330 variants removed (165,979 remaining)

plink --bfile ./input_data/genotypes_filtered_EUR_maf0.05 --extract ./input_data/LD_markers.prune.in --make-bed --out ./input_data/genotypes_filtered_EUR_maf0.05_LD
```
#### Run SAIGE to create sparse Genetic Relationship Matrix (GRM)
```
Rscript ../usr/local/bin/createSparseGRM.R \
	--plinkFile=./mydata/input_data/genotypes_filtered_EUR_maf0.05_LD \
	--nThreads=4  \
	--outputPrefix=./mydata/output_data/SAIGE/sparseGRM384	\
	--numRandomMarkerforSparseKin=5000	\
	--relatednessCutoff=0.05
```
Note: If the sample size is relatively small, e.g. N < 1000, may consider using sparse GRM...
sparseGRM--numRandomMarkerforSparseKin=5000--relatednessCutoff=0.05_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx
---------------------

## SINGLE VARIANT ASSOCIATION TEST
### Step 1: Fitting the null logistic/linear mixed model [2 separate analyses with different MAF cut-off]
#### Rename sparse GRM files
```
mv ./output_data/SAIGE/sparseGRM384--numRandomMarkerforSparseKin=5000--relatednessCutoff=0.05_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx ./output_data/SAIGE/sparseGRM384.mtx

mv ./output_data/SAIGE/sparseGRM384--numRandomMarkerforSparseKin=5000--relatednessCutoff=0.05_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt ./output_data/SAIGE/sparseGRM384_SampleIDs.txt
```
#### Fit the NULL model using two MAF cut-offs
```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.01 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/Encephalopathy/enceph_pheno.txt \
        --phenoCol=encephalopathy \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_MAF001 \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE

Leave-one-chromosome-out is not applied
384  samples have genotypes
formula is  encephalopathy~age+het+sex 
247  samples have non-missing phenotypes
384  samples are in the sparse GRM
247  samples who have non-missing phenotypes are also in the sparse GRM
137  samples in geno file do not have phenotypes
247  samples will be used for analysis

"isCovariateOffset=TRUE, so fixed effects coefficients won't be estimated."   
```
```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.05 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/Encephalopathy/enceph_pheno.txt \
        --phenoCol=encephalopathy \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_MAF005 \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
```
Note: --LOCO Whether to apply the leave-one-chromosome-out (LOCO) approach when fitting the null model using the full GRM [default=TRUE].
        LOCO should be TRUE if using all chromosomes to fit the NULL model (?!)

### Step 2: Performing single-variant association tests (accounting for relatedness)
```
Rscript ../usr/local/bin/step2_SPAtests.R \
        --bedFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.01.bed \
        --bimFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.01.bim \
        --famFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.01.fam \
        --AlleleOrder=alt-first \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
	--minMAF=0.01 \
	--GMMATmodelFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_MAF001.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_MAF001.varianceRatio.txt \
        --is_output_moreDetails=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=FALSE \
	--SAIGEOutputFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_SVA_sparseGRM_Firth_Enceph_RESULTS_maf001.txt

```
```
Rscript ../usr/local/bin/step2_SPAtests.R \
        --bedFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.05.bed \
        --bimFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.05.bim \
        --famFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.05.fam \
        --AlleleOrder=alt-first \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
	--minMAF=0.05 \
	--GMMATmodelFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_MAF005.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_MAF005.varianceRatio.txt \
        --is_output_moreDetails=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=FALSE \
	--SAIGEOutputFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_SVA_sparseGRM_Firth_Enceph_RESULTS_maf005.txt
```
Note: For binary traits, use –is_output_moreDetails=TRUE to output heterozygous and homozygous counts as well as allele frequencies in cases and controls
      --is_Firth_beta=TRUE and --pCutoffforFirth=0.05. The effect sizes of markers with p-value <= pCutoffforFirth will be estimated through the Firth’s Bias-Reduced Logistic Regression.
Firth's bias-reduced logistic regression in SAIGE (or any other software) should be used when dealing with complete or quasi-complete separation in logistic regression models, especially with small datasets or rare events, or when there is concern about bias in standard logistic regression estimates. It is also useful for analyses involving low-count variants, especially in balanced or unbalanced case-control studies PMID: 38298875, 39108101.
      When the sparse GRM was used for fitting the null model, set –LOCO=FALSE in Step 2
      --minMAF: Minimum minor allele frequency for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0]
      --minMAC: Minimum minor allele count for markers to be tested. The higher threshold between minMAC and minMAF will be used [default=0.5]
---------------------
## SET-BASED TEST
### Step 1: Fitting the null logistic/linear mixed model
If a a sparse GRM is used for fitting the null model and variance ratios were estimated with sparse and null GRMs, in Step 2, the sparse GRM (–sparseGRMFile, –sparseGRMSampleIDFile) and variance ratios (–varianceRatioFile) are used as input.
Use –is_output_markerList_in_groupTest=TRUE to output the markers used for each test.
```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/Encephalopathy/enceph_pheno.txt \
        --phenoCol=encephalopathy \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_NO_MAF \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE

Leave-one-chromosome-out is not applied
384  samples have genotypes
formula is  encephalopathy~age+het+sex 
247  samples have non-missing phenotypes
384  samples are in the sparse GRM
247  samples who have non-missing phenotypes are also in the sparse GRM
137  samples in geno file do not have phenotypes
247  samples will be used for analysis 
```

### Create VCF files per chromosome and their .csi index files (to use in step 2)
```
for i in {1..23}; do plink --bfile ./output_data/Geno_Imputed/merged_imputed_FINAL2 --chr "$i" --recode vcf-iid bgz --out ./output_data/Geno_Imputed/chr"$i"; done

for i in {1..23}; do tabix --csi -p vcf ./output_data/Geno_Imputed/chr$i.vcf.gz; done
```

### Check Variant Duplicates -- Group File Issue!!!!!
Run one chromosome at a time and take note of variants not found in VCFs [404 excel]
1. chr10: No markers in region NRBF2 are found in the VCF file
Search in Group_File_chr10 for NRBF2 variants: 10:64911911:A:C	10:64911911:A:G
zgrep '64911911' ./output_data/Geno_Imputed/chr10.vcf.gz | awk '{print $3,$4,$5}'
10:64911911 A C
10:64911911 A G
grep '10:64911911' ./output_data/QC_postImp/QC_postImp.frq | awk '{print $2,$3,$4,$5}'
10:64911911 C A 0.001225
10:64911911 G A 0.01348
Remove the one with the lowest MAF from group file.

2. chr19:No markers in region FEM1A are found in the VCF file
Search in Group_File_chr19 for FEM1A variants: 19:4792049:G:C	19:4792049:G:A
zgrep '4792049' ./output_data/Geno_Imputed/chr19.vcf.gz | awk '{print $3,$4,$5}'
19:4792049 G C
19:4792049 G A
grep '19:4792049' ./output_data/QC_postImp/QC_postImp.frq | awk '{print $2,$3,$4,$5}'
19:4792049 C G 0.001225
19:4792049 A G 0.003676
Remove the one with the lowest MAF from group file (annotation also).

### Step 2:  Perform region- or gene-based association tests

```
for i in {1..22}
do
Rscript ../usr/local/bin/step2_SPAtests.R \
        --vcfFile=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz \
        --vcfFileIndex=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz.csi \
        --vcfField=GT \
        --chrom=$i \
        --AlleleOrder=alt-first \
        --minMAF=0 \
        --minMAC=0.5 \
	--GMMATmodelFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_NO_MAF.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_NO_MAF.varianceRatio.txt \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --r.corr=0 \
        --groupFile=./mydata/output_data/Geno_Imputed/Group_File_chr$i \
        --annotation_in_groupTest="lof,missense:lof,missense:lof:synonymous" \
        --maxMAF_in_groupTest=0.5 \
        --is_output_markerList_in_groupTest=TRUE \
        --LOCO=FALSE \
        --SAIGEOutputFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_SBA_chr$i
done

--r.corr=R.CORR
                If r.corr = 1, only Burden tests will be performed. If r.corr = 0, SKAT-O tests will be performed and results for Burden tests and SKAT tests will be output too. [default = 0]
--is_output_markerList_in_groupTest=IS_OUTPUT_MARKERLIST_IN_GROUPTEST
                Whether to output the marker lists included in the set-based tests for each mask.[default=TRUE]
--maxMAF_in_groupTest=MAXMAF_IN_GROUPTEST
                Max MAF for markers tested in group test seperated by comma. e.g. 0.0001,0.001,0.01, [default=0.0001,0.001,0.01]
--maxMAC_in_groupTest=MAXMAC_IN_GROUPTEST
                Max MAC for markers tested in group test seperated by comma.The list will be combined with maxMAF_in_groupTest. e.g. 1,2 . By default, 0 and no maxMAC cutoff are applied. [default=0]
--groupFile=GROUPFILE
                Path to the file containing the group information for gene-based tests. Each gene/set has 2 or 3 lines in the group file. The first element is the gene/set name. The second element in the first line is to indicate whether this line contains variant IDs (var), annotations (anno), or weights (weight). The line for weights is optional. If not specified, the default weights will be generated based on beta(MAF, 1, 25). Use --weights.beta to change the parameters for the Beta distribution. The variant ids must be in the format chr:pos_ref/alt. Elements are seperated by tab or space.        

--annotation_in_groupTest=ANNOTATION_IN_GROUPTEST
                annotations of markers to be tested in the set-based tests seperated by comma. using ; to combine multiple annotations in the same test, e.g. lof,missense;lof,missense;lof;synonymous will test lof variants only, missense+lof variants, and missense+lof+synonymous variants. default: lof,missense;lof,missense;lof;synonymous
##      CHROM in VCF is chr1, --LOCO=FALSE so no --chrom is needed.
##      --maxMAF_in_groupTest=0.01,0.1,0.25 \  ##ROISIN
```

Burden tests assume that all causal variants have the same direction of effect (either all increase or all decrease the trait), while SKAT allows for different directions of effect
So, if Pvalue_Burden is not significanr (but Pvalue and Pvalue_SKAT are) then the variant should have a protective effect for the trait????

In SAIGE (and SAIGE-GENE+), the Cauchy combination method is not used for all regions because it requires a set of p-values, and some p-values might be equal to 1. In such cases, the Cauchy combination method fails, and a minimum p-value approach is used instead. Additionally, computational limitations can restrict the number of regions tested, especially when incorporating functional annotations and multiple MAF cutoffs. 
ALSO
For some genes (regions), no Cauchy is calculated because variants associated with the gene (from the group file) fall in a single category (missense|lof|synonymous).


Rscript ../usr/local/bin/step2_SPAtests.R \
        --vcfFile=./mydata/output_data/Geno_Imputed/chr19.vcf.gz \
        --vcfFileIndex=./mydata/output_data/Geno_Imputed/chr19.vcf.gz.csi \
        --vcfField=GT \
        --chrom=19 \
        --AlleleOrder=alt-first \
        --minMAF=0 \
        --minMAC=0.5 \
	--GMMATmodelFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_NO_MAF.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_NullLMM_Enceph_NO_MAF.varianceRatio.txt \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --r.corr=0 \
        --groupFile=./mydata/output_data/Geno_Imputed/Group_File_chr19 \
        --annotation_in_groupTest="lof,missense:lof,missense:lof:synonymous" \
        --maxMAF_in_groupTest=0.5 \
        --is_output_markerList_in_groupTest=TRUE \
        --LOCO=FALSE \
        --SAIGEOutputFile=./mydata/output_data/SAIGE/Encephalopathy/Results/SAIGE_SBA_chr19
