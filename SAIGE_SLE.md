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

plink --bfile ./input_data/genotypes_filtered --remove ./PCA/Non_EUR_IDs_PLINK_2.txt --make-bed --out ./input_data/genotypes_filtered_EUR

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

## SINGLE VARIANT ASSOCIATION TEST
### Step 1: Fitting the null logistic/linear mixed model [2 separate analyses with different MAF cut-off]
```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.01 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/SLE/sle_pheno.txt \
        --phenoCol=sle \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_MAF001 \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
```
Leave-one-chromosome-out is not applied
384  samples have genotypes
formula is  sle~age+het+sex 
247  samples have non-missing phenotypes
384  samples are in the sparse GRM
247  samples who have non-missing phenotypes are also in the sparse GRM
137  samples in geno file do not have phenotypes
247  samples will be used for analysis 

```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.05 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/SLE/sle_pheno.txt \
        --phenoCol=sle \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_MAF005 \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
```

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
	--GMMATmodelFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_MAF001.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_MAF001.varianceRatio.txt \
        --is_output_moreDetails=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=FALSE \
	--SAIGEOutputFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_SVA_sparseGRM_Firth_SLE_RESULTS_maf001.txt

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
	--GMMATmodelFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_MAF005.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_MAF005.varianceRatio.txt \
        --is_output_moreDetails=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=FALSE \
	--SAIGEOutputFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_SVA_sparseGRM_Firth_SLE_RESULTS_maf005.txt
```
---------------------
## SET-BASED TEST
### Step 1: Fitting the null logistic/linear mixed model
If a a sparse GRM is used for fitting the null model and variance ratios were estimated with sparse and null GRMs, in Step 2, the sparse GRM (–sparseGRMFile, –sparseGRMSampleIDFile) and variance ratios (–varianceRatioFile) are used as input.
```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/SLE/sle_pheno.txt \
        --phenoCol=sle \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_NO_MAF \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
```
Leave-one-chromosome-out is not applied
384  samples have genotypes
formula is  sle~age+het+sex 
247  samples have non-missing phenotypes
384  samples are in the sparse GRM
247  samples who have non-missing phenotypes are also in the sparse GRM
137  samples in geno file do not have phenotypes
247  samples will be used for analysis 


### Create VCF files per chromosome and their .csi index files (to use in step 2)
```
for i in {1..23}; do plink --bfile ./output_data/Geno_Imputed/merged_imputed_FINAL2 --chr "$i" --recode vcf-iid bgz --out ./output_data/Geno_Imputed/chr"$i"; done

for i in {1..23}; do tabix --csi -p vcf ./output_data/Geno_Imputed/chr$i.vcf.gz; done
```

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
	--GMMATmodelFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_NO_MAF.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_NullLMM_SLE_NO_MAF.varianceRatio.txt \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --r.corr=0 \
        --groupFile=./mydata/output_data/Geno_Imputed/Group_File_chr$i \
        --annotation_in_groupTest="lof,missense:lof,missense:lof:synonymous" \
        --maxMAF_in_groupTest=0.5 \
        --is_output_markerList_in_groupTest=TRUE \
        --LOCO=FALSE \
        --SAIGEOutputFile=./mydata/output_data/SAIGE/SLE/Results/SAIGE_SBA_chr$i
done
```