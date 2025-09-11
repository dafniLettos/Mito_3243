# Running a GWAS for Diabetes to be used in the PGS Calculation

```
docker container run --rm -it -v .:/app/mydata 'wzhou88/saige:1.3.0' /bin/bash
```
## SINGLE VARIANT ASSOCIATION TEST
### Fit the NULL model using two MAF cut-offs [MAF>0.01 and MAF>0.05]
```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.01 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/Diabetes/diabetes_pheno.txt \
        --phenoCol=diabetes \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_MAF001 \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
```

```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2_maf0.05 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/Diabetes/diabetes_pheno.txt \
        --phenoCol=diabetes \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_MAF005 \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
```

### Perform single-variant association tests
```
for i in {1..22}
do
Rscript ../usr/local/bin/step2_SPAtests.R \
        --vcfFile=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz \
        --vcfFileIndex=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz.csi \
        --vcfField=DS \
        --chrom=$i \
        --AlleleOrder=alt-first \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
	--minMAF=0.01 \
	--GMMATmodelFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_MAF001.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_MAF001.varianceRatio.txt \
        --is_output_moreDetails=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=FALSE \
	--SAIGEOutputFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_SVA_maf001_chr$i.txt
done
```
```
for i in {1..22}
do
Rscript ../usr/local/bin/step2_SPAtests.R \
        --vcfFile=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz \
        --vcfFileIndex=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz.csi \
        --vcfField=DS \
        --chrom=$i \
        --AlleleOrder=alt-first \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
	--minMAF=0.05 \
	--GMMATmodelFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_MAF005.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_MAF005.varianceRatio.txt \
        --is_output_moreDetails=TRUE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.05 \
        --LOCO=FALSE \
	--SAIGEOutputFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_SVA_maf005_chr$i.txt
done
```
---------------------
## SET-BASED TEST
### Fit the null logistic/linear mixed model [no MAF cut-off]
```
Rscript ../usr/local/bin/step1_fitNULLGLMM.R \
        --plinkFile=./mydata/output_data/Geno_Imputed/merged_imputed_FINAL2 \
        --useSparseGRMtoFitNULL=TRUE \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --phenoFile=./mydata/output_data/SAIGE/Diabetes/diabetes_pheno.txt \
        --phenoCol=diabetes \
        --covarColList=age,het,sex \
        --sampleIDColinphenoFile=pid \
        --traitType=binary \
        --outputPrefix=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_NO_MAF \
        --nThreads=4 \
        --LOCO=TRUE \
        --IsOverwriteVarianceRatioFile=TRUE
```
### Perform gene-based association tests
```
for i in {1..22}
do
Rscript ../usr/local/bin/step2_SPAtests.R \
        --vcfFile=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz \
        --vcfFileIndex=./mydata/output_data/Geno_Imputed/chr$i.vcf.gz.csi \
        --vcfField=DS \
        --chrom=$i \
        --AlleleOrder=alt-first \
        --minMAF=0 \
        --minMAC=0.5 \
	--GMMATmodelFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_NO_MAF.rda \
	--varianceRatioFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_NullLMM_Diabetes_NO_MAF.varianceRatio.txt \
        --sparseGRMFile=./mydata/output_data/SAIGE/sparseGRM384.mtx \
        --sparseGRMSampleIDFile=./mydata/output_data/SAIGE/sparseGRM384_SampleIDs.txt \
        --r.corr=0 \
        --groupFile=./mydata/output_data/Geno_Imputed/Group_File_chr$i \
        --annotation_in_groupTest="lof,missense:lof,missense:lof:synonymous" \
        --maxMAF_in_groupTest=0.5 \
        --is_output_markerList_in_groupTest=TRUE \
        --LOCO=FALSE \
        --SAIGEOutputFile=./mydata/output_data/SAIGE/Diabetes/Results/SAIGE_SBA_DS_chr$i
done
```