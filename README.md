## Genetic and environmental modulators of mitochondrial disease in carriers of the m.3243A>G mitochondrial variant

This repository holds instruction steps and scripts used to:

1. Pre-process and impute genotypes (nuclear variants of m.3243 carriers) (GenoImp.md; ATCG_QC_Report.R)
2. Post-imputation QC and ancestry check (GenoImp.md; Post_Imp_QC.R)
3. Variant annotation using ANNOVAR for subsequent gene burden testing (Annovar.md; CreateGroupFile_ANNOVAR.R)
4. R scripts to prepare covariate files for association analyses
5. Steps to perform association analyses using SAIGE-GENE+ for:
	* Encephalopathy
	* Stroke-like Episodes (SLE)
	* Diabetes
6. R scripts to handle results from association analyses
7. Steps to perform Polygenic Score (PGS) analyses (PGS.md)


There is also a Dockerfile to create a linux-based image including necessary tools for these analyses.

An attempt is being made to incorporate some steps into a Nextflow workflow.


