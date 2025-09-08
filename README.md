# Genetic and environmental modulators of mitochondrial diesease in carriers of the m.3243A>G mitochondrial variant

This repository holds instruction steps and scripts used to:

1. Pre-process and impute genotypes (nuclear variants of m.3243 carriers) [GenoImp.md]
2. Post-imputation QC and ancestry check [GenoImp.md]
3. Variant annotation using ANNOVAR for subsequent gene burden testing [Annovar.md]
4. Steps to perform association analyses using SAIGE-GENE+ for:
	* Encephalopathy
	* Stroke-like Episodes (SLE)
	* Diabetes

There is also a Dockerfile to create a linux-based image including necessary tools for these analyses.

An attempt is being made to incorporate dome steps into a Nextflow workflow.


