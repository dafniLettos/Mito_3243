#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */
// PLINK Genotype Files - Name (Not QCed)
params.input_genoS = "${projectDir}/input_data/genotypes"
//Input Directory
params.inputdir = "${projectDir}/input_data"
// Scripts Directory
params.scriptDir = "${projectDir}/Scripts"
// Output Directories
params.outdir = "${projectDir}/output_data"

//genofiles = Channel.fromFilePairs("${params.input_genoS}.{bed,bim,fam}",size:3).set {genofiles}
//def genofiles = listOfFiles("${params.inputdir}/genotypes.*")

/*
 * Generate Minor Allele Frequencies, Hardy–Weinberg Equilibrim P-value and Missingness Calculations
 */
process PLINK_QC {

    container 'genoimppgs'
    publishDir params.outdir, mode: 'symlink'

    input:
        path genofiles
        //tuple val(genofiles_name), files(pfiles) from genofiles
        //path pfiles from genofiles
        //tuple path(genofiles), file("${genofiles}.bed"), file("${genofiles}.bim"), file("${genofiles}.fam")
        //set file (genofiles) from genofiles
        //val genofilesgenotypes_name
        //set file (genofiles) from genofiles
        path outPrefix

    output:
        //path("${outPrefix}.*")// into plink_QC_files
        path "${outPrefix}.afreq"
        path "${outPrefix}.hardy"
        path "${outPrefix}.smiss"
        path "${outPrefix}.vmiss"
        path "${outPrefix}.het"
        path "${outPrefix}.log"

    script:
    //    outPrefix = 'QC'
    //    inPrefix =  "geno0types"
        inPrefix = genofiles.getBaseName()
    """
        plink2 --bfile "${params.inputdir}/${inPrefix}" --freq --hardy --missing --het --out "${outPrefix}"
    """
}

/*
 * Generate ATCG SNP list
 * ToDo --- HTML report
 *
process R_QC {

    container 'genoimppgs'
    publishDir "${params.outdir}/QC", mode: 'copy'

    input:
        file rscript_atcg_qc

    output:
        path "ATCG_SNPs.txt", emit: txt
    //    path "${params.outdir}/QC/ATCG_SNPs.txt"
    //    path("${params.outdir}/QC/ATCG_SNPs.txt"), emit:txt
        path "QC_Info_Plink_ALL.txt"

    script:
    """
    Rscript ${rscript_atcg_qc}
    """
}*/

/*
 * Apply PLINK filters - KEEP markers with P_hwe > 0.00000001, MAF > 0.01, Missingness < 0.03 and exclude ATCG (ambiguous) markers
 */
process PLINK_Filters{
    
    container 'genoimppgs'
    publishDir params.outdir, mode: 'symlink'

    input:
        path genofiles
        //tuple path(genofiles), file("${genofiles}.bed"), file("${genofiles}.bim"), file("${genofiles}.fam")
        //path atcg_txt
        path filtered_genofiles

    output:
        set path("${filtered_genofiles}.bed"),path("${filtered_genofiles}.bim"),path("${filtered_genofiles}.fam") into genofiles_qc
        // emit:"${out}.bed","${out}.bim", "${out}.fam"
        //path "${out}.bed"
        //path "${out}.bim"
        //path "${out}.fam"
        file "${filtered_genofiles}.log"
//plink2 --bfile "${genofiles}" --exclude "${atcg_txt}"" --make-bed --hwe 0.00000001 --maf 0.01 --geno 0.3 --out "${filtered_genofiles}"
    script:
    """
        plink2 --bfile "${genofiles}" --make-bed --hwe 0.00000001 --maf 0.01 --geno 0.3 --out "${filtered_genofiles}"
        
    """
}

/*
 * Convert from PLINK to VCF format
 *
process PLINK_VCF{
    
    container 'genoimppgs'
    publishDir params.outdir, mode: 'symlink'

    input:
        tuple path(genofiles), file("${genofiles}.bed"),file("${genofiles}.bim"),file("${genofiles}.fam")
        path out

    output:
        path "${out}.vcf"
        path "${out}.log"
       //plink2 --bfile "${params.inputdir}/${genofiles}" --recode vcf --out "${params.inputdir}/${out}"
    script:
    """
        plink2 --bfile "${genofiles}" --chr 1 --keep-allele-order --recode vcf bgz --out $out
        
    """
}

/*
 * Impute genotypes using the cloned repository of imputationserver2 (Michigan Imputation Server)
 *
process imputeWithMIS {
    
    publishDir "./output_data/Imputed_Geno/"

    input:
        path config_file_MIS

    script:
    """
        echo nextflow run ./imputationserver2/main.nf -c "${config_file_MIS}"
        
    """
}*/

workflow {

    genofiles = Channel.fromPath(params.input_genoS)
    genofiles.view()
    //geno_qc = "QC"
    //genofiles = Channel.fromFilePairs("${params.input_genoS}.{bed,bim,fam}",size:3)
    //genofiles.view()
    //genofiles.getBaseName()
    //genotypes_name = ${params.input_genoS}
    //genotypes = Channel.fromFilePairs("${params.input_genoS}.{bed,fam,bim}", size:3)
    //genofiles = Channel
    //       .fromFilePairs("${params.input_genoS}.{bed,fam,bim}", size:3)
    //       .ifEmpty {error "No matching plink files"}
    //        .set {plink_data}
    //genofiles.view()
    //genotypes_name = genotypes.getBaseName()
    qc_output_dir = "${params.outdir}/QC/QC"
    //atcg_txt = "${params.outdir}/QC/ATCG_SNPs.txt"
    //postQC_genotypes = "${params.inputdir}/${genotypes}_postQC"
    //postQC_genotypes_name = "${genotypes}_postQC"
    //MIS_config = "${params.inputdir}/MIS_args.config"

    // Generate PLINK QC files
    PLINK_QC(genofiles, qc_output_dir)

    // Generate ATCG SNP list
    //R_QC("${params.scriptDir}/ATCG_QC_Report.r")

    // Apply MAF, HWE and MISS filters using PLINK
    PLINK_Filters(genotypes, R_QC.out.txt, postQC_genotypes)

    // Recode to VCF using PLINK
    //PLINK_VCF(PLINK_Filters.out.filtered_genofiles, "${params.inputdir}/chr1")

    // Impute Genotypes
    //imputeWithMIS(MIS_config)

}