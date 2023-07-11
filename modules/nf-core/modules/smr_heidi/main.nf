process SMR_HEIDI {
    tag "${eQTL} - ${GWAS}"
    
    label 'process_medium'
    publishDir "${params.outdir}/smr_heidi/${GWAS}/${eQTL_path}", mode: "${params.copy_mode}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.coloc_container}"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "to be replaced"
    }

    input: 
        tuple val(variant_name), path(gwas_name), path(eQTL_path), val(bfile), path(plink_files)

    output:
        path('Done.tmp')
        path('*.jpg') optional true 
        path('*.tsv') optional true 

    script:

    """
        #smr --bfile ${bfile} --gwas-summary /scratch/vasccell/cs806/colocalization/zhu_SMR/GWAS_Large_Artery_Stroke_Eur_Mishra_2022_Nature_hg38.txt_hg38_smrHEIDI.txt --beqtl-summary tenQTLs/HUVEC_ATACseq_Cis_eqtls --out smrResults/tenHUVEC_ATAC-SEQ_QTL_Mishra2022_Stroke_LAS_SMR_HEIDI
        smr --bfile ${bfile} --gwas-summary ${gwas_name} --beqtl-summary ${eQTL_path} --out smrResults/smr
        Rscript auto_process_SMR_HEIDI_Results.R smrResults/smr.smr
    """
}