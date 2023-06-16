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
        tuple path(GWAS), path(eQTL)

    output:
        path('Done.tmp')
        path('*.jpg') optional true 
        path('*.tsv') optional true 

    script:

    """
        smr --bfile /scratch/vasccell/cs806/colocalization/1000Genome/euroSamps1kGMerge --gwas-summary /scratch/vasccell/cs806/colocalization/zhu_SMR/GWAS_Large_Artery_Stroke_Eur_Mishra_2022_Nature_hg38.txt_hg38_smrHEIDI.txt --beqtl-summary tenQTLs/HUVEC_ATACseq_Cis_eqtls --out smrResults/tenHUVEC_ATAC-SEQ_QTL_Mishra2022_Stroke_LAS_SMR_HEIDI
        Rscript auto_process_SMR_HEIDI_Results.R smrResults/tenHUVEC_ATAC-SEQ_QTL_Mishra2022_Stroke_LAS_SMR_HEIDI.smr
    """
}