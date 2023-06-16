process eCAVIAR {
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
        Rscript autoRunCAVIAR.R -e ../huvecEQTL/huvecImputeGenoEQTL/tensorEQTL_All_pval005.tsv.gz -g /scratch/vasccell/cs806/colocalization/cleanGWAS_Summary_Stats/GWAS_Stroke_Eur_Mishra_2022_Nature_hg38.txt -a tenHuvecEQTL_MishraStroke2022_gp5e5

    """
}