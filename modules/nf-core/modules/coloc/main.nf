
process COLOC_RUN {
    tag "${eQTL} - ${GWAS}"
    
    label 'process_medium'
    publishDir "${params.outdir}/coloc/${GWAS}/${eQTL_path}", mode: "${params.copy_mode}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.coloc_container}"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "to be replaced"
    }

    input: 
        tuple path(GWAS), path(eQTL)

    output:
        // path('Done.tmp')
        path('*.jpg') optional true 
        path('*.tsv') optional true 

    script:
    eQTL_path = "${eQTL}".minus(".gz")
            .split("\\.\\d")[0]
    """
        echo ${eQTL_path} > eQTL_path.txt
        coloc_GWAS.R ${eQTL} ${GWAS} ${params.bfile}
        echo Done > Done.tmp
    """
}
