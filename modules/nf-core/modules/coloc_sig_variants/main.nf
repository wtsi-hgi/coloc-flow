
process COLOC_ON_SIG_VARIANTS {
    tag "${variant_name}"
    
    label 'process_medium'
    publishDir "${params.outdir}/coloc/${gwas_name}/${eQTL_path}", mode: "${params.copy_mode}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.coloc_container}"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "to be replaced"
    }

    input:
        tuple val(variant_name), path(gwas_name), path(eQTL_path), val(bfile), path(plink_files)

    output:
        path 'coloc_results.csv', emit: done optional true 
        path('*.jpg') optional true 
        path('*.pdf') optional true 
    script:
    """
    cp $projectDir/bin/dataIO.R ./dataIO.R
    cp $projectDir/bin/cojo.R ./cojo.R
    cp $projectDir/bin/helpers.R ./helpers.R
        coloc_sig_variants.R \
            --gwas ${gwas_name} \
            --rs ${variant_name} \
            --bfile ${bfile} \
            --eqtl ${eQTL_path} \
            --eqtl_snps ${params.eqtl_snps} \
            --config ${params.yaml}
        echo ${variant_name} > Done.tmp
    """
}
