
process COLOC_ON_SIG_VARIANTS {
    tag "${variant}"
    
    label 'process_medium'
    // publishDir "${params.outdir}/coloc/${GWAS}/${eQTL_path}", mode: "${params.copy_mode}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.coloc_container}"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "to be replaced"
    }
    container "/lustre/scratch123/hgi/mdt1/projects/bhf_finemap/coloc/pipeline_ip13/coloc-ip13.sif"
    memory '15 GB'

    input:
        tuple val(variant_name), path(gwas_name), path(eQTL_path), val(bfile), path(plink_files)

    output:
        path('Done.tmp', emit: done)

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
