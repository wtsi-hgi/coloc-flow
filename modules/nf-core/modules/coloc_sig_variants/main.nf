
process COLOC_ON_SIG_VARIANTS {
    tag "${variant}"
    
    label 'process_tiny'
    // publishDir "${params.outdir}/coloc/${GWAS}/${eQTL_path}", mode: "${params.copy_mode}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/bhf_finemap/coloc/coloc.img"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "to be replaced"
    }
    container "/lustre/scratch123/hgi/mdt1/projects/bhf_finemap/coloc/pipeline_ip13/coloc-ip13.sif"
    input:
        each variant
        path(GWAS)
        tuple val(bfile), path(plink_files)

    output:
        path('Done.tmp', emit: done)

    script:
    gwas_name ="${variant}".split("--")[1]
    eQTL_path ="${variant}".split("--")[2]
    variant_name ="${variant}".split("--")[0]
    """
    cp $projectDir/bin/dataIO.R ./
    cp $projectDir/bin/cojo.R ./
    cp $projectDir/bin/helpers.R ./
        coloc_sig_variants.R \
            --gwas ${gwas_name} \
            --rs ${variant_name} \
            --bfile ${bfile} \
            --eqtl ${eQTL_path} \
            --eqtl_snps ${params.eqtl_snps}
        echo ${variant_name} > Done.tmp
    """
}
