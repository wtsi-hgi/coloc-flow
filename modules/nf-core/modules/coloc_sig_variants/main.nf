
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
    container "/lustre/scratch123/hgi/projects/bhf_finemap/coloc/coloc.img"
    input:
        each variant
        path(ped_file_prefix)
        path("${ped_file_prefix}.bed")
        path("${ped_file_prefix}.bim")
        path("${ped_file_prefix}.fam")

    output:
        path('Done.tmp', emit: done)

    script:
    gwas_name ="${variant}".split("--")[1]
    eQTL_path ="${variant}".split("--")[2]
    variant_name ="${variant}".split("--")[0]
    """
        Rscript coloc_sig_variants.R \
            --gwas ${gwas_name} \
            --rs ${variant_name} \
            --bfile ${ped_file_prefix} \
            --eqtl ${eQTL_path} \
            --eqtl_snps ${params.eqtl_snps}
        echo ${variant_name} > Done.tmp
    """
}
