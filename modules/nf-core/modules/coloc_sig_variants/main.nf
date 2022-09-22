
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
        path(sig_signals)
        path(ped_file_folder)
        path(frx_file)
        

    output:
        // path('Done.tmp')
        path('Done.tmp', emit: done)

    script:
    // gwas_name = "${frx_file}".minus(".random").split(/\./)[0]
    gwas_name ="${variant}".split("-")[1]
    variant_name ="${variant}".split("-")[0]
    """
        echo ${variant_name} > Done.tmp
        echo ${gwas_name}
    """
}
