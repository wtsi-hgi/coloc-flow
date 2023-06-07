
process COLOC_FREQ_AND_SNPS {
    tag "${GWAS} frq"
    
    label 'process_low'
    // publishDir "${params.outdir}/coloc/${GWAS}/${eQTL_path}", mode: "${params.copy_mode}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.coloc_container}"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "to be replaced"
    }
    container "/lustre/scratch123/hgi/projects/bhf_finemap/coloc/coloc.img"
    input: 
        path(GWAS)
        path(bfile)

    output:
        // path('Done.tmp')
        path('*.sig_signals.list', emit: sig_signals)
        path('*.frqx', emit: frqx)
        path("Filtered_${gwas_name}", emit: bed_file)
        path('*all_signals.tsv',emit: sig_signals_eqtls)
        path(GWAS,emit: GWAS)


    script:
    
    gwas_name = "${GWAS}".minus(".random").split(/\./)[0]
            // outfile = "${file__anndata}".minus(".h5ad")
            // .split("-").drop(1).join("-")
    """
        
       
        coloc_GWAS_frq.R ${GWAS} ${bfile} ${params.input_table}
        mkdir Filtered_${gwas_name}
        mv ${gwas_name}.bed Filtered_${gwas_name}/Filtered_${gwas_name}.bed
        mv ${gwas_name}.bim Filtered_${gwas_name}/Filtered_${gwas_name}.bim
        mv ${gwas_name}.fam Filtered_${gwas_name}/Filtered_${gwas_name}.fam
        echo ${gwas_name} > Done.tmp
    """
}


