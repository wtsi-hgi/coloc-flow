
process COLOC_FREQ_AND_SNPS {
    tag "${GWAS} frq"
    
    label 'process_low'
    // publishDir "${params.outdir}/coloc/${GWAS}/${eQTL_path}", mode: "${params.copy_mode}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/bhf_finemap/coloc/coloc.img"
    } else {
        container "to be replaced"
    }
    container "/lustre/scratch123/hgi/projects/bhf_finemap/coloc/pipeline_ip13/coloc-ip13.sif"
    input: 
        path(GWAS)
        path(eqtl_snps)

    output:
        // path('Done.tmp')
        path('*.sig_signals.list', emit: sig_signals)
        path('*all_signals.tsv', emit: sig_signals_eqtls)
        path(GWAS, emit: GWAS)


    script:
    
    gwas_name = "${GWAS}".minus(".random").split(/\./)[0]
            // outfile = "${file__anndata}".minus(".h5ad")
            // .split("-").drop(1).join("-")
    """
        cp $projectDir/bin/dataIO.R ./
        cp $projectDir/bin/helpers.R ./
        coloc_GWAS_frq.R \
            --gwas ${GWAS} \
            --input ${params.input_table} \
            --eqtl_snps ${eqtl_snps}
        echo ${gwas_name} > Done.tmp
    """
}

process GWAS_FREQ {
    cpus 1
    memory '6 GB'
    container "/lustre/scratch123/hgi/projects/bhf_finemap/coloc/coloc.img"
    input:
        path(bfile)
    output:
        path(${outfile}, emit: bfile)
        path("*.frq", emit: frq)
    script:
    outfile = ${bfile.name}.filtered
    """
        plink --bfile ${bfile} --maf 0.0001 --make-bed --freq --out ${outfile}
    """
}
