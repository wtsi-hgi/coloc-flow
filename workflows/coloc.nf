/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// // Validate input parameters
// WorkflowColoc.initialise(params, log)

// // TODO nf-core: Add all file path parameters for the pipeline to the list below
// // Check input path parameters to see if they exist
// def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
// for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

// ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
// ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { COLOC_FREQ_AND_SNPS } from '../modules/nf-core/modules/coloc_frq/main'
// include { GWAS_FREQ } from '../modules/nf-core/modules/coloc_frq/main'
include { COLOC_ON_SIG_VARIANTS } from '../modules/nf-core/modules/coloc_sig_variants/main'
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []
workflow COLOC {

    //
    input_channel = Channel.fromPath(params.input_table, followLinks: true, checkIfExists: true)

    // input_table
    input_channel
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
        .map{row->tuple(row.GWAS, row.eQTL)}
        .set{input_gwas_eqtl}

    input_channel
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
        .map{row->row.GWAS}.unique()
        .set{input_gwas}

    input_channel
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
        .map{row->row.eQTL}.unique()
        .set{input_eQTL}

    // Calculate frequencies and extract number of significant GWAS hits for each input GWAS sum stats.
    COLOC_FREQ_AND_SNPS(input_gwas, params.eqtl_snps)

    // Then for each of the GWAS independent SNPs and each of the corresponding eQTLs we generate a new job - we can split this up lated on even more.
    COLOC_FREQ_AND_SNPS.out.sig_signals_eqtls.splitCsv(header: true, sep: '\t', by: 1)
        .map { row -> row.gwas_name2}
        .set { variant_id }

    // Have to run this on each of the eQTL files separately.
    COLOC_ON_SIG_VARIANTS(
        variant_id,
        COLOC_FREQ_AND_SNPS.out.GWAS.collect(),
        Channel.fromFilePairs("${params.bfile}.{bed,bim,fam}", checkIfExists: true, size: 3)
    )
    // variant_id.view()
    // variant_id
    //   .subscribe onNext: {println "variant_id: $it"},
    //   onComplete: {println "variant_id: done"}

    // COLOC_RUN(input_gwas_eqtl)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
