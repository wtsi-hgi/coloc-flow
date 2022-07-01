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

include { COLOC_RUN } from '../modules/nf-core/modules/coloc/main' 
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

    COLOC_RUN(input_gwas_eqtl)
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
