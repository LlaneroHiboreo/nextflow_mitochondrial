// Define version to use
nextflow.enable.dsl = 2


include { GATK_MITOCHONDRIA } from './workflows/gatk_mt'

//
// WORKFLOW: Run main nf-core/raredisease analysis pipeline
//
workflow GATK_MT {
    GATK_MITOCHONDRIA ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    GATK_MT ()
}