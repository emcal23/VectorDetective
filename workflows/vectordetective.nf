/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowVectordetective.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reference FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Mosquito ACE2
// https://www.ncbi.nlm.nih.gov/gene/?term=Mosquito+ACE-2
// https://www.ncbi.nlm.nih.gov/datasets/gene/id/6031419/

Channel.fromPath("$projectDir/ref/Culicidae_ACE2.fna", checkIfExists: true)
    .map{ it -> tuple("ACE2", file(it[0]))}
    .view()
    ch_culicidae_ACE2

// Homo sapiens COI
// https://www.ncbi.nlm.nih.gov/datasets/gene/id/4512/
Channel.fromPath("$projectDir/ref/Hsap_MTCOI.fna", checkIfExists: true)
    .map{ it -> tuple("COI", file(it[0]))}
    .view()
    ch_nonculicidae_COI

ch_genomes = ch_culicidae_ACE2.mix(ch_nonculicidae_COI)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK                           } from '../subworkflows/local/input_check'
include { FASTQ_TRIM_FASTP_FASTQC as FASTQ_QC   } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { FASTQ_ALIGN_BWA as ALIGN_CULICIDAE    } from '../subworkflows/nf-core/fastq_align_bwa/main' 
include { FASTQ_ALIGN_BWA as ALIGN_NONCULICIDAE } from '../subworkflows/nf-core/fastq_align_bwa/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { BWA_INDEX as PREPARE_GENOME           } from '../modules/nf-core/bwa/index/main' 
include { MULTIQC                               } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS           } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_READS                         } from '../modules/local/extract_reads'
include { CREATE_CONSENSUS                      } from '../modules/local/create_consensus'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VECTORDETECTIVE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Short reads QC and trim adapters
    //
    FASTQ_QC (
        ch_shortreads,
        [],
        params.save_trimmed_fail,
        params.save_merged,
        params.skip_fastp,
        params.skip_fastqc
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.fastqc_raw_zip)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.fastqc_trim_zip)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.trim_json)
    ch_versions = ch_versions.mix(FASTQ_QC.out.versions.ifEmpty(null))

    //
    // Module: Prepare input fasta files
    //

    PREPARE_GENOME (
        ch_genomes
    )

    //
    // SUBWORKFLOW: Align reads to Culicidae genome
    //

    ALIGN_CULICIDAE (
        FASTQ_QC.out.reads,
        PREPARE_GENOME.out.bwa_index.filter{it =~ /ACE/},
        true,
        [],
    )
    ch_culicidae_bam = ALIGN_CULICIDAE.out.bam
    ch_culicidae_bai = ALIGN_CULICIDAE.out.bai
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_CULICIDAE.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_CULICIDAE.out.flagstat)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_CULICIDAE.out.idxstats)

    ch_versions = ch_versions.mix(ALIGN_CULICIDAE.out.versions.first())
    
    //
    // SUBWORKFLOW: Extract relevant reads
    //

    EXTRACT_READS(
        ch_culicidae_bam
    )

    //
    // SUBWORKFLOW: Align reads to additional gene(s): COI
    //

    ALIGN_NONCULICIDAE (
        EXTRACT_CULICIDAE_READS.out.unmapped,
        PREPARE_GENOME.out.bwa_index.filter{it =~ /COI/},
        true,
        [],
    )
    ch_nonculicidae_bam = ALIGN_NONCULICIDAE.out.bam
    ch_nonculicidae_bai = ALIGN_NONCULICIDAE.out.bai
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_NONCULICIDAE.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_NONCULICIDAE.out.flagstat)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_NONCULICIDAE.out.idxstats)

    ch_versions = ch_versions.mix(ALIGN_NONCULICIDAE.out.versions.first())

    //
    // MODULE: Create consensus
    //

    CREATE_CONSENSUS (
        ch_culicidae_bam.mix(ch_nonculicidae_bam)
    )

    ch_versions = ch_versions.mix(CREATE_CONSENSUS.out.versions.first())

    //
    // MODULE: Collect and format software versions
    //
    

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowVectordetective.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowVectordetective.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
