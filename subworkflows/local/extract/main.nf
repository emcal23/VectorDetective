//
// Get fasta files for Mosquito/Non mosquito reads
//

include { MAPPED        } from '../../../modules/local/mapped'
include { UNMAPPED      } from '../../../modules/local/unmapped'
include { CONSENSUS     } from '../../../modules/local/create_consensus'
include { SEQKIT_FQ2FA  } from '../../../modules/local/seqkit_fq2fa' 
include { SEQKIT_SAMPLE } from '../../../modules/local/seqkit_sample' 

workflow EXTRACT {
    take:
    ch_ref // channel: [ val(meta), path(fasta) ]
    ch_bam // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    ch_ref
        .map {it -> it[1]}
        .splitFasta( record: [id: true])
        .map{ it -> it.id }
        .set{ch_regions}

    //
    // MODULE: Get mapped reads
    //
    MAPPED(
        ch_bam.combine(ch_regions)
    )
    ch_versions = ch_versions.mix(MAPPED.out.versions.first())

    //
    // MODULE: Extract unmapped reads
    //
    UNMAPPED(
        ch_bam
    )
    ch_versions = ch_versions.mix(UNMAPPED.out.versions.first())

    //
    // MODULE: Create consensus for genes in reference
    //
    CONSENSUS (
        ch_bam.combine(ch_regions)
    )
    ch_versions = ch_versions.mix(CONSENSUS.out.versions.first())

    //
    // MODULE: Subsample Reads
    //
    SEQKIT_SAMPLE (
        UNMAPPED.out.fastq
    )
    ch_versions = ch_versions.mix(SEQKIT_SAMPLE.out.versions.first())    

    //
    // MODULE: Convert unmapped reads to fasta for blast
    //
    SEQKIT_FQ2FA (
        SEQKIT_SAMPLE.out.fastq
    )
    ch_versions = ch_versions.mix(SEQKIT_FQ2FA.out.versions.first())

    emit:
    fasta     = SEQKIT_FQ2FA.out.fasta          // channel: [ val(meta), path(fasta) ]
    consensus = CONSENSUS.out.fasta             // channel: [ val(meta), path(fasta) ]
    versions  = ch_versions                     // channel: [ versions.yml ]
}
