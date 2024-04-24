//
// Get fasta files for Mosquito/Non mosquito reads
//

include { BLOBTOOLS              } from '../../../modules/local/blobtools'
include { BLAST_BLASTN           } from '../../../modules/nf-core/blast/blastn/main'

workflow BLOBTOOLS {
    take:
    ch_reads     // channel: [ val(meta), [ fastq ] ]
    ch_consensus // channel: [ val(meta), [ fasta ] ]
    ch_blastdb   // channel: [ val(meta2), [ db ] ]

    main:

    ch_versions = Channel.empty()

    ch_consensus
        .grouptuple()
        .set { ch_fasta }

    ALIGN (
        ch_reads,
        ch_consensus
    )

    EXTRACT_MAPPED_READS(
        ch_bam.combine(ch_regions)
    )
    ch_versions = ch_versions.mix(EXTRACT_MAPPED_READS.out.versions.first())

    EXTRACT_UNMAPPED_READS(
        ch_bam
    )
    ch_versions = ch_versions.mix(EXTRACT_UNMAPPED_READS.out.versions.first())

    //
    // MODULE: Create consensus for genes in reference
    //

    CONSENSUS (
        ch_bam.combine(ch_regions)
    )

    ch_versions = ch_versions.mix(CONSENSUS.out.versions.first())

    //
    // MODULE: De novo assembly of unmapped reads
    //
    SPADES (
        EXTRACT_UNMAPPED_READS.out.fastq,
        [],
        []
    )

    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}
