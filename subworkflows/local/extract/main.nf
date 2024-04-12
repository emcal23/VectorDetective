//
// Get fasta files for Mosquito/Non mosquito reads
//

include { EXTRACT_MAPPED_READS   } from '../../../modules/local/extract_mapped_reads'
include { EXTRACT_UNMAPPED_READS } from '../../../modules/local/extract_unmapped_reads'
include { CONSENSUS              } from '../../../modules/local/create_consensus'
include { SPADES                 } from '../../../modules/nf-core/spades/main' 

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


    // //
    // // MODULE: Identify contigs
    // //
    // BLAST_BLASTN ( 
    //     SPADES.out.contigs,
    //     ch_blastdb
    // )

    emit:

    versions = ch_versions                     // channel: [ versions.yml ]
}
