//
// Alignment with BWA
//

include { BWA_MEM   } from '../../../modules/local/bwa_mem'
include { BWA_INDEX } from '../../../modules/nf-core/bwa/index/main' 

workflow ALIGN {
    take:
    ch_reads  // channel (mandatory): [ val(meta), [ path(reads) ] ]
    ch_ref    // channel (mandatory): [ val(meta2), path(fasta) ]

    main:
    ch_versions = Channel.empty()

    //
    // Prepare genome
    //

    BWA_INDEX (ch_ref )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    //
    // Map reads with BWA
    //

    BWA_MEM ( 
        ch_reads.combine(BWA_INDEX.out.index)
    )

    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    emit:
    bam      = BWA_MEM.out.bam    // channel: [ val(meta), path(bam) ]
    versions = ch_versions        // channel: [ path(versions.yml) ]
}
