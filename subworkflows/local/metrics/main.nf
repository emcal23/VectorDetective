//
// Run SAMtools stats, flagstat, coverage, and idxstats
//

include { SAMTOOLS_COVERAGE } from '../../../modules/local/samtools_coverage'
include { SAMTOOLS_DEPTH    } from '../../../modules/local/samtools_depth'
include { SAMTOOLS_STATS    } from '../../../modules/local/samtools_stats'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'

workflow METRICS {
    take:
    ch_bam   // channel: [ val(meta), [ bam ] ]

    main:

    ch_versions = Channel.empty()

    SAMTOOLS_INDEX ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0])
        .set { ch_bam_bai }

    SAMTOOLS_FLAGSTAT ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    SAMTOOLS_IDXSTATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.first())

    SAMTOOLS_STATS ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())

    SAMTOOLS_DEPTH ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    SAMTOOLS_COVERAGE ( ch_bam_bai )
    ch_versions = ch_versions.mix(SAMTOOLS_COVERAGE.out.versions.first())

    emit:
    bai      = SAMTOOLS_INDEX.out.bai         // channel: [ val(meta), [ bai ] ]

    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    coverage = SAMTOOLS_COVERAGE.out.coverage // channel: [ val(meta), [ coverage ] ]
    depth    = SAMTOOLS_DEPTH.out.depth       // channel: [ val(meta), [ depth ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
