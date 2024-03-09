process EXTRACT_READS {
    tag "$meta.id"
    label 'process_low'

    conda 'samtools=1.19.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*mapped*.fastq.gz"), emit: mapped_fastq
    tuple val(meta), path("*unmapped*.fastq.gz"), emit: unmapped_fastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: '-f 4'
    def argz   = task.ext.args ?: '-F 4'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mapped = meta.single_end ? "-o ${prefix}_mapped.fastq.gz"   : "-s ${prefix}_mapped_S.fastq.gz   -1 ${prefix}_mapped_R1.fastq.gz   -2 ${prefix}_mapped_R2.fastq.gz"
    def unmppd = meta.single_end ? "-o ${prefix}_unmapped.fastq.gz" : "-s ${prefix}_unmapped_S.fastq.gz -1 ${prefix}_unmapped_R1.fastq.gz -2 ${prefix}_unmapped_R2.fastq.gz"
    """
    samtools \\
        fastq \\
        -F 4 \\
        --threads ${task.cpus} \\
        $args \\
        $bam \\
        $mapped

    samtools \\
        fastq \\
        -f 4 \\
        --threads ${task.cpus} \\
        $argz \\
        $bam \\
        $unmppd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extractreads: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mapped.fastq.gz
    touch ${prefix}_unmapped.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        extract_reads: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
