process SEQKIT_SAMPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "seqkit=2.8.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.1--h9ee0642_0':
        'biocontainers/seqkit:2.8.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*sub_R*.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-n 5000 -s23'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    seqkit \\
        sample \\
        ${args} \\
        ${fastq[0]} \\
        -o ${prefix}_sub_R1.fastq.gz

    seqkit \\
        sample \\
        ${args} \\
        ${fastq[1]} \\
        -o ${prefix}_sub_R2.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_sub_R1.fa.gz
    touch ${prefix}_sub_R2.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}
