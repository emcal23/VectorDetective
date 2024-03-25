process CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda 'samtools=1.19.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), val(region)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-a -d 20'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools index $bam

    samtools \\
        consensus \\
        $args \\
        -r ${region} \\
        -o ${prefix}_${region}.fasta \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${prefix}_${region}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
