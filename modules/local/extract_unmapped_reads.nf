process EXTRACT_UNMAPPED_READS {
    tag "$meta.id"
    label 'process_low'

    conda 'samtools=1.19.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*unmapped*.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out    = meta.single_end ? "-o ${prefix}_unmapped.fastq.gz" : "-s ${prefix}_unmapped_S.fastq.gz -1 ${prefix}_unmapped_R1.fastq.gz -2 ${prefix}_unmapped_R2.fastq.gz"
    """
    samtools \\
        fastq \\
        -f 4 \\
        --threads ${task.cpus} \\
        $args \\
        $bam \\
        $out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unmapped.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
