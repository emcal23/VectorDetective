process EXTRACT_MAPPED_READS {
    tag "$meta.id"
    label 'process_low'

    conda 'samtools=1.19.2'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), val(region)

    output:
    tuple val(meta), path("*${region}*.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out    = meta.single_end ? "-o ${prefix}_${region}.fastq.gz" : "-s ${prefix}_${region}_S.fastq.gz -1 ${prefix}_${region}_R1.fastq.gz -2 ${prefix}_${region}_R2.fastq.gz"
    """
    samtools index $bam
    samtools view \\
        -h $bam \\
        $region | \\
    samtools \\
        fastq \\
        -F 4 \\
        --threads ${task.cpus} \\
        $args \\
        - \\
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
    touch ${prefix}_${region}.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
