process fastqc {
    tag { sample_id }
    label 'lowmem'
    conda 'envs/multiqc.yml'
    publishDir "Results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), path("*_fastqc.html")

    script:
    """
    fastqc -t 4 -o ./ ${reads}
    """
}


