process trimmomatic {
    tag { sample_id }
    label 'standard'
    conda 'envs/trimmomatic.yml'
    publishDir "./Results/1.trim", mode: 'copy'

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R1.trimmed.fastq.gz"), path("${sample_id}_R2.trimmed.fastq.gz")

    script:
    """
    trimmomatic PE -threads 4 \
      ${read1} ${read2} \
      ${sample_id}_R1.trimmed.fastq.gz /dev/null \
      ${sample_id}_R2.trimmed.fastq.gz /dev/null \
      ILLUMINACLIP:${params.adapters}:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
