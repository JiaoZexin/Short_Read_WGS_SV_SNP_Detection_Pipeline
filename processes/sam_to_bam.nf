process sam_to_bam {
    tag { sample_id }

    input:
    tuple val(sample_id), path(sam_file)

    output:
    path "4.bam/${sample_id}.bam"
    path "4.bam/${sample_id}_samtools.log"

    script:
    """
    mkdir -p 4.bam
    start_time=\$(date +%s)
    echo "[`date`] Converting SAM to BAM for ${sample_id}" > 4.bam/${sample_id}_samtools.log
    samtools view -Sb -o 4.bam/${sample_id}.bam ${sam_file} 2>> 4.bam/${sample_id}_samtools.log
    end_time=\$(date +%s)
    echo "SAM to BAM finished in \$((end_time - start_time)) seconds" >> 4.bam/${sample_id}_samtools.log
    """
}
