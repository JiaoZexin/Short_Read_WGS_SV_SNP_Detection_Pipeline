process align_bwa {
    tag { sample_id }

    input:
    tuple val(sample_id), path(read1), path(read2)
    path ref_fa

    output:
    path "3.sam/${sample_id}.sam"
    path "3.sam/${sample_id}_bwa.log"

    script:
    """
    mkdir -p 3.sam
    start_time=\$(date +%s)
    echo "[`date`] Running BWA MEM for ${sample_id}" > 3.sam/${sample_id}_bwa.log

    bwa mem -t 4 -R '@RG\\tID:${sample_id}\\tPL:illumina\\tDS:resequencing\\tSM:${sample_id}\\tLB:Library_1' \\
        ${ref_fa} ${read1} ${read2} > 3.sam/${sample_id}.sam 2>> 3.sam/${sample_id}_bwa.log

    end_time=\$(date +%s)
    echo "BWA MEM finished in \$((end_time - start_time)) seconds" >> 3.sam/${sample_id}_bwa.log
    """
}
