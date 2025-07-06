process alignment {
    tag { sample_id }
    label 'bigmem'
    conda 'envs/samtools.yml'
    publishDir "./Results/", mode: 'copy'
    input:
    tuple val(sample_id), path(read1), path(read2), path (ref_index_path)

    output:
    path "3.sam/${sample_id}/*"
    path "4.bam/${sample_id}/*"
    path "5.sorted/${sample_id}/*"
    path "${sample_id}_alignment_full.log"
    tuple val(sample_id), path("5.sorted/${sample_id}/"), emit: bam_path


    script:
    """
    mkdir -p 3.sam/${sample_id} 4.bam/${sample_id} 5.sorted/${sample_id}
    start_time=\$(date +%s)
    echo "[`date`] Running BWA MEM for ${sample_id}" > 3.sam/${sample_id}/${sample_id}_bwa.log

    echo "[\$(date)] Step 1: BWA MEM" >> ${sample_id}_alignment_full.log
    bwa mem -t 8 -R '@RG\\tID:${sample_id}\\tPL:illumina\\tDS:resequencing\\tSM:${sample_id}\\tLB:Library_1' \\
        ${ref_index_path}/ref.fa ${read1} ${read2} > 3.sam/${sample_id}/${sample_id}.sam 2>> ${sample_id}_alignment_full.log || exit 1

    echo "[\$(date)] Step 2: Convert SAM to BAM" >> ${sample_id}_alignment_full.log
    samtools view -Sb -o 4.bam/${sample_id}/${sample_id}.bam 3.sam/${sample_id}/${sample_id}.sam 2>> ${sample_id}_alignment_full.log || exit 1

    echo "[\$(date)] Step 3: Sort BAM" >> ${sample_id}_alignment_full.log
    samtools sort -o 5.sorted/${sample_id}/${sample_id}.sorted.bam 4.bam/${sample_id}/${sample_id}.bam 2>> ${sample_id}_alignment_full.log || exit 1

    echo "[\$(date)] Step 4: Index BAM" >> ${sample_id}_alignment_full.log
    samtools index 5.sorted/${sample_id}/${sample_id}.sorted.bam 2>> ${sample_id}_alignment_full.log || exit 1

    end_time=\$(date +%s)
    echo "[\$(date)] Alignment pipeline finished in \$((end_time - start_time)) seconds" >> ${sample_id}_alignment_full.log
    echo "[\$(date)] Alignment for ${sample_id} completed successfully." >> ${sample_id}_alignment_full.log

    """

}
