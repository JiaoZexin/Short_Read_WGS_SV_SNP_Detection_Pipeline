process sort_and_index_bam {
    tag { sample_id }

    input:
    tuple val(sample_id), path(bam_file)

    output:
    path "5.sorted/${sample_id}.sorted.bam"
    path "5.sorted/${sample_id}.sorted.bam.bai"
    path "5.sorted/${sample_id}_sort.log"

    script:
    """
    mkdir -p 5.sorted
    start_time=\$(date +%s)
    echo "[`date`] Sorting BAM for ${sample_id}" > 5.sorted/${sample_id}_sort.log
    samtools sort -o 5.sorted/${sample_id}.sorted.bam ${bam_file} 2>> 5.sorted/${sample_id}_sort.log
    samtools index 5.sorted/${sample_id}.sorted.bam 2>> 5.sorted/${sample_id}_sort.log
    end_time=\$(date +%s)
    echo "Sorting + Indexing finished in \$((end_time - start_time)) seconds" >> 5.sorted/${sample_id}_sort.log
    """
}
