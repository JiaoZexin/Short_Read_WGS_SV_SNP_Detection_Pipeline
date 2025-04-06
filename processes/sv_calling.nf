process call_smoove {
    tag { sample_id }

    input:
    tuple val(sample_id), path(bam), path(bai), path(ref_fa)

    output:
    path "6.smoove/${sample_id}/*"

    script:
    """
    mkdir -p 6.smoove/${sample_id}
    start_time=\$(date +%s)
    echo "[`date`] Running smoove for ${sample_id}" > 6.smoove/${sample_id}/${sample_id}_smoove.log

    smoove call \\
        --outdir 6.smoove/${sample_id}/ \\
        --name ${sample_id} \\
        --fasta ${ref_fa} \\
        -p 1 \\
        --genotype ${bam} 2>> 6.smoove/${sample_id}/${sample_id}_smoove.log

    end_time=\$(date +%s)
    echo "smoove call finished in \$((end_time - start_time)) seconds" >> 6.smoove/${sample_id}/${sample_id}_smoove.log
    """
}
