process SV_call_delly {
    tag { sample_id }
    label 'highmem'
    conda 'envs/delly.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_path), path(ref_index_path)

    output:
    path "6.delly/${sample_id}/*"

    script:
    """
    mkdir -p 6.delly/${sample_id}

    echo "[`date`] Running Delly for ${sample_id}" > 6.delly/${sample_id}/${sample_id}_delly.log

    echo "[DEBUG] Using delly at: \$(which delly)"
    ldd \$(which delly)

    delly call \\
        -o 6.delly/${sample_id}/${sample_id}.delly.bcf \\
        -g ${ref_index_path}/ref.fa \\
        ${bam_path}/${sample_id}.sorted.bam 2>> 6.delly/${sample_id}/${sample_id}_delly.log

    echo "[`date`] Delly finished for ${sample_id}" >> 6.delly/${sample_id}/${sample_id}_delly.log
    """
}

