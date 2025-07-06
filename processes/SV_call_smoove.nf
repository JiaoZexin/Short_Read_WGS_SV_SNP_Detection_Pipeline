process SV_call_smoove {
    tag { sample_id }
    label 'highmem'
    conda 'envs/smoove.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_path), path (ref_index_path)

    output:
    path "6.smoove/${sample_id}/*"

    script:
    """
    mkdir -p 6.smoove/${sample_id}

    
    echo "[`date`] Running smoove for ${sample_id}" > 6.smoove/${sample_id}/${sample_id}_smoove.log

    smoove call \\
        --outdir 6.smoove/${sample_id}/ \\
        --name ${sample_id} \\
        --fasta ${ref_index_path}/ref.fa \\
        -p 1 \\
        ${bam_path}/${sample_id}.sorted.bam 2>> 6.smoove/${sample_id}/${sample_id}_smoove.log

    echo "[`date`] smoove finished for ${sample_id}" >> 6.smoove/${sample_id}/${sample_id}_smoove.log
    """
}
