process SNP_call_freebayes {
    tag { sample_id }
    label 'bigmem'
    conda 'envs/freebayes.yml' 
    publishDir "./Results/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_path), path(ref_index_path)

    output:
    path "7.freebayes/${sample_id}/*"

    script:
    """
    mkdir -p 7.freebayes/${sample_id}

    echo "[`date`] Running FreeBayes for ${sample_id}" > 7.freebayes/${sample_id}/${sample_id}_freebayes.log

    freebayes \\
        -f ${ref_index_path}/ref.fa \\
        ${bam_path}/${sample_id}.sorted.bam \\
        > 7.freebayes/${sample_id}/${sample_id}.vcf 2>> 7.freebayes/${sample_id}/${sample_id}_freebayes.log

    echo "[`date`] FreeBayes finished for ${sample_id}" >> 7.freebayes/${sample_id}/${sample_id}_freebayes.log
    """
}
