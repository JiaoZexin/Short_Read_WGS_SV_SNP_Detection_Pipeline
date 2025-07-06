process SV_call_cnvpytor {
    tag { sample_id }
    label 'standard'
    conda 'envs/cnvpytor.yml'
    publishDir "./Results/", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_path), path(ref_index_path)

    output:
    path "6.cnvpytor/${sample_id}/*"
  
    script:
    """
    mkdir -p 6.cnvpytor/${sample_id}

    echo "[`date`] Running CNVpytor for ${sample_id}" > ${sample_id}_cnvpytor.log

    cnvpytor -root 6.cnvpytor/${sample_id}/${sample_id}.pytor -rd ${bam_path}/${sample_id}.sorted.bam >> 6.cnvpytor/${sample_id}/${sample_id}_cnvpytor.log 2>&1
    cnvpytor -root 6.cnvpytor/${sample_id}/${sample_id}.pytor -his 1000 10000 100000 >> 6.cnvpytor/${sample_id}/${sample_id}_cnvpytor.log 2>&1
    cnvpytor -root 6.cnvpytor/${sample_id}/${sample_id}.pytor -partition 1000 10000 100000 >> 6.cnvpytor/${sample_id}/${sample_id}_cnvpytor.log 2>&1
    cnvpytor -root 6.cnvpytor/${sample_id}/${sample_id}.pytor -call 1000 10000 100000 >> 6.cnvpytor/${sample_id}/${sample_id}_cnvpytor.log 2>&1

    echo "[`date`] CNVpytor finished for ${sample_id}" >> ${sample_id}_cnvpytor.log
    """
}
